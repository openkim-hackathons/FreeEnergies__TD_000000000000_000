import re
from enum import Enum
from pathlib import Path
from typing import Optional, Sequence

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sc
from ase import Atoms
from ase.data import atomic_masses, atomic_numbers
from ase.io import read
from kim_tools import (
    get_isolated_energy_per_atom,
    get_stoich_reduced_list_from_prototype,
)
from kim_tools.symmetry_util.core import reduce_and_avg
from kim_tools.test_driver import SingleCrystalTestDriver
from kim_tools.kimunits import convert_units
from scipy import integrate

from .helper_functions import run_lammps
from .lammps_template import LammpsTemplate
from .structure_utils import compute_supercell_for_target_size

# Physical constants from scipy
EV = sc.value("electron volt")
MU = sc.value("atomic mass constant")
HBAR = sc.value("Planck constant in eV/Hz") / (2 * np.pi)
KB = sc.value("Boltzmann constant in eV/K")


class LammpsStatus(Enum):
    """Status of LAMMPS simulation completion."""
    SUCCESS = "success"
    MELTED = "melted"
    NOT_FOUND = "not_found"
    INCOMPLETE = "incomplete"


class TestDriver(SingleCrystalTestDriver):
    """Test driver for computing Gibbs free energy using Frenkel-Ladd method."""

    # Unit conversion constants
    BAR_TO_PA = 1e5
    ANGSTROM3_TO_M3 = 1e-30
    JOULE_TO_EV = convert_units(
        from_value=1, from_unit="J", wanted_unit="eV", suppress_unit=True
        )

    # Accuracy thresholds
    DEFAULT_RELATIVE_ACCURACY = 0.01
    ORTHOGONAL_THRESHOLD_DEGREES = 0.1

    def _calculate(
        self,
        timestep_ps: float = 0.001,
        fl_switch_timesteps: int = 50000,
        fl_equil_timesteps: int = 10000,
        target_size: int = 200,
        repeat: Optional[Sequence[int]] = None,
        lammps_command: str = "lmp",
        msd_threshold: float = 0.1,
        msd_timesteps: int = 20000,
        thermo_sampling_period: int = 100,
        ave_pos_timesteps: int = 50000,
        random_seed: int = 101010,
        rlc_n_every: int = 10,
        rlc_initial_run_length: int = 10000,
        rlc_min_samples: int = 5,
        output_dir: str = "output",
        equilibration_plots: bool = True,
        fl_plots: bool = True,
        k_factor: float = 1.0,
        **kwargs) -> None:
        """Compute Gibbs free energy at constant temperature and hydrostatic pressure using a nonequilibrium
        thermodynamic integration method implemented in LAMMPS.

        This test driver repeats the unit cell of the zero-temperature crystal structure to build a supercell and then
        runs molecular-dynamics simulations in the NVT and NPT ensembles using LAMMPS.

        This test driver uses kim_convergence to detect equilibrated molecular-dynamics simulations. It checks for the
        convergence of the volume, temperature, and cell shape parameters at a set frequency.

        After the molecular-dynamics NVT simulation, the symmetry of the equilibrium time-averaged structure is checked
        to ensure that it did not change in comparison to the initial structure. Also, it is ensured that replicated atoms
        in replicated unit atoms are not too far away from the average atomic positions.

        The crystal might melt or vaporize during the simulations. In that case, kim-convergence would only detect
        equilibration after unnecessarily long simulations. Therefore, this test driver initially checks for melting or
        vaporization during short initial simulations. During these initial runs, the mean-squared displacement (MSD) of
        atoms during the simulations is monitored. If the average time-derivative of MSD exceeds a given threshold value
        (msd_threshold), an error is raised.
        
        Args:
            timestep_ps: MD timestep in picoseconds.
            fl_switch_timesteps: Number of timesteps for Frenkel-Ladd switching.
            fl_equil_timesteps: Number of timesteps for equilibration before switching.
            target_size: Target number of atoms for supercell. Uses cutoff-based expansion
                with target size constraint (good for non-cubic cells). The algorithm starts
                with an 8Å cutoff radius and recursively decreases it until the supercell has
                fewer atoms than target_size. Default: 10000. Ignored if repeat is specified.
            repeat: Explicit (nx, ny, nz) supercell repetitions. If specified, this takes
                precedence over target_size and the supercell is created with these exact
                repetitions. Default: None (uses target_size instead).
            lammps_command: Command to invoke LAMMPS executable.
            msd_threshold: MSD threshold for detecting melting/vaporization. Has units of angstrom^2 per
                "thermo_sampling_period" timesteps.
            thermo_sampling_period: Number of timesteps between thermodynamic data measurements.
            msd_timesteps: Number of timesteps to monitor MSD for melting detection.
            ave_pos_timesteps: Number of timesteps for averaging atomic positions after run_length_control/kim-convergence ends.
            random_seed: Random seed for velocity initialization and Langevin thermostat.
            rlc_n_every: Sampling frequency for run length control.
            rlc_initial_run_length: Initial run length for convergence checking.
            rlc_min_samples: Minimum number of independent samples for convergence.
            output_dir: Directory for output files.
            equilibration_plots: Whether to generate equilibration diagnostic plots.
            fl_plots: Whether to generate Frenkel-Ladd switching plots.
            k_factor: Multiply calculated spring constants by this value before running Frenkel-Ladd switching.
                Values greater than unity alleviate lost-atoms error caused by unrealistic coordination of atoms
                during the switching procedure.
            **kwargs: Additional keyword arguments (unused).
            
        Note:
            Supercell sizing: If repeat is specified, uses those exact repetitions. Otherwise,
                uses cutoff-based expansion with target size constraint. This method ensures
                the supercell stays below the target atom count while maintaining appropriate
                minimum image distances, making it suitable for both cubic and non-cubic
                crystal structures.
            
        Raises:
            ValueError: If input parameters are invalid.
            RuntimeError: If simulation fails (melting, incomplete, etc.).
            FileNotFoundError: If LAMMPS log file is not found.
        """
        # Initialize parameters and validate
        self._initialize_params(
            timestep_ps, fl_switch_timesteps, fl_equil_timesteps, target_size, repeat,
            lammps_command, msd_threshold, msd_timesteps, thermo_sampling_period, ave_pos_timesteps,
            random_seed, rlc_n_every, rlc_initial_run_length, rlc_min_samples, output_dir, k_factor
        )
        
        # Run LAMMPS simulation
        self._run_simulation()
        
        # Compute free energy
        free_energy_per_atom = self._free_energy()

        # Make diagnostic plots
        if equilibration_plots:
            self._plot_equilibration()
        if fl_plots:
            self._plot_frenkel_ladd()

        # Process and report results
        self._finalize_and_report_results(free_energy_per_atom)

    def _initialize_params(
        self,
        timestep_ps: float,
        fl_switch_timesteps: int,
        fl_equil_timesteps: int,
        target_size: int,
        repeat: Optional[Sequence[int]],
        lammps_command: str,
        msd_threshold: int,
        msd_timesteps: int,
        thermo_sampling_period: int,
        ave_pos_timesteps: int,
        random_seed: int,
        rlc_n_every: int,
        rlc_initial_run_length: int,
        rlc_min_samples: int,
        output_dir: str,
        k_factor: float
    ) -> None:
        """Initialize calculation parameters and instance variables."""
        # Get crystal structure info
        self.prototype_label = self._get_nominal_crystal_structure_npt()["prototype-label"]["source-value"]
        self.temperature_K = self._get_temperature(unit="K")
        self.cauchy_stress = self._get_cell_cauchy_stress(unit='bars')
        self.pressure = -self.cauchy_stress[0]
        self.atoms = self._get_atoms()
        
        # Store simulation parameters
        self.timestep_ps = timestep_ps
        self.fl_switch_timesteps = fl_switch_timesteps
        self.fl_equil_timesteps = fl_equil_timesteps
        self.target_size = target_size
        self.repeat = repeat
        # If repeat is None, it will be set by compute_supercell_for_target_size
        self.lammps_command = lammps_command
        self.msd_threshold = msd_threshold
        self.msd_timesteps = msd_timesteps
        self.thermo_sampling_period = thermo_sampling_period
        self.ave_pos_timesteps = ave_pos_timesteps
        self.random_seed = random_seed
        self.rlc_n_every = rlc_n_every
        self.rlc_initial_run_length = rlc_initial_run_length
        self.rlc_min_samples = rlc_min_samples
        self.output_dir = output_dir
        self.k_factor = k_factor
        
        self.atom_style = self._get_supported_lammps_atom_style()
        
        self._validate_inputs()

    def _run_simulation(self) -> None:
        """Set up and run the LAMMPS simulation."""
        # Set up initial structure
        self.supercell = self._setup_initial_structure(
            filename=f"{self.output_dir}/zero_temperature_crystal.data"
        )
        
        # Prepare LAMMPS input files
        self._modify_run_length_control()
        self.templates = LammpsTemplate(root=f"{self.output_dir}/")
        self.templates._write_lammps_file(
            nspecies=len(self.species), is_triclinic=self.is_triclinic
        )
        
        # Get species list for LAMMPS
        symbols = self.atoms.get_chemical_symbols()
        species = sorted(set(symbols))
        
        # Run LAMMPS
        _, _, self.spring_constants, self.volume, self.lammps_status = run_lammps(
            self.kim_model_name, self.temperature_K, self.pressure,
            self.timestep_ps, self.fl_switch_timesteps, self.fl_equil_timesteps,
            species, self.msd_threshold, self.msd_timesteps, self.thermo_sampling_period,
            self.ave_pos_timesteps, self.rlc_n_every, self.lammps_command, self.random_seed,
            self.output_dir, self.k_factor
        )
        
        # Verify simulation completed successfully
        self._verify_lammps_completion()

    def _verify_lammps_completion(self) -> None:
        """Verify that LAMMPS simulation completed successfully."""
        if self.lammps_status == LammpsStatus.MELTED:
            raise RuntimeError("Crystal melted or vaporized")
        elif self.lammps_status == LammpsStatus.NOT_FOUND:
            raise FileNotFoundError(
                f"LAMMPS log file not found: {self.output_dir}/free_energy.log"
            )
        elif self.lammps_status == LammpsStatus.INCOMPLETE:
            raise RuntimeError("LAMMPS simulation did not complete successfully")

    def _finalize_and_report_results(self, free_energy_per_atom: float) -> None:
        """Process results, compute derived quantities, and write properties.
        
        Args:
            free_energy_per_atom: Raw free energy per atom from Frenkel-Ladd integration.
        """
        # Reduce supercell and update crystal structure
        reduced_atoms = self._reduce_and_average_supercell(
            atoms_npt_path=f"{self.output_dir}/free_energy.data",
            reduced_atoms_save_path=f"{self.output_dir}/reduced_atoms.data"
        )
        self._update_nominal_parameter_values(reduced_atoms)
        
        # Write crystal structure property
        self._add_property_instance_and_common_crystal_genome_keys(
            "crystal-structure-npt", write_temp=True, write_stress=True
        )
        self._add_file_to_current_property_instance(
            "restart-file", f"{self.output_dir}/free_energy.restart"
        )

        # Subtract isolated atom energies
        isolated_atom_energy = self._collect_isolated_atom_energies(reduced_atoms)
        free_energy_per_atom = free_energy_per_atom - isolated_atom_energy

        # Compute derived quantities
        num_atoms_in_formula = sum(self.stoichiometry)
        free_energy_per_formula = free_energy_per_atom * num_atoms_in_formula

        atoms_per_cell = len(reduced_atoms.get_masses())
        mass_per_cell = sum(reduced_atoms.get_masses())
        free_energy_per_cell = free_energy_per_atom * atoms_per_cell
        specific_free_energy = free_energy_per_cell / mass_per_cell

        # Print results
        print("####################################")
        print("# Frenkel Ladd Free Energy Results #")
        print("####################################")
        print(f"G_FL = {free_energy_per_atom:.5f} (eV/atom)")
        print(f"Isolated atom energy = {isolated_atom_energy:.5f} (eV/atom)")

        # Write free energy properties
        self._add_property_instance_and_common_crystal_genome_keys(
            "gibbs-free-energy-crystal", write_stress=True, write_temp=self.temperature_K
        )
        self._add_key_to_current_property_instance(
            "gibbs-free-energy-per-atom", free_energy_per_atom, "eV"
        )
        self._add_key_to_current_property_instance(
            "gibbs-free-energy-per-formula", free_energy_per_formula, "eV"
        )
        self._add_key_to_current_property_instance(
            "specific-gibbs-free-energy", specific_free_energy, "eV/amu"
        )

    def _validate_inputs(self) -> None:
        """Validate all input parameters for the calculation.
        
        Checks that all simulation parameters have valid values and are
        physically meaningful. Ensures mutual exclusivity of supercell
        sizing options.
        
        Raises:
            ValueError: If any parameter has an invalid value, including:
                - Non-positive temperature, timestep, or random seed
                - Invalid stress tensor (must be hydrostatic pressure)
                - MSD parameters inconsistent with sampling frequency
        """
        if not self.temperature_K > 0.0:
            raise ValueError("Temperature must be greater than zero.")

        if not len(self.cauchy_stress) == 6:
            raise ValueError("Specify all six (x, y, z, xy, xz, yz) entries of the cauchy stress tensor.")

        if not (self.cauchy_stress[0] == self.cauchy_stress[1] == self.cauchy_stress[2]):
            raise ValueError("The diagonal entries of the stress tensor must be equal so that a hydrostatic "
                             "pressure is used.")

        if not (self.cauchy_stress[3] == self.cauchy_stress[4] == self.cauchy_stress[5] == 0.0):
            raise ValueError("The off-diagonal entries of the stress tensor must be zero so that a hydrostatic "
                             "pressure is used.")

        if not self.timestep_ps > 0.0:
            raise ValueError("Timestep must be greater than zero.")

        if not self.fl_switch_timesteps > 0:
            raise ValueError("The number of timesteps for Frenkel-Ladd switching must be greater than 0.")

        if not self.fl_equil_timesteps > 0:
            raise ValueError("The number of timesteps for equilibration before Frenkel-Ladd switching must be greater than 0.")

        # Validate repeat if provided
        if self.repeat is not None:
            if not len(self.repeat) == 3:
                raise ValueError("The repeat argument must be a tuple of three integers.")
            if not all(r > 0 for r in self.repeat):
                raise ValueError("All number of repeats must be greater than zero.")
        
        if not self.msd_threshold > 0.0:
            raise ValueError("The mean-squared displacement slope threshold must be greater than zero.")

        if not self.msd_timesteps > 0:
            raise ValueError("The number of timesteps to monitor the mean-squared displacement must be greater than zero.")
        
        if not self.thermo_sampling_period > 0:
            raise ValueError("The number of timesteps between thermo samples must be greater than zero.")

        if not self.msd_timesteps % self.thermo_sampling_period == 0:
            raise ValueError("The number of timesteps over which to monitor the mean-squared displacement must be "
                             "an integer multiple of the number of timesteps between mean-squared displacement samples.")
        
        if not self.ave_pos_timesteps > 0:
            raise ValueError("The number of timesteps for averaging atomic positions must be greater than zero.")
        
        if not self.rlc_n_every > 0:
            raise ValueError("The number of timesteps between rlc samples must be greater than zero.")
        
        if not self.rlc_initial_run_length > 0:
            raise ValueError("The number of timesteps between run_length_control equilbration checks must be greater than zero.")
        
        if not self.rlc_min_samples > 0:
            raise ValueError("The minimum number of uncorrelated samples required to determine equilibration must be greater than zero.")

        if not self.random_seed > 0:
            raise ValueError("The random seed must be greater than zero.")
        
        if not self.k_factor > 0:
            raise ValueError("The k_factor must be greater than zero.")

    def _setup_initial_structure(self, filename: str) -> Atoms:
        """Set up the initial supercell structure for the simulation.
        
        Creates a supercell from the primitive cell. If repeat is specified, uses
        those exact repetitions. Otherwise, uses cutoff-based expansion with target
        size constraint. Determines species, masses, and concentrations, then
        writes the structure to a LAMMPS data file.
        
        Args:
            filename: Path to write the LAMMPS data file.
            
        Returns:
            ASE Atoms object representing the supercell.
            
        Note:
            Sets instance variables: species, mass, concentration, zero_k_structure_path,
            is_triclinic, repeat.
        """
        # Copy original atoms so that their information does not get lost when the new atoms are modified.
        atoms_new = self.atoms.copy()

        # Build supercell
        if self.repeat is not None:
            # Use explicit repeat tuple
            atoms_new = atoms_new.repeat(self.repeat)
        else:
            # Use cutoff-based expansion with target size constraint
            # (good for non-cubic cells, ensures natoms >= target_size)
            atoms_new, self.repeat = compute_supercell_for_target_size(
                self.atoms, self.target_size
            )

        # This is how ASE obtains the species that are written to the initial configuration.
        # These species are passed to kim interactions.
        # See https://wiki.fysik.dtu.dk/ase/_modules/ase/io/lammpsdata.html#write_lammps_data
        symbols = atoms_new.get_chemical_symbols()
        self.species = np.array(sorted(set(symbols)))
        self.mass = np.array(
            [
                atomic_masses[atomic_numbers[element_symbol]]
                for element_symbol in self.species
            ]
        )

        symbols = self.atoms.get_chemical_symbols()
        self.concentration = np.array([
            symbols.count(s) / len(symbols) for s in self.species
        ])

        atoms_new.write(filename, format="lammps-data", masses=True, atom_style=self.atom_style)
        self.zero_k_structure_path = filename

        # Treat all cells as triclinic for simplicity
        self.is_triclinic = True

        return atoms_new
    
    def _free_energy(self) -> float:
        """Compute the Gibbs free energy per atom using Frenkel-Ladd integration.
        
        Integrates the work done during forward and backward switching between
        the interatomic potential and harmonic reference system. Includes corrections
        for harmonic free energy, center of mass constraint, and PV term.
        
        Returns:
            Gibbs free energy per atom in eV (before subtracting isolated atom energies).
        """
        Hi_f, Hf_f, lamb_f = np.loadtxt(
            f"{self.output_dir}/FL_switch1.dat", unpack=True, skiprows=1
        )
        w_forw = integrate.simpson(Hf_f - Hi_f, x=lamb_f)

        Hf_b, Hi_b, lamb_b = np.loadtxt(
            f"{self.output_dir}/FL_switch2.dat", unpack=True, skiprows=1
        )
        w_back = integrate.simpson(Hf_b - Hi_b, x=(1 - lamb_b))

        work = (w_forw - w_back) / 2

        # array of omegas, one per component
        omega = (
            np.sqrt(self.spring_constants * EV / (self.mass * MU)) * 1.0e10
        )  # [rad/s].

        natoms = len(self.supercell)

        # array of harmonic free energies, one per component
        f_harm = (
            self.concentration
            * 3
            * KB
            * self.temperature_K
            * np.log(HBAR * omega / (KB * self.temperature_K))
        )  # [eV/atom].

        total_mass = np.sum(natoms * self.concentration * self.mass)

        # Center of mass correction
        f_cm = (
            KB
            * self.temperature_K
            * np.log(
                (natoms / self.volume)
                * (
                    (2 * np.pi * KB * self.temperature_K)
                    / (np.sum(natoms * self.concentration * total_mass**2 * self.spring_constants * EV
                              / (self.mass*MU)**2))
                )
                ** (3 / 2)
            )
            / natoms
        )

        pv_term = (
            self.pressure * self.BAR_TO_PA * self.volume * self.ANGSTROM3_TO_M3 * self.JOULE_TO_EV
        ) / natoms

        free_energy = np.sum(f_harm) - work + f_cm + pv_term

        return free_energy

    def _plot_equilibration(self) -> None:
        """Generate equilibration diagnostic plots (volume, temperature, box dimensions)."""
        data = np.loadtxt(f"{self.output_dir}/equilibration.dat", unpack=True)
        step, vol, temp, lx, ly, lz, xy, xz, yz = data
        time_ps = step / 1000

        # Volume plot
        self._save_plot(
            time_ps, vol,
            xlabel="Time [ps]", ylabel=r"Volume [$\AA^3$]",
            filename="volume.pdf"
        )

        # Temperature plot
        self._save_plot(
            time_ps, temp,
            xlabel="Time [ps]", ylabel="T [K]",
            filename="temp.pdf"
        )

        # Box lengths plot
        self._save_multiline_plot(
            time_ps, [lx, ly, lz], ["lx", "ly", "lz"],
            xlabel="Time [ps]", ylabel=r"$\AA$",
            filename="box_lengths.pdf"
        )

        # Tilt factors plot
        self._save_multiline_plot(
            time_ps, [xy, xz, yz], ["xy", "xz", "yz"],
            xlabel="Time [ps]", ylabel="Tilt Factor",
            filename="tilt_factors.pdf"
        )

    def _plot_frenkel_ladd(self) -> None:
        """Generate Frenkel-Ladd switching plots."""
        # Forward switch (potential --> springs)
        pe1, _, lamb1 = np.loadtxt(
            f"{self.output_dir}/FL_switch1.dat", unpack=True
        )
        # Reverse switch (springs --> potential)
        pe2, _, lamb2 = np.loadtxt(
            f"{self.output_dir}/FL_switch2.dat", unpack=True
        )

        fig, ax = plt.subplots()
        ax.plot(lamb1, pe1, color='blue', marker='none', ls='-', label='Forward')
        ax.plot(lamb2, pe2, color='red', marker='none', ls='-', label='Reverse')
        ax.set_xlabel(r'$\lambda$')
        ax.set_ylabel('Potential Energy [eV/atom]')
        ax.legend(loc='best')
        fig.savefig(f'{self.output_dir}/E_vs_lambda.pdf', bbox_inches='tight')
        plt.close(fig)

    def _save_plot(
        self, x: np.ndarray, y: np.ndarray,
        xlabel: str, ylabel: str, filename: str
    ) -> None:
        """Save a simple line plot to file.
        
        Args:
            x: X-axis data.
            y: Y-axis data.
            xlabel: X-axis label.
            ylabel: Y-axis label.
            filename: Output filename (saved to output_dir).
        """
        fig, ax = plt.subplots()
        ax.plot(x, y, marker='none', ls='-')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        fig.savefig(f"{self.output_dir}/{filename}", bbox_inches='tight')
        plt.close(fig)

    def _save_multiline_plot(
        self, x: np.ndarray, ys: list, labels: list,
        xlabel: str, ylabel: str, filename: str
    ) -> None:
        """Save a multi-line plot to file.
        
        Args:
            x: X-axis data.
            ys: List of Y-axis data arrays.
            labels: List of labels for each line.
            xlabel: X-axis label.
            ylabel: Y-axis label.
            filename: Output filename (saved to output_dir).
        """
        fig, ax = plt.subplots()
        for y, label in zip(ys, labels):
            ax.plot(x, y, marker='none', ls='-', label=label)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.legend(loc='best')
        fig.savefig(f"{self.output_dir}/{filename}", bbox_inches='tight')
        plt.close(fig)
    
    def _modify_run_length_control(self) -> None:
        """Generate run_length_control.py with calculated accuracy values based on cell geometry."""
        # Start accuracy lists (temperature, volume, x, y, and z are normal)
        relative_accuracy = [self.DEFAULT_RELATIVE_ACCURACY] * 5
        absolute_accuracy = [None] * 5

        # Get cell parameters and add appropriate values to accuracy lists
        # get_cell_lengths_and_angles() returns angles (in degrees) in place of tilt factors.
        # Angle = 90 --> tilt factor = 0.0.
        [_, _, _, yz_angle, xz_angle, xy_angle] = self.supercell.get_cell_lengths_and_angles()
        
        # Process each cell angle and set appropriate accuracy values
        for angle in [xy_angle, xz_angle, yz_angle]:
            is_orthogonal = abs(90 - angle) < self.ORTHOGONAL_THRESHOLD_DEGREES
            relative_accuracy.append(None if is_orthogonal else self.DEFAULT_RELATIVE_ACCURACY)
            absolute_accuracy.append(self.DEFAULT_RELATIVE_ACCURACY if is_orthogonal else None)
        
        # Read the template content from the original run_length_control.py
        original_rlc_file = Path(__file__).parent / "run_length_control.py"
        with open(original_rlc_file, 'r', encoding='utf-8') as file:
            template_content = file.read()
        
        # Set RELATIVE_ACCURACY
        new_content = re.sub(
            r"RELATIVE_ACCURACY: Sequence\[Optional\[float\]\]\s*=\s*\[.*?\]",
            f"RELATIVE_ACCURACY: Sequence[Optional[float]] = {relative_accuracy}",
            template_content
        )

        # Set ABSOLUTE_ACCURACY
        new_content = re.sub(
            r"ABSOLUTE_ACCURACY: Sequence\[Optional\[float\]\]\s*=\s*\[.*?\]",
            f"ABSOLUTE_ACCURACY: Sequence[Optional[float]] = {absolute_accuracy}",
            new_content
        )

        # Set INITAL_RUN_LENGTH
        new_content = re.sub(
            r"INITIAL_RUN_LENGTH: int\s*=\s*.*",
            f"INITIAL_RUN_LENGTH: int = {self.rlc_initial_run_length}",
            new_content
        )

        # Set MINIMUM_NUMBER_OF_INDEPENDENT_SAMPLES
        new_content = re.sub(
            r"MINIMUM_NUMBER_OF_INDEPENDENT_SAMPLES: Optional\[int\]\s*=\s*.*",
            f"MINIMUM_NUMBER_OF_INDEPENDENT_SAMPLES: Optional[int] = {self.rlc_min_samples}",
            new_content
        )
        
        # Write the modified content to the temporary file
        with open(f"{self.output_dir}/run_length_control.py", 'w', encoding='utf-8') as file:
            file.write(new_content)

    def _reduce_and_average_supercell(self, atoms_npt_path: str, reduced_atoms_save_path: str) -> Atoms:
        """Reduce supercell to unit cell by averaging equivalent atom positions.
        
        Args:
            atoms_npt_path: Path to LAMMPS data file with average positions.
            reduced_atoms_save_path: Path to save the reduced atoms structure.
            
        Returns:
            Reduced ASE Atoms object.
        """
        # Read lammps data file of average positions
        atoms_npt = read(atoms_npt_path, format='lammps-data')

        # Reduce to unit cell
        reduced_atoms = reduce_and_avg(atoms_npt, self.repeat)

        # Print reduced_atoms for verification
        reduced_atoms.write(reduced_atoms_save_path, format='lammps-data', masses=True, atom_style=self.atom_style)

        return reduced_atoms
    
    def _collect_isolated_atom_energies(self, reduced_atoms: Atoms) -> float:
        """Compute the isolated atom energy per atom for the given structure.
        
        Args:
            reduced_atoms: ASE Atoms object with the reduced unit cell.
            
        Returns:
            Isolated atom energy per atom in eV.
        """
        # List of unique elements (strings) in the "reduced_atoms" ASE atoms object
        element_list = list(dict.fromkeys(reduced_atoms.get_chemical_symbols()))

        # stoichiometry is already set in _calculate, use it directly
        # (prototype_label was set at the start of _calculate)
        self.stoichiometry = get_stoich_reduced_list_from_prototype(self.prototype_label)
        
        # Compute total isolated energy for one formula unit
        total_isolated_energy = sum(
            stoich * get_isolated_energy_per_atom(self.kim_model_name, element)
            for element, stoich in zip(element_list, self.stoichiometry)
        )
        
        # Divide by total number of atoms in formula to get per-atom energy
        num_atoms_in_formula = sum(self.stoichiometry)
        isolated_atom_energy = total_isolated_energy / num_atoms_in_formula
        
        return isolated_atom_energy

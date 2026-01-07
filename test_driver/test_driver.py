import multiprocessing
import os
import re
from pathlib import Path
from typing import List, Tuple, Sequence

import numpy as np
import scipy.constants as sc
from ase import Atoms
from ase.data import atomic_masses, atomic_numbers
from ase.io import read
from kim_tools import KIMTestDriverError, get_isolated_energy_per_atom, get_stoich_reduced_list_from_prototype
from kim_tools.test_driver import SingleCrystalTestDriver
from kim_tools.symmetry_util.core import reduce_and_avg

import matplotlib.pyplot as plt

from .lammps_template import LammpsTemplate
from .helper_functions import run_lammps

EV = sc.value("electron volt")
MU = sc.value("atomic mass constant")
HBAR = sc.value("Planck constant in eV/Hz") / (2 * np.pi)
KB = sc.value("Boltzmann constant in eV/K")


class TestDriver(SingleCrystalTestDriver):
    def _calculate(
        self,
        timestep_ps: float = 0.001,
        FL_switch_timesteps: int = 50000,
        FL_equil_timesteps: int = 10000,
        number_sampling_timesteps: int = 100,
        target_size: int = 1000,
        repeat: Sequence[int] = (0, 0, 0),
        lammps_command = "lmp",
        msd_threshold_angstrom_squared_per_sampling_timesteps: float = 0.1,
        number_msd_timesteps: int = 10000,
        number_avePOS_timesteps: int = 30000,
        random_seed: int = 101010,
        rlc_N_every: int = 10,
        rlc_inital_run_length: int = 500,
        rlc_min_samples: int = 100,
        output_dir: str = "lammps_output",
        equilibration_plots: bool = True,
        FL_plots: bool = True,
        **kwargs) -> None:

        """
        Gibbs free energy of a crystal at constant temperature and stress using a nonequilibrium thermodynamic integration method implemented in LAMMPS (doi = 10.1016/j.commatsci.2015.10.050).

        Args:
        """

        # Set prototype label
        self.prototype_label = self._get_nominal_crystal_structure_npt()["prototype-label"]["source-value"]

        self.test_driver_directory = os.path.dirname(os.path.realpath(__file__))
        print(self.test_driver_directory)

        self.temperature_K = self._get_temperature(unit="K")
        self.cauchy_stress = self._get_cell_cauchy_stress(unit='bars')
        self.pressure = -self.cauchy_stress[0] 
        self.atoms = self._get_atoms()
        self.timestep_ps = timestep_ps
        self.FL_switch_timesteps = FL_switch_timesteps
        self.FL_equil_timesteps = FL_equil_timesteps
        self.number_sampling_timesteps = number_sampling_timesteps
        self.target_size = target_size
        self.repeat = repeat
        self.lammps_command = lammps_command
        self.rlc_N_every = rlc_N_every
        self.rlc_inital_run_length = rlc_inital_run_length
        self.rlc_min_samples = rlc_min_samples
        self.msd_threshold_angstrom_squared_per_sampling_timesteps = msd_threshold_angstrom_squared_per_sampling_timesteps
        self.number_msd_timesteps = number_msd_timesteps
        self.number_avePOS_timesteps = number_avePOS_timesteps
        self.random_seed = random_seed
        self.output_dir = output_dir

        self.atom_style = self._get_supported_lammps_atom_style()
        symbols = self.atoms.get_chemical_symbols()
        species = sorted(set(symbols))

        self._validate_inputs()

        self.supercell = self._setup_initial_structure(filename=f"{self.output_dir}/zero_temperature_crystal.data")

        # Modify the variables in run_length_control.py
        self._modify_run_length_control()

        # Write initial template file "in.lammps"
        self.templates = LammpsTemplate(root=f"{self.output_dir}/")
        self.templates._write_lammps_file(
            nspecies=len(self.species), is_triclinic=self.is_triclinic
        )

        # Run LAMMPS
        log_filename, restart_filename, self.spring_constants, self.volume = run_lammps(
            self.kim_model_name, self.temperature_K, self.pressure, self.timestep_ps, self.FL_switch_timesteps, self.FL_equil_timesteps, self.number_sampling_timesteps, species,
            self.msd_threshold_angstrom_squared_per_sampling_timesteps, self.number_msd_timesteps, self.number_avePOS_timesteps,
            self.rlc_N_every, self.lammps_command, self.random_seed, self.output_dir, self.test_driver_directory)

        # Check that LAMMPS ran to completion
        if self._check_if_lammps_ran_to_completion(lammps_log=f"{self.output_dir}/free_energy.log") == str("Crystal melted or vaporized"):
            raise Exception("Crystal melted or vaporized")

        # Compute free-energy
        free_energy_per_atom = self._free_energy()

        # Make equilibration plots
        if equilibration_plots:

            step, vol, temp, lx, ly, lz, xy, xz, yz = np.loadtxt(f"{self.output_dir}/equilibration.dat", unpack=True)

            # ============================== Volume ============================== #

            fig = plt.figure()

            plt.plot(step/1000, vol, marker='none', ls='-')

            plt.xlabel("Time [ps]")
            plt.ylabel(r"Volume [$\AA^3$]")

            plt.savefig(f"{self.output_dir}/volume.pdf", bbox_inches='tight')
            plt.close(fig)

            # ============================== Temp ============================== #

            fig = plt.figure()

            plt.plot(step/1000, temp, marker='none', ls='-')

            plt.xlabel("Time [ps]")
            plt.ylabel("T [K]")

            plt.savefig(f"{self.output_dir}/temp.pdf", bbox_inches='tight')
            plt.close(fig)

            # ============================== lx, ly, lz ============================== #

            fig = plt.figure()

            plt.plot(step/1000, lx, marker='none', ls='-', label='lx')
            plt.plot(step/1000, ly, marker='none', ls='-', label='ly')
            plt.plot(step/1000, lz, marker='none', ls='-', label='lz')

            plt.xlabel("Time [ps]")
            plt.ylabel(r"$\AA$")

            plt.legend(loc='best')

            plt.savefig(f"{self.output_dir}/box_lengths.pdf", bbox_inches='tight')
            plt.close(fig)

            # ============================== lx, ly, lz ============================== #

            fig = plt.figure()

            plt.plot(step/1000, xy, marker='none', ls='-', label='xy')
            plt.plot(step/1000, xz, marker='none', ls='-', label='xz')
            plt.plot(step/1000, yz, marker='none', ls='-', label='yz')

            plt.xlabel("Time [ps]")
            plt.ylabel("Tilt Factor")

            plt.legend(loc='best')

            plt.savefig(f"{self.output_dir}/tilt_factors.pdf", bbox_inches='tight')
            plt.close(fig)

        # Make Frenkel-Ladd plots
        if FL_plots:

            PE1, PE1_springs, lamb1 = np.loadtxt(f"{self.output_dir}/FL_switch1.dat", unpack=True) # Forward (potential --> springs)
            PE2, PE2_springs, lamb2 = np.loadtxt(f"{self.output_dir}/FL_switch2.dat", unpack=True) # Reverse (springs --> potential)

            fig = plt.figure()

            plt.plot(lamb1, PE1, color='blue', marker='none', ls='-', label='Forward')
            plt.plot(lamb2, PE2, color='red', marker='none', ls='-', label='Reverse')

            plt.xlabel(r'$\lambda$')
            plt.ylabel('Potential Energy [eV/atom]')

            plt.legend(loc='best')

            plt.savefig(f'{self.output_dir}/E_vs_lambda.pdf', bbox_inches='tight')
            plt.close(fig)

        self._add_property_instance_and_common_crystal_genome_keys("crystal-structure-npt", write_temp=True, write_stress=True)#, stress_unit="bars")
        self._add_file_to_current_property_instance("restart-file", f"{self.output_dir}/free_energy.restart")

        reduced_atoms = self._reduce_average_and_verify_symmetry(atoms_npt=f"{self.output_dir}/free_energy.data", reduced_atoms_save_path=f"{self.output_dir}/reduced_atoms.data")
        self._update_nominal_parameter_values(reduced_atoms)

        # Collect the energies of isolated atoms to subtract from final values
        isolated_atom_energy = self._collect_isolated_atom_energies(reduced_atoms)
        free_energy_per_atom = free_energy_per_atom - isolated_atom_energy

        # Convert to eV/formula (originally in eV/atom)
        num_atoms_in_formula = sum(self.stoichiometry)
        free_energy_per_formula = free_energy_per_atom * num_atoms_in_formula

        # Convert to eV/amu (originally in eV/atom)
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

        # Write keys to property
        self._add_property_instance_and_common_crystal_genome_keys(
            "free-energy", write_stress=True, write_temp=self.temperature_K
        )
        self._add_key_to_current_property_instance(
            "gibbs-free-energy-per-atom", free_energy_per_atom, "eV/atom"
        )
        self._add_key_to_current_property_instance(
            "gibbs-free-energy-per-formula", free_energy_per_formula, "eV/formula"
        )
        self._add_key_to_current_property_instance(
            "specific-gibbs-free-energy", specific_free_energy, "eV/amu"
        )

    def _validate_inputs(self):

        if not self.temperature_K > 0.0:
            raise ValueError("Temperature has to be larger than zero.")

        if not len(self.cauchy_stress) == 6:
            raise ValueError("Specify all six (x, y, z, xy, xz, yz) entries of the cauchy stress tensor.")

        if not (self.cauchy_stress[0] == self.cauchy_stress[1] == self.cauchy_stress[2]):
            raise ValueError("The diagonal entries of the stress tensor have to be equal so that a hydrostatic "
                             "pressure is used.")

        if not (self.cauchy_stress[3] == self.cauchy_stress[4] == self.cauchy_stress[5] == 0.0):
            raise ValueError("The off-diagonal entries of the stress tensor have to be zero so that a hydrostatic "
                             "pressure is used.")

        if not self.timestep_ps > 0.0:
            raise ValueError("Timestep has to be larger than zero.")

        if not self.number_sampling_timesteps > 0:
            raise ValueError("Number of timesteps between sampling in Lammps has to be bigger than zero.")

        if not len(self.repeat) == 3:
            raise ValueError("The repeat argument has to be a tuple of three integers.")

        if not all(r >= 0 for r in self.repeat):
            raise ValueError("All number of repeats must be bigger than zero.")

        if not self.msd_threshold_angstrom_squared_per_sampling_timesteps > 0.0:
            raise ValueError("The mean-squared displacement threshold has to be bigger than zero.")

        if not self.number_msd_timesteps > 0:
            raise ValueError("The number of timesteps to monitor the mean-squared displacement has to be bigger than "
                             "zero.")

        if not self.number_msd_timesteps % self.number_sampling_timesteps == 0:
            raise ValueError("The number of timesteps to monitor the mean-squared displacement has to be a multiple of "
                             "the number of sampling timesteps.")

        if not self.random_seed > 0:
            raise ValueError("The random seed has to be bigger than zero.")

    def _setup_initial_structure(self, filename: Path) -> Atoms:
        # Copy original atoms so that their information does not get lost when the new atoms are modified.

        atoms_new = self.atoms.copy()

        # Build supercell
        if self.repeat == (0,0,0):
            # Get a size close to 10K atoms (shown to give good convergence)
            x = int(np.ceil(np.cbrt(self.target_size / len(self.atoms))))
            self.repeat = (x,x,x)

        atoms_new = atoms_new.repeat(self.repeat)

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

        self.concentration = np.zeros(len(self.species))
        symbols = np.array(self.atoms.get_chemical_symbols())
        for i in range(len(self.species)):

            self.concentration[i] = (symbols == self.species[i]).sum()
        self.concentration *= 1 / len(symbols)

        atoms_new.write(filename, format="lammps-data", masses=True)
        self.zero_k_structure_path = filename

        # Check if the cell is triclinic by checking if any of the cell lengths or angles are different from the others.
        # tol = 1e-2
        # a, b, c, alpha, beta, gamma = cell_to_cellpar(self.atoms.cell)
        # self.is_triclinic = all(abs(x - y) > 1e-2 for x, y in [(a, b), (b, c), (a, c), (alpha, beta), (beta, gamma), (alpha, gamma)]) or all(abs(x - 90) > tol for x in [alpha, beta, gamma])

        # Treat all cells as triclinic for simplicity.
        self.is_triclinic = True
       
        return atoms_new
    
    def _free_energy(self) -> float:

        Hi_f, Hf_f, lamb_f = np.loadtxt(
            f"{self.output_dir}/FL_switch1.dat", unpack=True, skiprows=1
        )
        W_forw = np.trapz(Hf_f - Hi_f, lamb_f)

        Hf_b, Hi_b, lamb_b = np.loadtxt(
            f"{self.output_dir}/FL_switch2.dat", unpack=True, skiprows=1
        )
        W_back = np.trapz(Hf_b - Hi_b, 1 - lamb_b)

        Work = (W_forw - W_back) / 2
        Dissipation = (W_forw + W_back) / 2

        # array of omegas, one per component
        omega = (
            np.sqrt(self.spring_constants * EV / (self.mass * MU)) * 1.0e10
        )  # [rad/s].

        natoms = len(self.supercell)

        # array of harmoinc free energies, one per component
        F_harm = (
            self.concentration
            * 3
            * KB
            * self.temperature_K
            * np.log(HBAR * omega / (KB * self.temperature_K))
        )  # [eV/atom].

        total_mass = np.sum(natoms * self.concentration * self.mass)

        F_CM = (
            KB
            * self.temperature_K
            * np.log(
                (natoms / self.volume)
                * (
                    (2 * np.pi * KB * self.temperature_K)
                    #/ (natoms * self.concentration * self.mass * omega**2)
                    / (np.sum(natoms * self.concentration * total_mass**2 * self.spring_constants * EV
                              / (self.mass*MU)**2)) # Khanna 2021, J. Chem. Phys., eq. 10
                )
                ** (3 / 2)
            )
            / natoms
        )  # correction for fixed center of mass

        bar_to_Pa = 1e5
        A3_to_m3 = 1e-30
        J_to_eV = 6.2415e18

        PV_term = (
            self.pressure * bar_to_Pa * self.volume * A3_to_m3 * J_to_eV
        ) / natoms

        free_energy = np.sum(F_harm) - Work + F_CM + PV_term

        return free_energy

    def _check_if_lammps_ran_to_completion(self, lammps_log: Path):
        try:
            with open(lammps_log, "r") as file:
                lines = file.readlines()

                if not lines:
                    return False
                
                last_line = lines[-2].strip()

                if last_line.startswith("Crystal melted or vaporized"):
                    return str("Crystal melted or vaporized")

                last_line = lines[-1].strip()

                return last_line.startswith("Total wall time:")

        except FileNotFoundError:
            return str(f"LAMMPS log does not exist.")
    
    def _modify_run_length_control(self):
        # Start accuracy lists (temperature, volume, x, y, and z are normal)
        relative_accuracy = [0.01, 0.01, 0.01, 0.01, 0.01]
        absolute_accuracy = [None, None, None, None, None]

        # Get cell parameters and add appropriate values to accuracy lists (0.1 and None for non-zero tilt factors, vice-versa for zero)
        # get_cell_lengths_and_angles() returns angles (in degrees) in place of tilt factors. Angle = 90 --> tilt factor = 0.0.
        # The criterion for an orthogonal tilt factor ("abs(90-angle) < 0.1") can be modified, depending on how small of a tilt factor is too small for kim-convergence
        [X_cell, Y_cell, Z_cell, YZ_cell, XZ_cell, XY_cell] = self.supercell.get_cell_lengths_and_angles()
        
        # Define threshold for considering angles as orthogonal (in degrees).
        # This should probably be replaced with something like self._get_supported_lammps_atom_style(), but to check orthogonality
        ORTHOGONAL_THRESHOLD = 0.1
        
        # Process each cell angle and set appropriate accuracy values
        for angle in [XY_cell, XZ_cell, YZ_cell]:
            is_orthogonal = abs(90 - angle) < ORTHOGONAL_THRESHOLD
            relative_accuracy.append(None if is_orthogonal else 0.01)
            absolute_accuracy.append(0.01 if is_orthogonal else None)
        
        # Read the template content from the original run_length_control.py
        original_rlc_file = Path(__file__).parent / "run_length_control.py"
        with open(original_rlc_file, 'r') as file:
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
            f"INITIAL_RUN_LENGTH: int = {self.rlc_inital_run_length}",
            new_content
        )

        # Set MINIMUM_NUMBER_OF_INDEPENDENT_SAMPLES
        new_content = re.sub(
            r"MINIMUM_NUMBER_OF_INDEPENDENT_SAMPLES: Optional\[int\]\s*=\s*.*",
            f"MINIMUM_NUMBER_OF_INDEPENDENT_SAMPLES: Optional[int] = {self.rlc_min_samples}",
            new_content
        )
        
        # Write the modified content to the temporary file
        with open(original_rlc_file, 'w') as file:
            file.write(new_content)

    def _reduce_average_and_verify_symmetry(self, atoms_npt: Path, reduced_atoms_save_path: Path):

        # Read lammps data file of average positions
        atoms_npt = read(atoms_npt, format='lammps-data')

        # Reduce to unit cell
        reduced_atoms = reduce_and_avg(atoms_npt, self.repeat)

        # Print reduced_atoms for verification
        reduced_atoms.write(reduced_atoms_save_path, format='lammps-data', masses=True)

        # Check symmetry
        if not self._verify_unchanged_symmetry(reduced_atoms):
            raise ValueError("The symmetry of the atoms have changed.")

        return reduced_atoms
    
    def _collect_isolated_atom_energies(self,reduced_atoms):
        isolated_atom_energy_list = []

        # List of unique elements (strings) in the "reduced_atoms" ASE atoms object
        element_list = list(dict.fromkeys(reduced_atoms.get_chemical_symbols()))

        # List of stoichiometric coefficients (should line up with "element_list" if everything is in alphabetical order)
        self.prototype_label = self._get_nominal_crystal_structure_npt()["prototype-label"][
            "source-value"
        ]
        self.stoichiometry = get_stoich_reduced_list_from_prototype(self.prototype_label)
        
        # Append correct number of each isolated atom energy according to stoichiometric coefficients
        for element,stoich in zip(element_list,self.stoichiometry):
            isolated_atom_energy_list.append(stoich * get_isolated_energy_per_atom(self.kim_model_name, element)) # appends as many as exist in one formula unit

        # Take the mean of the energies to subtract from the per-atom free energy
        isolated_atom_energy = np.mean(isolated_atom_energy_list)
        return isolated_atom_energy

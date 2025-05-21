import re
import subprocess
from typing import List, Tuple

import numpy as np
import scipy.constants as sc
from ase import Atoms
from ase.data import atomic_masses, atomic_numbers
from ase.io import read
from kim_tools import (
    SingleCrystalTestDriver,
    get_isolated_energy_per_atom,
    get_stoich_reduced_list_from_prototype,
)

from .helper_functions import reduce_and_avg, test_reduced_distances
from .lammps_templates import LammpsTemplates

EV = sc.value("electron volt")
MU = sc.value("atomic mass constant")
HBAR = sc.value("Planck constant in eV/Hz") / (2 * np.pi)
KB = sc.value("Boltzmann constant in eV/K")


class TestDriver(SingleCrystalTestDriver):
    def _calculate(
        self,
        size: Tuple[int, int, int] = (0,0,0),
        **kwargs,
    ) -> None:
        """Gibbs free energy of a crystal at constant temperature and pressure using Frenkel-Ladd Hamiltonian integration algorithm. Computed through one equilibrium NPT simulation ('preFL') and one NONequilibrium NVT simulation ('FL').

        Args:
            size (Tuple[int, int, int]): system size. By default, the system size is computed to have ~10,000 atoms.
        """

        self.temperature_K = self._get_temperature(unit="K")
        self.cauchy_stress = self._get_cell_cauchy_stress(unit='atm')
        self.pressure = -self.cauchy_stress[0] 
        self.atoms = self._get_atoms()
        self.size = size
        
        self._validate_inputs()

        
        self.supercell = self._setup_initial_structure()

        self._modify_accuracies()

        # Write initial template file
        self.templates = LammpsTemplates(root="output/lammps_templates/")
        self.templates._write_pre_fl_lammps_templates(
            nspecies=len(self.species), is_triclinic=self.is_triclinic
        )

        # preFL computes the equilibrium lattice parameter and spring constants for a given temperature and pressure.
        equilibrium_cell, self.spring_constants, self.volume, self.atom_style =self._run_preFL()
        

        # If the crystal melts or vaporizes, kim-convergence may run indefinately.
        # If lammps detects diffusion, it prints a specific string, then quits.
        # The 'if' statement below handles that case.
        # This is a first implementation. It's probably cleaner to leave the responsibility of quitting to some standardized function that is common across finite-temperature test-drivers
        self._check_diffusion(lammps_log="output/lammps_preFL.log")

        reduced_atoms_preFL = self._reduce_average_and_verify_symmetry(atoms_npt="output/lammps_preFL.data",  reduced_atoms_save_path="output/reduced_atoms_preFL.data")
        
        self._update_nominal_parameter_values(reduced_atoms_preFL,max_resid=1e-5)

        # crystal-structure-npt
        self._add_property_instance_and_common_crystal_genome_keys("crystal-structure-npt", write_temp=True, write_stress=True)
        self._add_file_to_current_property_instance("restart-file","output/lammps_preFL.restart")
    
        
        # Rescaling 0K supercell to have equilibrium lattice constant.
        # equilibrium_cell is 3x3 matrix or can also have [len(a), len(b), len(c), angle(b,c), angle(a,c), angle(a,b)]
        self.supercell.set_cell(equilibrium_cell, scale_atoms=True)
        self.supercell.write(
            "output/equilibrium_crystal.data",
            format="lammps-data",
            masses=True,
            atom_style=self.atom_style,
        )
        # self._add_file_to_current_property_instance("data-file","output/equilibrium_crystal.data")

        # Collect the energies of isolated atoms to subtract from final values
        isolated_atom_energy = self._collect_isolated_atom_energies(reduced_atoms_preFL)

        # FL computes the free energy at a given pressure and temperature.
        self.templates._write_fl_lammps_templates(
            spring_constants=self.spring_constants
        )
        
        free_energy_per_atom = self._FL() - isolated_atom_energy

        
        self._check_diffusion(lammps_log="output/lammps_FL.log")

        reduced_atoms_FL = self._reduce_average_and_verify_symmetry(atoms_npt="output/lammps_FL.data",  reduced_atoms_save_path="output/reduced_atoms_FL.data")

        self._update_nominal_parameter_values(reduced_atoms_FL,max_resid=1e-5)

        self._add_property_instance_and_common_crystal_genome_keys("crystal-structure-npt", write_temp=True, write_stress=True)
        self._add_file_to_current_property_instance("restart-file","output/lammps_FL.restart")

        # Convert to eV/formula (originally in eV/atom)
        # get_stoich_reduced_list_from_prototype returns a list corresponding to the stoichiometry of the prototype label, e.g. "A2B_hP9_152_c_a" -> [2,1]
        num_atoms_in_formula = sum(self.stoichiometry)
        free_energy_per_formula = free_energy_per_atom * num_atoms_in_formula

        # Convert to eV/amu (originally in eV/atom)
        atoms_per_cell = len(reduced_atoms_preFL.get_masses())
        mass_per_cell = sum(reduced_atoms_preFL.get_masses())
        free_energy_per_cell = free_energy_per_atom * atoms_per_cell
        specific_free_energy = free_energy_per_cell / mass_per_cell

        # Print results
        print("####################################")
        print("# Frenkel Ladd Free Energy Results #")
        print("####################################")

        print(f"G_FL = {free_energy_per_atom:.5f} (eV/atom)")

        # Write keys to property
        self._add_property_instance_and_common_crystal_genome_keys(
            "free-energy", write_stress=True, write_temp=True
        )
        self._add_key_to_current_property_instance(
            "gibbs_free_energy_per_atom", free_energy_per_atom, "eV/atom"
        )
        self._add_key_to_current_property_instance(
            "gibbs_free_energy_per_formula", free_energy_per_formula, "eV/formula"
        )
        self._add_key_to_current_property_instance(
            "specific_gibbs_free_energy", specific_free_energy, "eV/amu"
        )

    def _validate_inputs(self):
        if not self.temperature_K > 0.0:
            raise ValueError("Temperature has to be larger than zero.")

        if not self.pressure >= 0.0:
            raise ValueError("Pressure has to be greater than or equal to zero.")

        if not np.allclose(self.cauchy_stress[:3], self.cauchy_stress[0], rtol=1e-5, atol=1e-8):
            raise ValueError("The first three components of the Cauchy stress tensor must be equal (calculation are run at constant isotropic stress).")

    def _setup_initial_structure(
        self,
        filename: str = "output/zero_temperature_crystal.data",
    ) -> Atoms:
        # Copy original atoms so that their information does not get lost when the new atoms are modified.

        atoms_new = self.atoms.copy()

        # Build supercell
        if self.size == (0,0,0):
            # Get a size close to 10K atoms (shown to give good convergence)
            x = int(np.ceil(np.cbrt(10_000 / len(self.atoms))))
            self.size = (x,x,x)

        atoms_new = atoms_new.repeat(self.size)

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

        angles = self.atoms.get_cell().angles()
        self.is_triclinic = (
            all(angle != 90 for angle in angles) and len(set(angles)) == 3
        )

        return atoms_new

    def _run_preFL(self):
        equilibrium_cell, spring_constants, volume = self._preFL()

        # Some models want atom_style="charge", others want "atomic"
        # We tried with 'atomic', if it fails, try 'charge'
        atom_style = "atomic"
        if not self._check_if_lammps_run_to_completiton(
            lammps_log="output/lammps_preFL.log"
        ):
            atom_style = "charge"
            self.supercell.write(
                self.zero_k_structure_path,
                format="lammps-data",
                masses=True,
                atom_style=atom_style,
            )
            equilibrium_cell, spring_constants, volume = self._preFL()

        assert len(self.species) == len(spring_constants)
        return equilibrium_cell, spring_constants, volume, atom_style
    
    def _preFL(self) -> Tuple[List[float], List[float]]:
        variables = {
            "modelname": self.kim_model_name,
            "temperature": self.temperature_K,
            "temperature_damping": 0.1, # ps
            "temperature_seed": np.random.randint(low=100000, high=999999, dtype=int),
            "pressure": self.pressure,
            "pressure_damping": 1.0, # ps
            "timestep": 0.001,  # ps
            "species": " ".join(self.species),
            "output_filename": "output/lammps_preFL.dat",
            "write_restart_filename": "output/lammps_preFL.restart",
            "write_data_filename": "output/lammps_preFL.data",
        }
        # TODO: Possibly run MPI version of Lammps if available.
        command = (
            "lammps "
            + " ".join(f"-var {key} '{item}'" for key, item in variables.items())
            + " -log output/lammps_preFL.log"
            + f" -in {self.templates.root}preFL_template.lmp"
        )
        try:
            subprocess.run(command, check=True, shell=True)
        except:
            return (None, None, None)

        # handling the case where it did not create lammps_preFL (e.g, when atomic style needs to be charge)

        # Analyse lammps outputs
        data = np.loadtxt("output/lammps_preFL.dat", unpack=True)
        # xx, xy, xz, yx, yy, yz, zx, zy, zz, spring_constants = data

        (lx, ly, lz, xy, yz, xz, volume), spring_constants = (
            data[: -len(self.species)],
            data[-len(self.species) :],
        )

        if isinstance(spring_constants, float):
            spring_constants = [spring_constants]

        equilibrium_cell = np.array([[lx, 0, 0], [xy, ly, 0], [xz, yz, lz]])

        return equilibrium_cell, np.array(spring_constants), volume

    def _FL(
        self,
    ) -> float:

        variables = {
            "modelname": self.kim_model_name,
            "temperature": self.temperature_K,
            "pressure": self.pressure,
            "species": " ".join(self.species),
            "t_switch": 10000,
            "temperature_damping": 0.1,
            "temperature_seed": np.random.randint(low=100000, high=999999, dtype=int),
            "t_equil": 10000,
            "timestep": 0.001,  # ps
            "spring_constant": self.spring_constants,
            "output_filename": "output/lammps_FL.dat",
            "write_restart_filename": "output/lammps_FL.restart",
            "switch1_output_file": "output/FL_switch1.dat",
            "switch2_output_file": "output/FL_switch2.dat",
        }
        # TODO: Possibly run MPI version of Lammps if available.
        command = (
            "lammps "
            + " ".join(f"-var {key} '{item}'" for key, item in variables.items())
            + " -log output/lammps_FL.log"
            + f" -in {self.templates.root}FL_template.lmp"
        )
        subprocess.run(command, check=True, shell=True)

        return self._compute_free_energy()

    def _compute_free_energy(self) -> float:
        """Compute free energy via integration of FL path"""

        Hi_f, Hf_f, lamb_f = np.loadtxt(
            "output/FL_switch1.dat", unpack=True, skiprows=1
        )
        W_forw = np.trapz(Hf_f - Hi_f, lamb_f)

        Hf_b, Hi_b, lamb_b = np.loadtxt(
            "output/FL_switch2.dat", unpack=True, skiprows=1
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

 
    def _check_if_lammps_run_to_completiton(self, lammps_log: str):
        try:
            with open(lammps_log, "r") as file:
                lines = file.readlines()

                if not lines:
                    return False

                last_line = lines[-1].strip()

                return last_line.startswith("Total wall time:")

        except FileNotFoundError:
            print(f"The file {lammps_log} does not exist.")
            return False
    
    def _check_diffusion(self, lammps_log: str):
     
        try:
            with open(lammps_log, "r") as file:
                lines = file.readlines()

                if not lines:
                    return

                last_line = lines[-2].strip()

                if last_line.startswith("Crystal melted or vaporized"):
                    raise ValueError("Crystal melted or vaporized")

        except FileNotFoundError:
            print(f"The file {lammps_log} does not exist.")
            raise

    
    def _modify_accuracies(self):
         # Start accuracy lists (volume, x, y, and z are normal)
        relative_accuracy = [0.01, 0.01, 0.01, 0.01]
        absolute_accuracy = [None, None, None, None]

        # Get cell parameters and add appropriate values to accuracy lists (0.01 and None for non-zero tilt factors, vice-versa for zero)
        # get_cell_lengths_and_angles() returns angles in place of tilt factors. Angle = 90 --> tilt factor = 0.0.
        # The criterion for an orthogonal tilt factor ("abs(90-angle) < 0.1") can be modified, depending on how small of a tilt factor is too small for kim-convergence
        [X_cell, Y_cell, Z_cell, YZ_cell, XZ_cell, XY_cell] = self.supercell.get_cell_lengths_and_angles()
        
        # Define threshold for considering angles as orthogonal (in degrees)
        ORTHOGONAL_THRESHOLD = 0.1
        
        # Process each cell angle and set appropriate accuracy values
        for angle in [XY_cell, XZ_cell, YZ_cell]:
            is_orthogonal = abs(90 - angle) < ORTHOGONAL_THRESHOLD
            relative_accuracy.append(None if is_orthogonal else 0.01)
            absolute_accuracy.append(0.01 if is_orthogonal else None)

        # Replace lists in "accuracies_general.py"
        new_accuracies = {
                        "RELATIVE_ACCURACY: Sequence[Optional[float]]": relative_accuracy,
                        "ABSOLUTE_ACCURACY: Sequence[Optional[float]]": absolute_accuracy
                        }
        
        with open("test_driver/accuracies.py", 'r') as file:
            content = file.read()
        
        pattern = r"RELATIVE_ACCURACY: Sequence\[Optional\[float\]\].s*=.s*\[.*?\]"
        replacement = f"RELATIVE_ACCURACY: Sequence[Optional[float]] = {new_accuracies['RELATIVE_ACCURACY: Sequence[Optional[float]]']}"
        content = re.sub(pattern, replacement, content)

        pattern = r"ABSOLUTE_ACCURACY: Sequence\[Optional\[float\]\].s*=.s*\[.*?\]"
        replacement = f"ABSOLUTE_ACCURACY: Sequence[Optional[float]] = {new_accuracies['ABSOLUTE_ACCURACY: Sequence[Optional[float]]']}"
        content = re.sub(pattern, replacement, content)

        with open("test_driver/accuracies.py", 'w') as file:
            file.write(content)

    def _reduce_average_and_verify_symmetry(self, atoms_npt: str,reduced_atoms_save_path: str):
        # Read lammps dump file of average positions
        atoms_npt = read(atoms_npt, format='lammps-data')
        # Reduce to unit cell
        reduced_atoms, reduced_distances = reduce_and_avg(atoms_npt, self.size)
        test_reduced_distances(reduced_distances)
        # Print reduced_atoms for verification
        reduced_atoms.write(reduced_atoms_save_path, format='lammps-data',masses=True)

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
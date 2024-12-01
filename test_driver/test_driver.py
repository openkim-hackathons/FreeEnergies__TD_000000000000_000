import os
import sys
import re
import subprocess
import shutil
from typing import List, Tuple
from ase.build import bulk
from ase import Atoms
from ase import io
from ase.io import read, write
import numpy as np
from ase.data import atomic_masses, atomic_numbers
from kim_tools import CrystalGenomeTestDriver, get_stoich_reduced_list_from_prototype
import scipy.constants as sc
from .lammps_templates import LammpsTemplates
from .helper_functions import reduce_and_avg
from kim_python_utils.ase import get_isolated_energy_per_atom

EV = sc.value("electron volt")
MU = sc.value("atomic mass constant")
HBAR = sc.value("Planck constant in eV/Hz") / (2 * np.pi)
KB = sc.value("Boltzmann constant in eV/K")


class TestDriver(CrystalGenomeTestDriver):
    def _calculate(
        self,
        temperature: float = 20.0,
        pressure: float = 0.0,
        size: Tuple[int, int, int] = (0,0,0),
        **kwargs,
    ) -> None:
        """Gibbs free energy of a crystal at constant temperature and pressure using Frenkel-Ladd Hamiltonian integration algorithm. Computed through one equilibrium NPT simulation ('preFL') and one NONequilibrium NVT simulation ('FL').

        Args:
            temperature (float): system temperature in K.
            pressure (float): system pressure in atm.
            size (Tuple[int, int, int]): system size.
        """
        # Check arguments
        self.temperature = temperature
        self.pressure = pressure
        self._validate_inputs()

        # Write initial atomic structure to lammps dump file
        if size == (0,0,0):
            # Get a size close to 10K atoms (shown to give good convergence)
            x = int(np.ceil(np.cbrt(10_000 / len(self.atoms))))
            size = (x,x,x)
        self.supercell = self._setup_initial_structure(size)

        # Start accuracy lists (volume, x, y, and z are normal)
        relative_accuracy = [0.01, 0.01, 0.01, 0.01]
        absolute_accuracy = [None, None, None, None]

        # Get cell parameters and add appropriate values to accuracy lists (0.01 and None for non-zero tilt factors, vice-versa for zero)
        # get_cell_lengths_and_angles() returns angles in place of tilt factors. Angle = 90 --> tilt factor = 0.0.
        # The criterion for an orthogonal tilt factor ("abs(90-angle) < 0.1") can be modified, depending on how small of a tilt factor is too small for kim-convergence
        [X_cell, Y_cell, Z_cell, YZ_cell, XZ_cell, XY_cell] = self.supercell.get_cell_lengths_and_angles()
        
        for angle in [XY_cell, XZ_cell, YZ_cell]:
            if abs(90-angle) < 0.1:
                relative_accuracy.append(None)
                absolute_accuracy.append(0.01)
            else:
                relative_accuracy.append(0.01)
                absolute_accuracy.append(None)

        # Replace lists in "accuracies_general.py"
        new_accuracies = {
                        "RELATIVE_ACCURACY: Sequence[Optional[float]]": relative_accuracy,
                        "ABSOLUTE_ACCURACY: Sequence[Optional[float]]": absolute_accuracy
                        }
        self._modify_accuracies(new_accuracies)

        # Write initial template file
        self.templates = LammpsTemplates(root="lammps_templates/")
        self.templates._write_pre_fl_lammps_templates(
            nspecies=len(self.species), is_triclinic=self.is_triclinic
        )

        # preFL computes the equilibrium lattice parameter and spring constants for a given temperature and pressure.
        # TODO: This should probably be replaced with its own test driver, which compute equilibrium lattice constants, and which can handles arbitrary crystal structures. Then we can get spring constants.
        equilibrium_cell, self.spring_constants, self.volume = self._preFL()

        # If the crystal melts or vaporizes, kim-convergence may run indefinately.
        # If lammps detects diffusion, it prints a specific string, then quits.
        # The 'if' statement below handles that case.
        # This is a first implementation. It's probably cleaner to leave the responsibility of quitting to some standardized function that is common across finite-temperature test-drivers
        if self._check_diffusion(
            lammps_log="output/lammps_preFL.log"
        ):
            sys.exit("Crystal melted or vaporized")

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
            equilibrium_cell, self.spring_constants, self.volume = self._preFL()

        assert len(self.species) == len(self.spring_constants)

        # Read lammps dump file of average positions
        atoms_npt = io.read("output/lammps_preFL.data", format='lammps-data')

        # Reduce to unit cell
        reduced_atoms = reduce_and_avg(atoms_npt, size)

        # Print reduced_atoms for verification
        write('output/reduced_atoms.data', reduced_atoms, format='lammps-data')

        # Check symmetry
        crystal_genome_designation = self._get_crystal_genome_designation_from_atoms_and_verify_unchanged_symmetry(
                reduced_atoms, loose_triclinic_and_monoclinic=False)

        # crystal-structure-npt
        self._add_property_instance_and_common_crystal_genome_keys("crystal-structure-npt", write_temp=True, write_stress=True)
        self._add_file_to_current_property_instance("restart-file","output/lammps_preFL.restart")
        self._add_key_to_current_property_instance("temperature", temperature, "K")
        self._add_key_to_current_property_instance("cell-cauchy-stress", [-pressure, -pressure, -pressure, 0.0, 0.0, 0.0], "atm")
    
        # Rescaling 0K supercell to have equilibrium lattice constant.
        # equilibrium_cell is 3x3 matrix or can also have [len(a), len(b), len(c), angle(b,c), angle(a,c), angle(a,b)]

        self.supercell.set_cell(equilibrium_cell, scale_atoms=True)
        self.supercell.write(
            # os.path.join(
            #     os.path.dirname(os.path.realpath(__file__)),
            #     "output/equilibrium_crystal.dump",
            # ),
            "output/equilibrium_crystal.dump",
            format="lammps-data",
            masses=True,
            atom_style=atom_style,
        )

        # Collect the energies of isolated atoms to subtract from final values
        isolated_atom_energy_list = []

        # List of unique elements (strings) in the "reduced_atoms" ASE atoms object
        element_list = list(dict.fromkeys(reduced_atoms.get_chemical_symbols()))

        # List of stoichiometric coefficients (should line up with "element_list" if everything is in alphabetical order)
        stoichiometry = get_stoich_reduced_list_from_prototype(self.prototype_label)
        
        # Append correct number of each isolated atom energy according to stoichiometric coefficients
        for element,stoich in zip(element_list,stoichiometry):
            isolated_atom_energy_list.append(stoich * get_isolated_energy_per_atom(self.kim_model_name, element)) # appends as many as exist in one formula unit

        # Take the mean of the energies to subtract from the per-atom free energy
        isolated_atom_energy = np.mean(isolated_atom_energy_list)

        # FL computes the free energy at a given pressure and temperature.
        self.templates._write_fl_lammps_templates(
            spring_constants=self.spring_constants
        )
        free_energy_per_atom = self._FL() - isolated_atom_energy

        # Convert to eV/formula (originally in eV/atom)
        # get_stoich_reduced_list_from_prototype returns a list corresponding to the stoichiometry of the prototype label, e.g. "A2B_hP9_152_c_a" -> [2,1]
        num_atoms_in_formula = sum(stoichiometry)
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

        # KIM tries to save some coordinate file, disabling it.
        self.poscar = None

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
        self._add_key_to_current_property_instance(
            "temperature", temperature, "K"
        )
        self._add_key_to_current_property_instance(
            "cell-cauchy-stress", [-pressure, -pressure, -pressure, 0.0, 0.0, 0.0], "atm"
        )

    def _validate_inputs(self):
        if not self.temperature > 0.0:
            raise ValueError("Temperature has to be larger than zero.")

        if not self.pressure >= 0.0:
            raise ValueError("Pressure has to be greater than or equal to zero.")

    def _setup_initial_structure(
        self,
        size: Tuple[int, int, int],
        filename: str = "output/zero_temperature_crystal.dump",
    ) -> Atoms:
        # Copy original atoms so that their information does not get lost when the new atoms are modified.
        atoms_new = self.atoms.copy()

        # Build supercell
        atoms_new = atoms_new.repeat(size)

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

        # Write lammps file.
        # structure_file = os.path.join(
        #     os.path.dirname(os.path.realpath(__file__)), filename
        # )
        # os.makedirs(os.path.dirname(structure_file),exist_ok=True)
        structure_file = filename

        atoms_new.write(structure_file, format="lammps-data", masses=True)
        self.zero_k_structure_path = structure_file

        angles = self.atoms.get_cell().angles()
        self.is_triclinic = (
            all(angle != 90 for angle in angles) and len(set(angles)) == 3
        )

        return atoms_new

    def _preFL(self) -> Tuple[List[float], List[float]]:
        variables = {
            "modelname": self.kim_model_name,
            "temperature": self.temperature,
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
            "temperature": self.temperature,
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
            * self.temperature
            * np.log(HBAR * omega / (KB * self.temperature))
        )  # [eV/atom].

        total_mass = np.sum(natoms * self.concentration * self.mass)

        F_CM = (
            KB
            * self.temperature
            * np.log(
                (natoms / self.volume)
                * (
                    (2 * np.pi * KB * self.temperature)
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
                    return False

                last_line = lines[-2].strip()

                return last_line.startswith("Crystal melted or vaporized")

        except FileNotFoundError:
            print(f"The file {lammps_log} does not exist.")
            return False

    def _modify_accuracies(self, new_accuracies):

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

if __name__ == "__main__":
    model_name = "LJ_Shifted_Bernardes_1958MedCutoff_Ar__MO_126566794224_004"
    subprocess.run(f"kimitems install {model_name}", shell=True, check=True)
    subprocess.run("mkdir -p output", shell=True, check=True)
    test_driver = TestDriver(model_name)
    test_driver(
        bulk("Ar", "fcc", a=5.248),
        size=(3, 3, 3),
        temperature=20.0,
        pressure=0.0,
    )

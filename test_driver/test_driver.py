import os
import subprocess
from typing import List, Tuple
from ase.build import bulk
from ase import Atoms
import numpy as np
from ase.data import atomic_masses, atomic_numbers
from kim_tools import CrystalGenomeTestDriver
import scipy.constants as sc
from .lammps_templates import LammpsTemplates

EV = sc.value("electron volt")
MU = sc.value("atomic mass constant")
HBAR = sc.value("Planck constant in eV/Hz") / (2 * np.pi)
KB = sc.value("Boltzmann constant in eV/K")


class TestDriver(CrystalGenomeTestDriver):
    def _calculate(
        self,
        temperature: float = 20.0,
        pressure: float = 0.0,
        size: Tuple[int, int, int] = (3, 3, 3),
        **kwargs,
    ) -> None:
        """Gibbs free energy of a crystal at constant temperature and pressure using Frenkel-Ladd Hamiltonian integration algorithm. Computed through one equilibrium NPT simulation ('preFL') and one NONequilibrium NVT simulation ('FL').

        Args:
            temperature (float): system temperature in K.
            pressure (float): system pressure in ??.
            size (Tuple[int, int, int]): system size.
        """
        # Check arguments
        self.temperature = temperature
        self.pressure = pressure
        self._validate_inputs()

        # Write initial atomic structure to lammps dump file
        self.supercell = self._setup_initial_structure(size)

        # Write initial template file
        self.templates = LammpsTemplates(root="lammps_templates/")
        self.templates._write_pre_fl_lammps_templates(
            nspecies=len(self.species), is_triclinic=self.is_triclinic
        )

        # preFL computes the equilibrium lattice parameter and spring constants for a given temperature and pressure.
        # TODO: This should probably be replaced with its own test driver, which compute equilibrium lattice constants, and which can handles arbitrary crystal structures. Then we can get spring constants.
        equilibrium_cell, self.spring_constants, self.volume = self._preFL()

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

        # crystal-structure-npt
        self._add_property_instance_and_common_crystal_genome_keys("crystal-structure-npt", write_temp=True, write_stress=True)
        self._add_file_to_current_property_instance("restart-file","output/lammps_preFL.restart")
    
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

        # FL computes the free energy at a given pressure and temperature.
        self.templates._write_fl_lammps_templates(
            spring_constants=self.spring_constants
        )
        free_energy = self._FL()

        # Convert to eV/cell
        free_energy = free_energy * len(self.supercell) / np.prod(size)

        # Print results
        print("####################################")
        print("# Frenkel Ladd Free Energy Results #")
        print("####################################")

        print(f"G_FL = {free_energy:.5f} (eV/cell)")

        # KIM tries to save some coordinate file, disabling it.
        self.poscar = None

        # Write free energy
        self._add_property_instance_and_common_crystal_genome_keys(
            "free-energy", write_stress=True, write_temp=True
        )
        self._add_key_to_current_property_instance(
            "free_energy", free_energy, "eV/cell"
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
            "temperature_damping": 0.1,
            "temperature_seed": np.random.randint(low=100000, high=999999, dtype=int),
            "pressure": self.pressure,
            "pressure_damping": 1.0,
            "timestep": 0.001,  # ps
            "species": " ".join(self.species),
            "output_filename": "output/lammps_preFL.dat",
            "write_restart_filename": "output/lammps_preFL.restart",
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

        F_CM = (
            KB
            * self.temperature
            * np.log(
                (natoms / self.volume)
                * (
                    (2 * np.pi * KB * self.temperature)
                    / (natoms * self.concentration * self.mass * omega**2)
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

        free_energy = np.sum(F_harm) - Work + np.sum(F_CM) + PV_term
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

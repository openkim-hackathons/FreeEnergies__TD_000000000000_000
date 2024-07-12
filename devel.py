import os
import subprocess
from typing import List, Dict, Tuple
from ase.build import bulk
from ase import Atoms
import numpy as np
from ase.data import atomic_masses, atomic_numbers
from kim_tools import CrystalGenomeTestDriver
import scipy.constants as sc

EV = sc.value("electron volt")
MU = sc.value("atomic mass constant")
HBAR = sc.value("Planck constant in eV/Hz") / (2 * np.pi)
KB = sc.value("Boltzmann constant in eV/K")


class FrenkelLaddFreeEnergies(CrystalGenomeTestDriver):

    def _calculate(
        self,
        temperature: float,
        pressure: float,
        size: Tuple[int, int, int],
        **kwargs,
    ) -> None:
        """Gibbs free energy of a crystal at constant temperature and pressure using Frenkel-Ladd Hamiltonian integration algorithm. Computed through one equilibrium NPT simulation ('preFL') and one NONequilibrium NVT simulation ('FL').

        Args:
            temperature (float): system temperature.
            pressure (float): system pressure.
            size (Tuple[int, int, int]): system size.
        """
        # Check arguments

        self.temperature = temperature
        self.pressure = pressure
        self._validate_inputs()

        # Write initial atomic structure to lammps dump file
        supercell = self._setup_initial_structure(size)

        # preFL computes the equilibrium lattice parameter and spring constants for a given temperature and pressure.
        # TODO: This should probably be replaced with its own test driver, which compute equilibrium lattice constants, and which can handles arbitrary crystal structures (right now only works for cubic crystals). Then we can get spring constants.
        equilibrium_cell, self.spring_constants, self.volume = self._preFL()

        # Rescaling 0K supercell to have equilibrium lattice constant.
        # equilibrium_cell is 3x3 matrix.
        # TODO: Divide cell by system size?
        supercell.set_cell(equilibrium_cell, scale_atoms=True)
        supercell.write(
            os.path.join(
                os.path.dirname(os.path.realpath(__file__)),
                "output/equuilibrium_crystal.dump",
            ),
            format="lammps-data",
            masses=True,
        )

        # FL computes the free energy at a given pressure and temperature.
        free_energy = self._FL()

        # Print results
        print("####################################")
        print("# Frenkel Ladd Free Energy Results #")
        print("####################################")

        print(r"$G_{FL} =$" + f" {free_energy:.5f} (eV/atom)")

        # KIM tries to save some coordinate file, disabling it.
        self.poscar = None

        # Write property
        self._add_key_to_current_property_instance(
            "free_energy", free_energy, "eV/atom"
        )

    def _validate_inputs(self):
        if not self.temperature > 0.0:
            raise ValueError("Temperature has to be larger than zero.")

        if not self.pressure > 0.0:
            raise ValueError("Pressure has to be larger than zero.")

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
        self.species = sorted(set(symbols))
        self.mass = [
            atomic_masses[atomic_numbers[element_symbol]]
            for element_symbol in self.species
        ]

        # Write lammps file.
        structure_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)), filename
        )
        atoms_new.write(structure_file, format="lammps-data", masses=True)

        return atoms_new

    def _preFL(self) -> Tuple[List[float], List[float]]:

        variables = {
            "modelname": self.kim_model_name,
            "temperature": self.temperature,
            "temperature_damping": 0.1,
            "temperature_seed": np.random.randint(low=100000, high=999999, dtype=int),
            "pressure": self.pressure,
            "pressure_damping": 1.0,
            "species": " ".join(self.species),
            "output_filename": "output/lammps_preFL.dat",
            "write_restart_filename": "output/lammps_preFL.restart",
        }
        # TODO: Possibly run MPI version of Lammps if available.
        command = (
            "lammps "
            + " ".join(f"-var {key} '{item}'" for key, item in variables.items())
            + " -log output/lammps_preFL.log"
            + " -in lammps_templates/preFL_template.lmp"
        )
        subprocess.run(command, check=True, shell=True)

        # Analyse lammps outputs
        data = np.loadtxt("output/lammps_preFL.log", unpack=True)
        xx, xy, xz, yx, yy, yz, zx, zy, zz, spring_constants = data
        equilibrium_cell = np.array([[xx, xy, xz], [yx, yy, yz], [zx, zy, zz]])

        # TODO: read in volume
        volume = 0

        return equilibrium_cell, spring_constants, volume

    def _FL(
        self,
    ) -> float:

        variables = {
            "modelname": self.kim_model_name,
            "temperature": self.temperature,
            "pressure": self.pressure,
            "species": " ".join(self.species),
            "tswitch": 100000,
            "temperature_damping": 0.01,
            "t_equil": 50000,
            "timestep": 0.001,  # ps
            "spring_constant": self.spring_constant,
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
            + " -in lammps_templates/FL_template.lmp"
        )
        subprocess.run(command, check=True, shell=True)

        return self.compute_free_energy()

    
    def compute_free_energy(self) -> float:
        """Compute free energy via integration of FL path
        """

        Hi_f, Hf_f, lamb_f = np.loadtxt(
            "output/FL_switch1.dat", unpack=True, skiprows=1
        )
        W_forw = np.trapz(Hf_f - Hi_f, lamb_f)

        Hf_b, Hi_b, lamb_b = np.loadtxt("output/FL_switch2.dat.", unpack=True, skiprows=1)
        W_back = np.trapz(Hf_b - Hi_b, 1 - lamb_b)

        Work = (W_forw - W_back) / 2

        omega = (
            np.sqrt(self.spring_constants * EV / (self.mass * MU)) * 1.0e10
        )  # [rad/s].

        F_harm = (
            3 * KB * self.temperature * np.log(HBAR * omega / (KB * self.temperature))
        )  # [eV/atom].

        natoms = len(self.atoms)

        F_CM = (
            KB
            * self.temperature
            * np.log(
                (natoms / self.volume)
                * (2 * np.pi * KB * self.temperature / (natoms * self.mass * omega**2))
                ** (3 / 2)
            )
            / natoms
        )  # correction for fixed center of mass

        free_energy = F_harm - Work + F_CM
        return free_energy


if __name__ == "__main__":
    model_name = "LJ_Shifted_Bernardes_1958MedCutoff_Ar__MO_126566794224_004"
    subprocess.run(f"kimitems install {model_name}", shell=True, check=True)
    test_driver = FrenkelLaddFreeEnergies(model_name)
    test_driver(
        bulk("Ar", "fcc", a=5.248),
        size=(3, 3, 3),
        temperature=20.0,
        pressure=1.0,
    )

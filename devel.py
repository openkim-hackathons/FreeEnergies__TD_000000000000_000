import os
import subprocess
from typing import List, Dict, Tuple
from ase.build import bulk
from ase import Atoms
import numpy as np

from kim_tools import CrystalGenomeTestDriver


class FrenkelLaddFreeEnergies(CrystalGenomeTestDriver):
    def _calculate(
        self,
        temperature: float,
        pressure: float,
        size: Tuple[int, int, int],
        **kwargs,
    ) -> None:
        """
        # TODO: Docstring
        """
        # Check arguments
        if not temperature > 0.0:
            raise ValueError("Temperature has to be larger than zero.")

        if not pressure > 0.0:
            raise ValueError("Pressure has to be larger than zero.")

        # Write initial atomic structure to lammps dump file
        supercell = self._write_initial_structure(size)

        # preFL computes the equilibrium lattice parameter and spring constants for a given temperature and pressure.
        # TODO: This should probably be replaced with its own test driver, which compute equilibrium lattice constants, and which can handles arbitrary crystal structures. Then we can get spring constants.
        equilibrium_cell, spring_constants = self._preFL(
            pressurs=pressure, temperature=temperature
        )

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
        free_energy = self._FL(pressures=pressure, temperature=temperature)

        # Print Result
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

    def _write_initial_structure(
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

        # Write lammps file.
        structure_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)), filename
        )
        atoms_new.write(structure_file, format="lammps-data", masses=True)

        return atoms_new

    def _preFL(
        self, pressure: float, temperature: float
    ) -> Tuple[List[float], List[float]]:

        variables = {
            "modelname": self.kim_model_name,
            "temperature": temperature,
            "temperature_seed": seed,
            "temperature_damping": tdamp,
            "pressure": pressure,
            "pressure_damping": pdamp,
            "timestep": timestep,
            "number_sampling_timesteps": number_sampling_timesteps,
            "species": " ".join(self.species),
            "output_filename": "output/lammps_preFL.dat",
            "write_restart_filename": "output/lammps_preFL_restart.restart",
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
        equilibrium_cell = np.array(
            [[xx, xy, xz], [yx, yy, yz], [zx, zy, zz]]
        )

        return equilibrium_cell, spring_constants

    def _FL(self, pressure: float, temperature: float):

        variables = {
            "modelname": self.kim_model_name,
            "temperature": temperature,
            "temperature_seed": seed,
            "temperature_damping": tdamp,
            "pressure": pressure,
            "pressure_damping": pdamp,
            "timestep": timestep,
            "number_sampling_timesteps": number_sampling_timesteps,
            "species": " ".join(self.species),
            "average_position_filename": "output/average_position_equilibration.dump.*",
            "average_cell_filename": "output/average_cell_equilibration.dump",
            "write_restart_filename": "output/final_configuration_equilibration.restart",
        }
        # TODO: Possibly run MPI version of Lammps if available.
        command = (
            "lammps "
            + " ".join(f"-var {key} '{item}'" for key, item in variables.items())
            + " -log output/lammps_FL.log"
            + " -in lammps_templates/FL_template.lmp"
        )
        subprocess.run(command, check=True, shell=True)

        # TODO: Analyse lammps outputs
        free_energies_vs_pressure_at_temperature = []
        return free_energies_vs_pressure_at_temperature


if __name__ == "__main__":
    model_name = "LJ_Shifted_Bernardes_1958MedCutoff_Ar__MO_126566794224_004"
    subprocess.run(f"kimitems install {model_name}", shell=True, check=True)
    test_driver = FrenkelLaddFreeEnergies(model_name)
    test_driver(
        bulk("Ar", "fcc", a=5.248),
        size=(3, 3, 3),
        temperature=100.0,
        pressure=1.0,
    )

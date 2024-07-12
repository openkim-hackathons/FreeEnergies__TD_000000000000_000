import os
import subprocess
from typing import List, Dict
from ase.build import bulk
import matplotlib.pyplot as plt
import numpy as np

from kim_tools import CrystalGenomeTestDriver


class FreeEnergies(CrystalGenomeTestDriver):
    def _calculate(
        self,
        temperatures: List[float],
        pressures: List[float],
        lammps_args: Dict[str],
        **kwargs,
    ) -> None:
        """ """
        # Check arguments
        if not np.all(temperatures > 0.0):
            raise ValueError("Temperature has to be larger than zero.")

        if not np.all(pressures > 0.0):
            raise ValueError("Pressure has to be larger than zero.")

        # Write initial atomic structure to lammps dump file
        self._write_initial_structure()

        # preFL computes the lattice parameter and spring constants as functions of pressure at the starting temperature (T = 100K).
        
        # TODO: If we notice that this takes too long, maybe use an initial temperature ramp.
        variables = {
            "modelname": self.kim_model_name,
            "temperature": temperature,
            "temperature_seed": seed,
            "temperature_damping": tdamp,
            "pressure": pressure,
            "pressure_damping": pdamp,
            "timestep": timestep,
            "number_sampling_timesteps": number_sampling_timesteps,
            "species": " ".join(species),
            "average_position_filename": "output/average_position_equilibration.dump.*",
            "average_cell_filename": "output/average_cell_equilibration.dump",
            "write_restart_filename": "output/final_configuration_equilibration.restart",
        }
        # TODO: Possibly run MPI version of Lammps if available.
        command = (
            "lammps "
            + " ".join(f"-var {key} '{item}'" for key, item in variables.items())
            + " -log output/lammps_equilibration.log"
            + " -in npt_equilibration.lammps"
        )
        subprocess.run(command, check=True, shell=True)

        # Analyse equilibration run.

        # Run first NPT simulation at higher temperature.

        # Print Result
        print("####################################")
        print("# Frenkel Ladd Free Energy Results #")
        print("####################################")

        print()

        # I have to do this or KIM tries to save some coordinate file
        self.poscar = None

        # Write property

    def _write_initial_structure(
        self, filename: str = "output/zero_temperature_crystal.lmp"
    ):

        # Copy original atoms so that their information does not get lost when the new atoms are modified.
        atoms_new = self.atoms.copy()

        # This is how ASE obtains the species that are written to the initial configuration.
        # These species are passed to kim interactions.
        # See https://wiki.fysik.dtu.dk/ase/_modules/ase/io/lammpsdata.html#write_lammps_data
        symbols = atoms_new.get_chemical_symbols()
        species = sorted(set(symbols))

        # Write lammps file.
        structure_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)), filename
        )
        atoms_new.write(structure_file, format="lammps-data", masses=True)

        return atoms_new


if __name__ == "__main__":
    model_name = "LJ_Shifted_Bernardes_1958MedCutoff_Ar__MO_126566794224_004"
    subprocess.run(f"kimitems install {model_name}", shell=True, check=True)
    test_driver = FreeEnergies(model_name)
    test_driver(
        bulk("Ar", "fcc", a=5.248).repeat((3, 3, 3)),
        temperature=10.0,
        pressure=1.0,
    )

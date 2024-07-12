from math import ceil
import os
import random
import re
import subprocess
from typing import Iterable, List, Optional, Tuple, Dict
from ase import Atoms
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
        # TODO: Check all arguments.

        # Copy original atoms so that their information does not get lost when the new atoms are modified.
        atoms_new = self.atoms.copy()

        # This is how ASE obtains the species that are written to the initial configuration.
        # These species are passed to kim interactions.
        # See https://wiki.fysik.dtu.dk/ase/_modules/ase/io/lammpsdata.html#write_lammps_data
        symbols = atoms_new.get_chemical_symbols()
        species = sorted(set(symbols))

        # Write lammps file.
        TDdirectory = os.path.dirname(os.path.realpath(__file__))
        structure_file = os.path.join(
            TDdirectory, "output/zero_temperature_crystal.lmp"
        )
        atoms_new.write(structure_file, format="lammps-data", masses=True)


        # TODO: Move damping factors to argument.
        pdamp = timestep * 100.0
        tdamp = timestep * 1000.0

        # Run NPT simulation for equilibration.
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
        print("# NPT Phonon Heat Capacity Results #")
        print("####################################")
        print(f"C_p:\t{c}")
        print(f"C_p Error:\t{c_err}")

        # I have to do this or KIM tries to save some coordinate file
        self.poscar = None

        # Write property
        self._add_property_instance_and_common_crystal_genome_keys(
            "heat-capacity-phonon-npt", write_stress=True, write_temp=True
        )  # last two default to False
        self._add_key_to_current_property_instance(
            "constant_pressure_heat_capacity", c, "eV/Kelvin"
        )
        self._add_key_to_current_property_instance(
            "constant_pressure_heat_capacity_err", c_err, "eV/Kelvin"
        )
        self._add_key_to_current_property_instance(
            "pressure", variables["pressure"], "bars"
        )



if __name__ == "__main__":
    model_name = "LJ_Shifted_Bernardes_1958MedCutoff_Ar__MO_126566794224_004"
    subprocess.run(f"kimitems install {model_name}", shell=True, check=True)
    test_driver = FreeEnergies(model_name)
    test_driver(
        bulk("Ar", "fcc", a=5.248).repeat((3, 3, 3)),
        temperature=10.0,
        pressure=1.0,
    )

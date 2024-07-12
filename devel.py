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
        temperatures: List[float],
        pressures: List[float],
        size: Tuple[int, int, int],
        **kwargs,
    ) -> None:
        """
        # TODO: Docstring
        """
        # Check arguments
        if not np.all(np.array(temperatures) > 0.0):
            raise ValueError("Temperature has to be larger than zero.")

        if not np.all(np.array(pressures) > 0.0):
            raise ValueError("Pressure has to be larger than zero.")

        # Write initial atomic structure to lammps dump file
        supercell = self._write_initial_structure(size)

        # preFL computes the lattice parameter and spring constants as functions of pressure at the starting temperature.
        # TODO: Does it work for anythign else than cubic structures? 
        equilibrium_lattice_parameters, spring_constants = self._preFL(
            pressures=pressures, temperature=np.min(temperatures)
        )

        # TODO: Need to scale self.atoms by the equilibrium lattice constant, write it to a dump file and load it in FL
        breakpoint()

        # FL computes the free energy as a function of pressure at the starting temperature (list of free energies vs P at T = starting temperature).
        free_energies_vs_pressure_at_temperature = self._FL(
            pressures=pressures, temperature=np.min(temperatures)
        )

        # RS computes the free energy at each pressure over a mesh of temperatures.
        free_energies_vs_pressure_vs_temperature = self._RS(
            pressures=pressures, temperature=np.min(temperatures)
        )

        # Print Result
        print("####################################")
        print("# Frenkel Ladd Free Energy Results #")
        print("####################################")

        print(
            r"$G_{FL} =$" + f" {free_energies_vs_pressure_vs_temperature:.5f} (eV/atom)"
        )

        # I have to do this or KIM tries to save some coordinate file
        self.poscar = None

        # Write property
        # TODO: Write them using kim utils helps for writting kim properties.
        self._add_key_to_current_property_instance("constant_pressure_free_energy", XX, "eV/atom")

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
        self, pressures: List[float], temperature: float
    ) -> Tuple[List[float], List[float]]:

        for pressure in pressures:
        

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
                + " -log output/lammps_preFL.log"
                + " -in lammps_templates/preFL_template.lmp"
            )
            subprocess.run(command, check=True, shell=True)

        # Analyse lammps outputs
        data = np.loadtxt("output/lammps_preFL.log", unpack=True)
        _, equilibrium_lattice_parameters, spring_constants = (
            data[:, 0],
            data[:, 1],
            data[:, 2],
        )

        return equilibrium_lattice_parameters, spring_constants

    def _FL(self, pressures: List[float], temperature: float) -> List[float]:

        for pressure in pressures:

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

    def _RS(
        self,
        pressures: List[float],
    ):
        for pressure in pressures:

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
                + " -log output/lammps_RS.log"
                + " -in lammps_templates/RS_template.lmp"
            )
            subprocess.run(command, check=True, shell=True)

        # TODO: Analyse lammps outputs
        free_energies_vs_pressure_vs_temperature = []
        return free_energies_vs_pressure_vs_temperature


if __name__ == "__main__":
    model_name = "LJ_Shifted_Bernardes_1958MedCutoff_Ar__MO_126566794224_004"
    subprocess.run(f"kimitems install {model_name}", shell=True, check=True)
    test_driver = FrenkelLaddFreeEnergies(model_name)
    test_driver(
        bulk("Ar", "fcc", a=5.248),
        size=(3, 3, 3),
        temperatures=[10.0, 20.0, 30.0],
        pressures=[1.0, 2.0, 3.0],
    )

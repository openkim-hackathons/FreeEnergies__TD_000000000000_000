from math import ceil
import os
import re
import subprocess
from typing import Iterable, List, Tuple
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt

def run_lammps(modelname: str, temperature_K: float, pressure_bar: float, timestep_ps: float,
               number_sampling_timesteps: int, species: List[str],
               msd_threshold_angstrom_squared_per_sampling_timesteps: float, number_msd_timesteps: int, number_avePOS_timesteps: int,
               rlc_N_every: int, lammps_command: str, random_seed: int, output_dir: str, test_driver_dir: str) -> Tuple[str, str, list]:
    
    pdamp = timestep_ps * 1000.0
    tdamp = timestep_ps * 100.0

    log_filename = f"{output_dir}/free_energy.log"
    restart_filename = f"{output_dir}/free_energy.restart"

    variables = {
        "modelname": modelname,
        "temperature": temperature_K,
        "temperature_damping": tdamp, # picoseconds
        "velocity_seed": random_seed,
        "langevin_seed": int(7*(random_seed)/5), # if random_seed is 101010 (default), this is 141414
        "pressure": pressure_bar,
        "t_switch": 10000,
        "t_equil": 10000,
        "number_msd_timesteps": number_msd_timesteps,
        "number_avePOS_timesteps": number_avePOS_timesteps,
        "rlc_N_every": rlc_N_every,
        "pressure_damping": pdamp, # picoseconds
        "msd_threshold": 0.1, # Angstrom^2 per 100 timesteps
        "timestep": timestep_ps,  # picoseconds
        "species": " ".join(species),
        "output_dir": f"{output_dir}",
        "output_filename": f"{output_dir}/free_energy.dat",
        "write_restart_filename": f"{restart_filename}",
        "write_data_filename": f"{output_dir}/free_energy.data",
        "write_dump_filename": f"{output_dir}/free_energy.dump",
        "zero_temperature_crystal": f"{output_dir}/zero_temperature_crystal.data",
        "melted_crystal_output": f"{output_dir}/melted_crystal.dump",
        "switch1_output_file": f"{output_dir}/FL_switch1.dat",
        "switch2_output_file": f"{output_dir}/FL_switch2.dat",
        "run_length_control": f"{test_driver_dir}/run_length_control.py",
    }
    
    # Construct base LAMMPS command
    command = (
        f"{lammps_command} "
        + " ".join(f"-var {key} '{item}'" for key, item in variables.items())
        + f" -log {log_filename}"
        + f" -in {output_dir}/in.lammps")

    # Run LAMMPS
    subprocess.run(command, check=True, shell=True)

    # Read spring constants from file written by lammps
    spring_constants = np.zeros(len(species))
    for i in range(len(species)):
        spring_constants[i] = np.loadtxt(f"{output_dir}/k{i}.dat")
    
    volume  = np.loadtxt(f"{output_dir}/volume.dat")
    
    return log_filename, restart_filename, spring_constants, volume
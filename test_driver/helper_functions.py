import subprocess
from typing import List, Tuple
import numpy as np
from pathlib import Path
from enum import Enum

class LammpsStatus(Enum):
    """Status of LAMMPS simulation completion."""
    SUCCESS = "success"
    MELTED = "melted"
    NOT_FOUND = "not_found"
    INCOMPLETE = "incomplete"

def run_lammps(modelname: str, temperature_K: float, pressure_bar: float, timestep_ps: float, fl_switch_timesteps: int, fl_equil_timesteps: int,
               species: List[str], msd_threshold: float, msd_timesteps: int, thermo_sampling_period: int, ave_pos_timesteps: int,
               rlc_n_every: int, lammps_command: str, random_seed: int, output_dir: str) -> Tuple[str, str, list, float]:
    
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
        "t_switch": fl_switch_timesteps,
        "t_equil": fl_equil_timesteps,
        "msd_timesteps": msd_timesteps,
        "ave_pos_timesteps": ave_pos_timesteps,
        "rlc_n_every": rlc_n_every,
        "pressure_damping": pdamp, # picoseconds
        "msd_threshold": msd_threshold, # Angstrom^2 per "msd_sampling_period" timesteps
        "thermo_sampling_period": thermo_sampling_period,
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
        "run_length_control": f"{output_dir}/run_length_control.py",
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

    # Verify simulation completed successfully
    lammps_status = _check_lammps_completion(lammps_log=f"{output_dir}/free_energy.log")
    
    return log_filename, restart_filename, spring_constants, volume, lammps_status

def _check_lammps_completion(lammps_log: Path) -> LammpsStatus:
        """Check if LAMMPS simulation completed successfully.
        
        Args:
            lammps_log: Path to the LAMMPS log file.
            
        Returns:
            LammpsStatus indicating the outcome of the simulation.
        """
        try:
            with open(lammps_log, "r", encoding="utf-8") as file:
                lines = file.readlines()

                if not lines:
                    return LammpsStatus.INCOMPLETE
                
                # Check second-to-last line for melting indicator
                if len(lines) >= 2 and lines[-2].strip().startswith("Crystal melted or vaporized"):
                    return LammpsStatus.MELTED

                # Check last line for successful completion
                if lines[-1].strip().startswith("Total wall time:"):
                    return LammpsStatus.SUCCESS
                    
                return LammpsStatus.INCOMPLETE

        except FileNotFoundError:
            return LammpsStatus.NOT_FOUND

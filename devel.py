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


class LammpsTemplates:
    def __init__(self, root):
        self.root = root
        os.makedirs(self.root, exist_ok=True)
        self.pre_fl = """
        kim init ${modelname} metal unit_conversion_mode

        # periodic boundary conditions along all three dimensions
        boundary p p p

        # Read crystal with 0K lattice parameter.
        read_data output/zero_temperature_crystal.dump

        # Convert box and all atomic positions to the correct units.
        #change_box all x scale ${_u_distance} &
        #               y scale ${_u_distance} &
        #               z scale ${_u_distance} &
        #               xy final $(xy*v__u_distance) &
        #               xz final $(xz*v__u_distance) &
        #               yz final $(yz*v__u_distance) &
        #               remap

        # Interatomic potential and neighbor settings
        kim           interactions ${species}

        # set the time step to 0.001 picoseconds
        variable      timestep_converted equal ${timestep}*${_u_time}
        timestep      ${timestep_converted}

        # Leaving pressure variables just in case we need to compute lattice parameters
        variable      temp_converted equal ${temperature}*${_u_temperature}
        variable      Tdamp_converted equal ${temperature_damping}*${_u_time}
        variable      press_converted equal ${pressure}*${_u_pressure}
        variable      Pdamp_converted equal ${pressure_damping}*${_u_time}

        # Initialize velocities.
        velocity      all create ${temp_converted} ${temperature_seed}

        # Set thermodynamic ensemble
        # if need to compute new box matrix, use uncomment second line
        fix           ensemble all nvt temp ${temp_converted} ${temp_converted} ${Tdamp_converted}
        #fix          ensemble all npt temp ${temp_converted} ${temp_converted} ${Tdamp_converted} tri ${press_converted} ${press_converted} ${Pdamp_converted}

        # Initialize measurement of box vectors
        #thermo_style custom avecx avecy avecz bvecx bvecy bvecz cvecx cvecy cvecz temp press vol etotal step
        thermo_style custom lx ly lz xy yz xz temp press vol etotal step
        thermo 1000

        # compute box information
        #variable      avecx_metal equal avecx/${_u_distance}
        #variable      avecy_metal equal avecy/${_u_distance}
        #variable      avecz_metal equal avecz/${_u_distance}

        #variable      bvecx_metal equal bvecx/${_u_distance}
        #variable      bvecy_metal equal bvecy/${_u_distance}
        #variable      bvecz_metal equal bvecz/${_u_distance}

        #variable      cvecx_metal equal cvecx/${_u_distance}
        #variable      cvecy_metal equal cvecy/${_u_distance}
        #variable      cvecz_metal equal cvecz/${_u_distance}

        variable       lx_metal equal lx/${_u_distance}
        variable       ly_metal equal ly/${_u_distance}
        variable       lz_metal equal lz/${_u_distance}

        # compute mean squared displacement
        compute       MSD all msd com yes average yes

        # Temperature may be off because of rigid bodies or SHAKE constraints. See https://docs.lammps.org/velocity.html
        run 0
        velocity all scale $(v_temperature*v__u_temperature)

        # Set up convergence check with kim-convergence.
        # Alternative is to set a safe general equilibration time
        # More important for expensive potentials
        #python run_length_control input 6 SELF 1 variable msd_metal format pissss file run_length_control_preFL.py

        # Run until converged.
        #python run_length_control invoke

        # Equilibration
        run 10000

        # compute averages of above variables
        fix           AVG all ave/time 100 100 10000 v_lx_metal v_ly_metal v_lz_metal c_MSD[4] ave running
                                                    #v_bvecx_metal v_bvecy_metal v_bvecz_metal &
                                                    #v_cvecx_metal v_cvecy_metal v_cvecz_metal &
                                                    #c_MSD[4] ave running

        # Run steps in the converged regime.
        run 10000

        # compute spring constant
        variable      kB equal 8.617333262e-5*(v__u_energy/v__u_temperature) # eV/K unless converted
        variable      MSD equal f_AVG[4]*(v__u_distance)^2
        variable      spring_constant equal $(3*v_kB*v_temp_converted/(v_MSD)^2)

        # Write final averages and spring constant
        #print "# xx | xy | xz | yx | yy | yz | zx | zy | zz | spring constant [eV/Ang^2]" file ${output_filename}
        #print "$(f_AVG[1]) $(f_AVG[2]) $(f_AVG[3]) $(f_AVG[4]) $(f_AVG[5]) $(f_AVG[6]) $(f_AVG[7]) $(f_AVG[8]) $(f_AVG[9]) ${spring_constant}" screen no append ${output_filename}
        print "# lx | ly | lz [Ang] | vol [Ang^3] | spring constant [eV/Ang^2]" file ${output_filename}
        print "$(f_AVG[1]) $(f_AVG[2]) $(f_AVG[3]) $(vol) ${spring_constant}" screen no append ${output_filename}

        # Reset.
        unfix ensemble
        unfix AVG
        #unfix cr_fix  # From run_length_control.py
        reset_timestep 0

        # Write the initial starting file for a true simulation.
        write_restart ${write_restart_filename}



        """

        self.fl = """

        kim init ${modelname} metal unit_conversion_mode

        # periodic boundary conditions along all three dimensions
        boundary p p p

        # Read crystal equilibrium crystal.
        read_data output/equilibrium_crystal.dump


        # Convert box and all atomic positions to the correct units.
        #change_box all x scale ${_u_distance} &
        #               y scale ${_u_distance} &
        #               z scale ${_u_distance} &
        #               xy final $(xy*v__u_distance) &
        #               xz final $(xz*v__u_distance) &
        #               yz final $(yz*v__u_distance) &
        #               remap

        # Interatomic potential and neighbor settings
        kim           interactions ${species}

        # set the time step to 0.001 picoseconds
        variable timestep_converted equal ${timestep}*${_u_time}
        timestep ${timestep_converted}

        variable temp_converted equal ${temperature}*${_u_temperature}
        variable Tdamp_converted equal ${temperature_damping}*${_u_time}

        # create initial velocities consistent with the chosen temperature
        velocity      all create ${temp_converted} ${temperature_seed}

        # Frenkel-Ladd fix
        {fix_springs}

        # set NVT ensemble for all atoms
        fix           ensemble all nvt temp ${temp_converted} ${temp_converted} ${Tdamp_converted}
        compute       cl all temp/com
        fix_modify    ensemble temp cl

        # print data to logfile every 1000 timesteps
        variable      etotal_metal equal etotal/${_u_energy}
        variable      pe_metal equal pe/${_u_energy}
        variable      T_metal equal temp/${_u_temperature}
        variable      P_metal equal press/${_u_pressure}

        thermo_style  custom step etotal v_etotal_metal pe v_pe_metal &
                    temp v_T_metal press v_P_metal
        thermo        1000

        # run equil 1
        run ${t_equil}

        # run switch 1
        {record_template}
                    screen no file ${switch1_output_file}
        run ${t_switch}
        unfix record

        # run equil 2
        run ${t_equil}

        # run switch 2
        {record_template}
                    screen no file ${switch2_output_file}
        run ${t_switch}
        unfix record

        print "LAMMPS calculation completed"
        quit 0
        """

    def _add_fl_fix_for_multicomponent(self, spring_constants: List[float]):

        fix_springs_template = """
        group {group} type {group}
        fix           {fix_name} {group} ti/spring {spring_constant} {t_switch} {t_equil} function {function}
        """

        fix_entries = [
            {
                "fix_name": f"FL{i}",
                "group": f"{i+1}",
                "spring_constant": spring_constants[i],
                "t_switch": "${t_switch}",
                "t_equil": "${t_equil}",
                "function": 2,
            }
            for i in range(len(spring_constants))
        ]

        fix_springs = "".join(
            [fix_springs_template.format(**entry) for entry in fix_entries]
        )

        self.fl = self.fl.replace("{fix_springs}", fix_springs)

        record_template = """
        fix           record all print 1 "$(pe/atoms) {data}" &
                     title "{title}" &
        """

        terms = [f"f_FL{i}" for i in range(len(self.spring_constants))]
        terms = "+".join(terms)
        final_sum = f"$(({terms})/atoms) $(f_FL0[1])"

        fix_entries = [
            {
                "fix_name": "record",
                "group": "all",
                "data": final_sum,
                "title": "# PE_potential [eV/atom] | PE_FL [eV/atom] | lambda",
            }
        ]

        record_template = "".join(
            [record_template.format(**entry) for entry in fix_entries]
        )

        self.fl = self.fl.replace("{record_template}", record_template)

    def _add_msd_fix_for_multicomponent(self, nspecies: int):

        compute_msd_template = """
        group {group} type {group}
        compute           {compute_name} {group} msd com yes average yes
        """

        compute_entries = [
            {
                "compute_name": f"MSD{i}",
                "group": f"{i+1}",
            }
            for i in range(nspecies)
        ]

        compute_msd = "".join(
            [compute_msd_template.format(**entry) for entry in compute_entries]
        )

        self.pre_fl = self.pre_fl.replace("{compute_msd}", compute_msd)

        avg_template = """
        fix          AVG all ave/time 100 100 10000 v_lx_metal v_ly_metal v_lz_metal {msd_data} ave running"
        """

        terms = [f"c_MSD{i}[4]" for i in range(nspecies)]
        terms = " ".join(terms)
        msd_data = f"{terms}"

        fix_entries = [
            {
                "fix_name": "AVG",
                "group": "all",
                "msd_data": msd_data,
            }
        ]

        avg_template = "".join([avg_template.format(**entry) for entry in fix_entries])

        self.pre_fl = self.pre_fl.replace("{avg_template}", avg_template)

        k_template = """
        variable      {variable_name_1} equal f_AVG[{avg_i}]*(v__u_distance)^2
        variable      {variable_name_2} equal $(3*v_kB*v_temp_converted/(v_{variable_name_1})^2)
        """

        k_entries = [
            {
                "variable_name_1": f"MSD{i}",
                "variable_name_2": f"spring_constant_{i}",
                "avg_i": f"{i+4}",
                "group": f"{i+1}",
            }
            for i in range(nspecies)
        ]

        k_lines = "".join(
            [k_template.format(**entry) for entry in k_entries]
        )

        self.pre_fl = self.pre_fl.replace("{compute_msd}", k_lines)

    def _write_pre_fl_lammps_templates(self, nspecies: int):
        self._add_msd_fix_for_multicomponent(nspecies)
        open(self.root + "preFL_template.lmp", "w").write(self.pre_fl)

    def _write_fl_lammps_templates(self, spring_constants: List[float]):
        self._add_fl_fix_for_multicomponent(spring_constants)
        open(self.root + "FL_template.lmp", "w").write(self.fl)


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

        # Write initial template file
        self.templates = LammpsTemplates(root="lammps_templates/")
        self.templates._write_pre_fl_lammps_templates(nspecies=len(self.species))

        # preFL computes the equilibrium lattice parameter and spring constants for a given temperature and pressure.
        # TODO: This should probably be replaced with its own test driver, which compute equilibrium lattice constants, and which can handles arbitrary crystal structures (right now only works for cubic crystals). Then we can get spring constants.
        equilibrium_cell, self.spring_constants, self.volume = self._preFL()
        assert len(self.species) == len(self.spring_constants)

        # Rescaling 0K supercell to have equilibrium lattice constant.
        # equilibrium_cell is 3x3 matrix or can also have [len(a), len(b), len(c), angle(b,c), angle(a,c), angle(a,b)]
        # TODO: Divide cell by system size?
        supercell.set_cell(equilibrium_cell, scale_atoms=True)
        supercell.write(
            os.path.join(
                os.path.dirname(os.path.realpath(__file__)),
                "output/equilibrium_crystal.dump",
            ),
            format="lammps-data",
            masses=True,
        )

        # FL computes the free energy at a given pressure and temperature.
        self.templates._write_fl_lammps_templates(spring_constans=self.spring_constants)
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

    ####################################

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
        subprocess.run(command, check=True, shell=True)

        # Analyse lammps outputs
        data = np.loadtxt("output/lammps_preFL.dat", unpack=True)
        # xx, xy, xz, yx, yy, yz, zx, zy, zz, spring_constants = data
        lx, ly, lz, volume, spring_constants = data
        equilibrium_cell = np.array([[lx, 0, 0], [0, ly, 0], [0, 0, lz]])

        return equilibrium_cell, spring_constants, volume

    def _FL(
        self,
    ) -> float:

        self._add_fl_fix_for_multicomponent()
        variables = {
            "modelname": self.kim_model_name,
            "temperature": self.temperature,
            "pressure": self.pressure,
            "species": " ".join(self.species),
            "t_switch": 10000,
            "temperature_damping": 0.01,
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
            "output/FL_switch2.dat.", unpack=True, skiprows=1
        )
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

        free_energy = np.sum(F_harm) - Work + F_CM
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

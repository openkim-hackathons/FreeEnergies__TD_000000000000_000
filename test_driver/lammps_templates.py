import os
from typing import List


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

        # Increase maximum number of neighbors allowed per atom (for potentials that support many species)
        #neigh_modify one 100

        # Change to triclinic box.
        change_box all triclinic

        # Interatomic potential and neighbor settings
        kim           interactions ${species}

        # set the time step to 0.001 picoseconds
        variable      timestep_converted equal ${timestep}*${_u_time}
        timestep      ${timestep_converted}

        # Convert temperature and pressure units
        variable      temp_converted equal ${temperature}*${_u_temperature}
        variable      Tdamp_converted equal ${temperature_damping}*${_u_time}
        variable      press_converted equal ${pressure}*${_u_pressure}
        variable      Pdamp_converted equal ${pressure_damping}*${_u_time}

        # Initialize velocities.
        velocity      all create ${temp_converted} ${temperature_seed}

        # Run NPT ensemble (barostat type depends on box type)
        fix          ensemble all npt temp ${temp_converted} ${temp_converted} ${Tdamp_converted} tri ${press_converted} ${press_converted} ${Pdamp_converted}
        compute      cl all temp/com
        fix_modify   ensemble temp cl

        # Temperature may be off because of rigid bodies or SHAKE constraints. See https://docs.lammps.org/velocity.html
        run 0
        velocity all scale ${temp_converted}

        # compute box information
        variable       lx_metal equal lx/${_u_distance}
        variable       ly_metal equal ly/${_u_distance}
        variable       lz_metal equal lz/${_u_distance}
        variable       xy_metal equal xy/${_u_distance}
        variable       yz_metal equal yz/${_u_distance}
        variable       xz_metal equal xz/${_u_distance}
        variable       vol_metal equal vol/(${_u_distance}^3)

        # Short run to equilibrate MSD
        compute msd all msd com yes
        thermo_style custom lx ly lz xy yz xz temp press vol etotal c_msd[4] step
        thermo 1000
        run 5000
        reset_timestep 0

        # Compute slope of mean squared displacement to detect diffusion
        fix msd_vector all vector 100 c_msd[4]
        variable msd_slope equal slope(f_msd_vector)

        # Thermodynamic output
        thermo_style custom lx ly lz xy yz xz temp press vol etotal c_msd[4] v_msd_slope step
        thermo 1000
        
        # Before kim-convergence, perform a short run and decide whether or not to quit
        run 5000
        if "${msd_slope} > 1e-2" then "write_dump all atom output/melted_crystal.dump" &
                                      "print 'Crystal melted or vaporized'" &
                                      "quit"
        unfix msd_vector
        reset_timestep 0
        thermo_style custom lx ly lz xy yz xz temp press vol etotal step

        # Set up convergence check with kim-convergence.
        python run_length_control input 16 SELF 1 variable vol_metal variable lx_metal variable ly_metal variable lz_metal variable xy_metal variable xz_metal variable yz_metal format pissssssssssssss file test_driver/run_length_control_preFL.py

        # Run until converged (minimum runtime 30000 steps)
        python run_length_control invoke

        unfix cr_fix # From run_length_control.py
        reset_timestep 0

        # Compute mean squared displacement
        #set group all image 0 0 0
        {compute_msd}

        # New set of values to print to log file
        {thermo_template}

        # Define variables for fix ave/time
        variable N_every equal 100 # sample msd at intervals of this many steps
        variable run_time equal 50000 # can be an input variable
        variable N_repeat equal v_run_time/(2*v_N_every)
        variable N_start equal v_run_time/2

        # compute averages of box vectors and msd
        {avg_template}

        # Compute unwrapped atom positions
        compute unwrapped all property/atom xu yu zu

        variable xu atom c_unwrapped[1]
        variable yu atom c_unwrapped[2]
        variable zu atom c_unwrapped[3]

        # Average unwrapped atom positions
        fix avePos all ave/atom ${N_every} $(v_run_time/v_N_every) ${run_time} v_xu v_yu v_zu

        # Run steps in the converged regime.
        run ${run_time}

        # compute spring constant
        variable      kB equal (8.617333262e-5)*(v__u_energy/v__u_temperature) # eV/K unless converted
        {k_lines}

        # Write final averages and spring constant
        print "# lx | ly | lz [Ang] | vol [Ang^3] | spring constant [eV/Ang^2]" file ${output_filename}
        {print_template}

        # Write the initial starting file for a true simulation.
        write_restart ${write_restart_filename}

        # Create per-atom variables for averaged positions
        variable ave_x atom f_avePos[1]
        variable ave_y atom f_avePos[2]
        variable ave_z atom f_avePos[3]

        # Set box dimensions to time-averaged dimensions
        change_box all x scale $(f_AVG[1]/lx) y scale $(f_AVG[2]/ly) z scale $(f_AVG[3]/lz) xy final $(f_AVG[4]) yz final $(f_AVG[5]) xz final $(f_AVG[6]) remap

        # Set atom positions to time-averaged positions
        set group all x v_ave_x y v_ave_y z v_ave_z

        # Reset.
        unfix ensemble
        unfix AVG
        #unfix ave_x
        #unfix ave_y
        #unfix ave_z
        unfix avePos
        #unfix cr_fix  # From run_length_control.py
        reset_timestep 0

        # Write data file of averaged positions and box dimensions
        write_data ${write_data_filename}
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

        # set NVE ensemble for all atoms
        fix           ensemble all nve

        # Frenkel-Ladd fix
        {fix_springs}

        # set langevin thermostat (NVE-Langevin samples spring system more accurately than NVT)
        fix           thermostat all langevin ${temp_converted} ${temp_converted} ${Tdamp_converted} ${temperature_seed} zero yes

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
                     title "{title}" &"""

        terms = [f"f_FL{i}" for i in range(len(spring_constants))]
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

    def _add_msd_fix_for_multicomponent(self, nspecies: int, is_triclinic: bool):
        compute_msd_template = """
        group {group} type {group}
        compute           {compute_name} {group} msd com yes
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

        thermo_template = """
        thermo_style custom lx ly lz xy yz xz temp press vol etotal {msd_data} step
        """

        terms = [f"c_MSD{i}[4]" for i in range(nspecies)]
        terms = " ".join(terms)
        msd_data = f"{terms}"

        fix_entries = [
            {
                "msd_data": msd_data,
            }
        ]

        thermo_template = "".join(
            [thermo_template.format(**entry) for entry in fix_entries]
        )

        self.pre_fl = self.pre_fl.replace("{thermo_template}", thermo_template)

        avg_template = """
        fix          AVG all ave/time ${{N_every}} ${{N_repeat}} ${{run_time}} v_lx_metal v_ly_metal v_lz_metal v_xy_metal v_yz_metal v_xz_metal {msd_data} ave running start ${{N_start}}
        """

        fix_entries = [
            {
                "msd_data": msd_data,
            }
        ]

        avg_template = "".join([avg_template.format(**entry) for entry in fix_entries])

        self.pre_fl = self.pre_fl.replace("{avg_template}", avg_template)

        k_template = """
        variable      {variable_name_1} equal f_AVG[{avg_i}]*(v__u_distance)^2
        variable      {variable_name_2} equal $(3*v_kB*v_temp_converted/(v_{variable_name_1}))
        """

        k_entries = [
            {
                "variable_name_1": f"MSD{i}",
                "variable_name_2": f"spring_constant_{i}",
                "avg_i": f"{i+7}",
            }
            for i in range(nspecies)
        ]

        k_lines = "".join([k_template.format(**entry) for entry in k_entries])

        self.pre_fl = self.pre_fl.replace("{k_lines}", k_lines)

        print_template = """
        print "$(f_AVG[1]) $(f_AVG[2]) $(f_AVG[3]) $(f_AVG[4]) $(f_AVG[5]) $(f_AVG[6]) $(vol) {spring_constants}" screen no append ${{output_filename}}
        """

        terms = [f"$(v_spring_constant_{i})" for i in range(nspecies)]
        terms = " ".join(terms)
        spring_constants = f"{terms}"

        print_entries = [
            {
                "spring_constants": spring_constants,
            }
        ]

        print_template = "".join(
            [print_template.format(**entry) for entry in print_entries]
        )

        self.pre_fl = self.pre_fl.replace("{print_template}", print_template)

        # determine barostat type based on self.is_triclinic
        '''
        self.pre_fl = (
            self.pre_fl.replace("{cell_type}", "tri")
            if is_triclinic
            else self.pre_fl.replace("{cell_type}", "aniso")
        )
        '''
    
    def _write_pre_fl_lammps_templates(self, nspecies: int, is_triclinic: bool):
        
        self._add_msd_fix_for_multicomponent(nspecies, is_triclinic)
        open(self.root + "preFL_template.lmp", "w").write(self.pre_fl)

    def _write_fl_lammps_templates(self, spring_constants: List[float]):
        self._add_fl_fix_for_multicomponent(spring_constants)
        open(self.root + "FL_template.lmp", "w").write(self.fl)

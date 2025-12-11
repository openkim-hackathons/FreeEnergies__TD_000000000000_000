import os
from typing import List


class LammpsTemplate:
    def __init__(self, root):
        self.root = root
        os.makedirs(self.root, exist_ok=True)
        
        self.free_energy = """
        kim init ${modelname} metal unit_conversion_mode

        # periodic boundary conditions along all three dimensions
        boundary p p p

        # Read crystal with 0K lattice parameter.
        read_data ${zero_temperature_crystal}

        # Increase maximum number of neighbors allowed per atom (for potentials that support many species)
        neigh_modify delay 0 every 1 check yes one 4000

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
        velocity      all create ${temp_converted} ${velocity_seed}

        # Run NPT ensemble (barostat type depends on box type)
        fix          ensemble all npt temp ${temp_converted} ${temp_converted} ${Tdamp_converted} tri ${press_converted} ${press_converted} ${Pdamp_converted} flip no

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
        variable       temp_metal equal temp/${_u_temperature}

        thermo_style custom lx ly lz xy yz xz temp press vol etotal step
        thermo 1000

        #==================================== Check melting with MSD ==================================#

        # Compute mean squared displacement
        compute msd all msd com yes

        # Perform a short run to equilibrate MSD
        run 5000
        
        # Calculate slope of MSD to detect diffusion
        fix msd_vector all vector 100 c_msd[4]
        variable msd_slope equal slope(f_msd_vector)

        # Thermodynamic output
        thermo_style custom lx ly lz xy yz xz temp press vol etotal c_msd[4] v_msd_slope step
        thermo 1000
        
        # Perform a short run and decide whether or not to quit
        run 5000
        if "${msd_slope} > $(v_msd_threshold*(v__u_distance^2)/v__u_time)" then "write_dump all atom ${melted_crystal_output}" &
                                      "print 'Crystal melted or vaporized'" &
                                      "quit"
        unfix msd_vector

        thermo_style custom lx ly lz xy yz xz temp press vol etotal step
        thermo 1000

        print '#================================== Kim-convergence ==================================#'

        # Set up convergence check with kim-convergence.
        python run_length_control input 18 SELF 1 variable vol_metal variable temp_metal variable lx_metal variable ly_metal variable lz_metal variable xy_metal variable xz_metal variable yz_metal format pissssssssssssssss file ${run_length_control}

        # Run until converged (minimum runtime 30000 steps)
        python run_length_control invoke

        unfix cr_fix # From run_length_control.py
        uncompute msd
        reset_timestep 0

        print '#================================== Thermal expansion ==================================#'

        # Thermodynamic output
        thermo_style custom lx ly lz xy yz xz temp press vol etotal step
        thermo 1000

        # Define variables for fix ave/time
        variable N_every equal 100 # sample every this many steps
        variable run_time equal 20000 # can be an input variable
        variable N_repeat equal v_run_time/v_N_every

        # compute averages of box vectors
        fix AVG all ave/time ${N_every} ${N_repeat} ${run_time} v_lx_metal v_ly_metal v_lz_metal v_xy_metal v_yz_metal v_xz_metal ave running

        # Compute unwrapped atom positions
        compute unwrapped all property/atom xu yu zu
        variable xu atom c_unwrapped[1]
        variable yu atom c_unwrapped[2]
        variable zu atom c_unwrapped[3]

        # Average unwrapped atom positions
        fix avePos all ave/atom ${N_every} ${N_repeat} ${run_time} v_xu v_yu v_zu

        # Run for thermal expansion
        run ${run_time}

        # Create per-atom variables for averaged positions
        variable ave_x atom f_avePos[1]
        variable ave_y atom f_avePos[2]
        variable ave_z atom f_avePos[3]

        # Set atom positions to time-averaged positions
        set group all x v_ave_x y v_ave_y z v_ave_z

        # Write data file of averaged positions and box dimensions
        write_data ${write_data_filename}

        unfix AVG
        unfix avePos
        unfix ensemble
        uncompute unwrapped
        reset_timestep 0

        print '#================================== Mean-squared displacement for spring constants ==================================#'

        # set NVT ensemble
        fix ensemble all nvt temp ${temp_converted} ${temp_converted} ${Tdamp_converted}

        # Compute mean squared displacement
        {compute_msd}

        # New set of variables to print to log file
        {msd_thermo_template}

        # compute average msd
        {avg_template}

        # Run for MSD spring constants
        run ${run_time}

        # compute spring constants
        variable      kB equal (8.617333262e-5)*(v__u_energy/v__u_temperature) # eV/K unless converted
        {k_lines}

        # Write final averages and spring constant
        print "# lx | ly | lz [Ang] | vol [Ang^3] | spring constant [eV/Ang^2]" file ${output_filename}
        {print_template}

        unfix ensemble
        unfix AVG
        reset_timestep 0

        print '#================================== Frenkel-Ladd NETI ==================================#'

        # Compute mean squared displacement
        compute msd all msd com yes

        fix msd_vector all vector 100 c_msd[4]
        variable msd_slope equal slope(f_msd_vector)

        # print data to logfile every 1000 timesteps
        variable      etotal_metal equal etotal/${_u_energy}
        variable      pe_metal equal pe/${_u_energy}
        variable      P_metal equal press/${_u_pressure}

        # Thermodynamic output
        thermo_style custom temp press vol etotal c_msd[4] v_msd_slope step
        thermo 1000

        # set NVE ensemble
        fix ensemble all nve

        # Frenkel-Ladd fix
        {fix_springs}
        
        # set langevin thermostat (NVE-Langevin samples spring system more accurately than NVT)
        fix           thermostat all langevin ${temp_converted} ${temp_converted} ${Tdamp_converted} ${langevin_seed} zero yes

        # run equil 1
        run ${t_equil}

        if "${msd_slope} > $(v_msd_threshold*(v__u_distance^2)/v__u_time)" then "write_dump all atom ${melted_crystal_output}" &
                                      "print 'Crystal melted or vaporized'" &
                                      "quit"

        # run switch 1
        {record_template}
                    screen no file ${switch1_output_file}
        run ${t_switch}
        unfix record

        # run equil 2
        run ${t_equil}

        variable msd_slope equal slope(f_msd_vector)

        if "${msd_slope} > $(v_msd_threshold*(v__u_distance^2)/v__u_time)" then "write_dump all atom ${melted_crystal_output}" &
                                      "print 'Crystal melted or vaporized'" &
                                      "quit"

        # run switch 2
        {record_template}
                    screen no file ${switch2_output_file}
        run ${t_switch}
        unfix record

        print '#================================== LAMMPS finished ==================================#'
        quit 0
        """

    def _add_fl_fix_for_multicomponent(self, nspecies: int):
        fix_springs_template = """
        group {group} type {group}
        fix           {fix_name} {group} ti/spring {spring_constant} {t_switch} {t_equil} function {function}
        """

        fix_entries = [
            {
                "fix_name": f"FL{i}",
                "group": f"{i+1}",
                "spring_constant": f"$(v_spring_constant_{i})",
                "t_switch": "${t_switch}",
                "t_equil": "${t_equil}",
                "function": 2,
            }
            for i in range(nspecies)
        ]

        fix_springs = "".join(
            [fix_springs_template.format(**entry) for entry in fix_entries]
        )

        self.free_energy = self.free_energy.replace("{fix_springs}", fix_springs)

        record_template = """
        fix           record all print 1 "$(pe/atoms) {data}" &
                     title "{title}" &"""

        terms = [f"f_FL{i}" for i in range(nspecies)]
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

        self.free_energy = self.free_energy.replace("{record_template}", record_template)

    def _add_msd_fix_for_multicomponent(self, nspecies: int, is_triclinic: bool):

        #===================================#

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

        self.free_energy = self.free_energy.replace("{compute_msd}", compute_msd)

        #===================================#

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

        self.free_energy = self.free_energy.replace("{msd_thermo_template}", thermo_template)

        #===================================#

        avg_template = """
        fix          AVG all ave/time ${{N_every}} ${{N_repeat}} ${{run_time}} v_lx_metal v_ly_metal v_lz_metal v_xy_metal v_yz_metal v_xz_metal {msd_data} ave running
        """

        fix_entries = [
            {
                "msd_data": msd_data,
            }
        ]

        avg_template = "".join([avg_template.format(**entry) for entry in fix_entries])

        self.free_energy = self.free_energy.replace("{avg_template}", avg_template)

        #===================================#

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

        self.free_energy = self.free_energy.replace("{k_lines}", k_lines)

        #===================================#

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

        self.free_energy = self.free_energy.replace("{print_template}", print_template)

        #===================================#
    
    def _write_lammps_file(self, nspecies: int, is_triclinic: bool):
        
        self._add_msd_fix_for_multicomponent(nspecies, is_triclinic)

        self._add_fl_fix_for_multicomponent(nspecies)

        open(self.root + "in.lammps", "w").write(self.free_energy)

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
        fix          AVG all ave/time 100 100 10000 v_lx_metal v_ly_metal v_lz_metal {msd_data} ave running
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
            }
            for i in range(nspecies)
        ]

        k_lines = "".join([k_template.format(**entry) for entry in k_entries])

        self.pre_fl = self.pre_fl.replace("{k_lines}", k_lines)

        print_template = """
        print "$(f_AVG[1]) $(f_AVG[2]) $(f_AVG[3]) $(vol) {spring_constants}" screen no append ${{output_filename}}
        """

        terms = [f"v_spring_constant_{i}" for i in range(nspecies)]
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

    def _write_pre_fl_lammps_templates(self, nspecies: int):
        self._add_msd_fix_for_multicomponent(nspecies)
        open(self.root + "preFL_template.lmp", "w").write(self.pre_fl)

    def _write_fl_lammps_templates(self, spring_constants: List[float]):
        self._add_fl_fix_for_multicomponent(spring_constants)
        open(self.root + "FL_template.lmp", "w").write(self.fl)

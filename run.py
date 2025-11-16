#!/usr/bin/python

from test_driver.test_driver import TestDriver
from ase.build import bulk

kim_model_name = "LJ_ElliottAkerson_2015_Universal__MO_959249795837_003"
test_driver = TestDriver(kim_model_name)
atoms = bulk("Au")
test_driver(atoms, temperature_K=300, lmp_cmd="mpirun -np 2 lmp", output_dir='lammps_output')

test_driver.write_property_instances_to_file()

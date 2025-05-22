#!/usr/bin/python
import subprocess
import sys

sys.path.append("../../")

import numpy as np
from kim_tools import query_crystal_structures
from test_driver.test_driver import TestDriver

# Notes:
# The reference results used Stillinger-Weber potential for silicon, with hexagonal diamond structure
# The size was ~20k atoms. On a single core, this is not feasible. This test driver will operate at ~2k atoms.

kim_model_name = "SW_StillingerWeber_1985_Si__MO_405512056662_006" # Stillinger-Weber
size = (10,10,6) # 2.4k atoms while trying to keep the system as cubic as possible
T = [1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000] # Temperatures at which to calculate the free-energy

save_file = "./output/results/Free-Energy.dat"

np.savetxt(save_file, [], delimiter=" ", header='T [K] | F [eV/atom]')

#======================================================================================#

subprocess.run(f"kimitems install {kim_model_name}", shell=True, check=True)

test_driver = TestDriver(kim_model_name)

list_of_queried_structures = query_crystal_structures(
    kim_model_name=kim_model_name,
    stoichiometric_species=['Si'],
    prototype_label="A_hP4_194_f", # hexagonal diamond
)

#======================================================================================#

for t in range(len(T)):

    print(f"\nRunning T = {T[t]}K\n")
    computed_property_instances = test_driver(
        list_of_queried_structures[0],
        size=size,
        temperature_K=T[t]
    )

    with open(save_file, 'a') as f:
        np.savetxt(f, np.column_stack([T[t], computed_property_instances[-1]]), delimiter=" ", fmt='%.3f %.7f')

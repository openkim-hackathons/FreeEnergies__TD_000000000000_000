#!/usr/bin/python
import subprocess
import sys

sys.path.append("../../")

import numpy as np
from kim_tools import query_crystal_structures
from test_driver.test_driver import TestDriver

kim_model_name = "SW_StillingerWeber_1985_Si__MO_405512056662_006" # Stillinger-Weber
size = (10,10,6)
T = []

save_file = "./output/results/Free-Energy.dat"

Ref_T = [1100, 1200, 1300, 1400, 1500.7142857142858, 1600.7142857142858, 1700.7142857142858, 1800.7142857142858, 1900.7142857142858, 2002.142857142857] # Temperatures from doi:10.1038/s41467-020-16892-4, supporting information (Fig 11b)
Ref_F = [-4.57429718875502, -4.624096385542169, -4.675502008032129, -4.7301204819277105, -4.785542168674699, -4.844176706827309, -4.903614457831325, -4.965461847389558, -5.028915662650602, -5.093172690763052] # Free-energy from doi:10.1038/s41467-020-16892-4, supporting information (Fig 11b)

F_results = []

# Notes:
# The reference results used Stillinger-Weber potential for silicon, with hexagonal diamond structure
# The size was ~20k atoms. On a single core, this is not feasible. This test driver will operate at ~2k atoms.

np.savetxt(save_file, )

subprocess.run(f"kimitems install {kim_model_name}", shell=True, check=True)

test_driver = TestDriver(kim_model_name)

list_of_queried_structures = query_crystal_structures(
    kim_model_name=kim_model_name,
    stoichiometric_species=['Si'],
    prototype_label="A_hP4_194_f", # hexagonal diamond
)

#======================================================================================#

for t in range(len(T)):

    print(f"\nRunning T = {T[t]}\n")
    computed_property_instances = test_driver(
        list_of_queried_structures[0],
        size=size,
        temperature_K=T[t]
    )

    with open(save_file, 'a') as f:
        np.savetxt(f, np.column_stack([T[t], computed_property_instances[-1]]), delimiter=" ", fmt='%.3f %.7f')

#!/usr/bin/python
import os
import subprocess
import sys

os.makedirs("output/results", exist_ok=True)
sys.path.append("../../")


from cte import TEMPERATURES
from kim_tools import query_crystal_structures
from test_driver.test_driver import TestDriver

kim_model_name = "SW_StillingerWeber_1985_Si__MO_405512056662_006" # Stillinger-Weber
subprocess.run(f"kimitems install {kim_model_name}", shell=True, check=True)


test_driver = TestDriver(kim_model_name)

list_of_queried_structures = query_crystal_structures(
    kim_model_name=kim_model_name,
    stoichiometric_species=['Si'],
    prototype_label="A_hP4_194_f", # hexagonal diamond
)

for t in range(len(TEMPERATURES)):

    print(f"\nRunning T = {TEMPERATURES[t]}K\n")
    computed_property_instances = test_driver(
        list_of_queried_structures[0],
        size=(10,10,6), # The reference results used Stillinger-Weber potential for silicon, with hexagonal diamond structure.The size was ~20k atoms. On a single core, this is not feasible. This test driver will operate at ~2k atoms, while trying to keep the system as cubic as possible.
        temperature_K=TEMPERATURES[t]
    )

    test_driver.write_property_instances_to_file(filename = f"output/results/results_{TEMPERATURES[t]}K.edn")









#!/usr/bin/python
import sys
from test_driver.test_driver import TestDriver
import subprocess
from kim_tools import query_crystal_structures

kim_model_name = "SW_StillingerWeber_1985_Si__MO_405512056662_006" # Stillinger-Weber
#kim_model_name = "EDIP_JustoBazantKaxiras_1998_Si__MO_958932894036_002" # EDIP
#kim_model_name = "SNAP_ZuoChenLi_2019_Si__MO_869330304805_000" # SNAP
#kim_model_name = "SNAP_ZuoChenLi_2019quadratic_Si__MO_721469752060_000" # qSNAP

subprocess.run(f"kimitems install {kim_model_name}", shell=True, check=True)

test_driver = TestDriver(kim_model_name)

species = "Si"

#prototype_label="A_cF8_227_a" # cubic diamond
prototype_label="A_hP4_194_f" # hexagonal diamond

list_of_queried_structures = query_crystal_structures(
    kim_model_name=kim_model_name,
    stoichiometric_species=['Si'],
    #prototype_label="A_cF8_227_a", # cubic diamond
    prototype_label="A_hP4_194_f", # hexagonal diamond
)

#======================================================================================#

print("\nRUNNING TEST DRIVER ON QUERIED STRUCTURE\n")
computed_property_instances = test_driver(
    list_of_queried_structures[0],
    size=(3,3,3),
    pressure=0.0,
    temperature_K=1684,
)

#======================================================================================#

print("\n--------------------------------------")
print(
    "Free energy G_FL (eV/atom): %f %s"
    % (
        computed_property_instances[1]["gibbs_free_energy_per_atom"]["source-value"],
        computed_property_instances[1]["gibbs_free_energy_per_atom"]["source-unit"],
    )
)
print("--------------------------------------\n")
print("Done")

#test_driver.write_property_instances_to_file()
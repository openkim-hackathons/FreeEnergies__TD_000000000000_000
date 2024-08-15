#!/usr/bin/python
from test_driver.test_driver import TestDriver
import subprocess
from kim_tools import query_crystal_genome_structures

kim_model_name ="EAM_Dynamo_ErcolessiAdams_1994_Al__MO_123629422045_005"
#kim_model_name = "SW_StillingerWeber_1985_Si__MO_405512056662_006" # Stillinger-Weber
#kim_model_name = "EDIP_JustoBazantKaxiras_1998_Si__MO_958932894036_002" # EDIP
#kim_model_name = "SNAP_ZuoChenLi_2019_Si__MO_869330304805_000" # SNAP
#kim_model_name = "SNAP_ZuoChenLi_2019quadratic_Si__MO_721469752060_000" # qSNAP
subprocess.run(f"kimitems install {kim_model_name}", shell=True, check=True)
test_driver = TestDriver(kim_model_name)

list_of_queried_structures = query_crystal_genome_structures(
    kim_model_name=kim_model_name,
    stoichiometric_species=['Al'],
    prototype_label="A_cF4_225_a",
    #prototype_label="A_cF8_227_a", # cubic diamond
    #prototype_label="A_hP4_194_f", # hexagonal diamond
)

for queried_structure in list_of_queried_structures:
    print("\nRUNNING TEST DRIVER ON QUERIED STRUCTURE\n")
    test_driver(
        **queried_structure,
        size=(5,5,5),
        temperature=100,
        pressure=0.0,
    )

print("\n--------------------------------------")
print(
    "Free energy G_FL (eV/cell): %f %s"
    % (
        test_driver.property_instances[1]["free_energy"]["source-value"],
        test_driver.property_instances[1]["free_energy"]["source-unit"],
    )
)
print("--------------------------------------\n")

test_driver.write_property_instances_to_file()

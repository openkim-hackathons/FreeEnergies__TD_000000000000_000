#!/usr/bin/python
from test_driver.test_driver import TestDriver
import subprocess
from kim_tools import query_crystal_genome_structures

kim_model_name = "Sim_LAMMPS_Buckingham_FreitasSantosColaco_2015_SiCaOAl__SM_154093256665_000"
subprocess.run(f"kimitems install {kim_model_name}", shell=True, check=True)
test_driver = TestDriver(kim_model_name)

list_of_queried_structures = query_crystal_genome_structures(
    kim_model_name=kim_model_name,
    stoichiometric_species=['Ca', 'O', 'Si'],
    prototype_label="AB3C_aP30_2_3i_9i_3i",
)

for queried_structure in list_of_queried_structures:
    print("\nRUNNING TEST DRIVER ON QUERIED STRUCTURE\n")
    test_driver(
        **queried_structure,
        size=(3,3,3),
        temperature=100.0,
        pressure=0.0,
    )

print("\n--------------------------------------")
print(
    "Free energy G_FL (eV/cell): %f %s"
    % (
        test_driver.property_instances[0]["free_energy"]["source-value"],
        test_driver.property_instances[0]["free_energy"]["source-unit"],
    )
)
print("--------------------------------------\n")

test_driver.write_property_instances_to_file()
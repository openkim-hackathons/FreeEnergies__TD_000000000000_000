#!/usr/bin/python

from test_driver.test_driver import TestDriver
from kim_tools import crystal_input_from_test_generator_line
from json import dump
import subprocess

calcs = [
    ["Sim_LAMMPS_Buckingham_MatsuiAkaogi_1991_TiO__SM_690504433912_000"],
    ["MEAM_LAMMPS_JeongParkDo_2018_PdAl__MO_616482358807_002"],
    [
        "Tersoff_LAMMPS_MunetohMotookaMoriguchi_2007_SiO__MO_501246546792_000",
        "Sim_LAMMPS_Vashishta_BroughtonMeliVashishta_1997_SiO__SM_422553794879_000",
    ],
    [
        "EAM_Dynamo_AcklandMendelevSrolovitz_2004_FeP__MO_884343146310_005",
        "MEAM_LAMMPS_KoJimLee_2012_FeP__MO_179420363944_002",
    ],
    ["Sim_LAMMPS_BOP_MurdickZhouWadley_2006_GaAs__SM_104202807866_001"],
    ["Sim_LAMMPS_ReaxFF_BrugnoliMiyataniAkaji_SiCeNaClHO_2023__SM_282799919035_000"],
    ["EAM_Dynamo_ErcolessiAdams_1994_Al__MO_123629422045_005"],
]

results = []
with open("test_generator.json") as f:
    for calc_set, line in zip(calcs, f):
        for calc in calc_set:            
            subprocess.run(f"kim-api-collections-management install system {calc}", shell=True, check=True)
            free_e = TestDriver(calc)
            inputs = crystal_input_from_test_generator_line(line, calc)
            for input in inputs:
                try:                
                    results += free_e(**input)
                    with open("runs.json", "w") as f:
                        dump(results, f, indent=2)
                except Exception as e:
                    print(f"Got exception {repr(e)}")
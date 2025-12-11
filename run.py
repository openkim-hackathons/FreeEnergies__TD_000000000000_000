#!/usr/bin/python

from test_driver.test_driver import TestDriver
from ase.build import bulk

# ===============================================================================================#

from kim_tools import crystal_input_from_test_generator_line
from json import dump
import subprocess

calcs = [
    ["EAM_Dynamo_ErcolessiAdams_1994_Al__MO_123629422045_006"],
]

results = []
with open("test_generator.json") as f:
    for calc_set, line in zip(calcs, f):
        for calc in calc_set:
            subprocess.run(f"kimitems install -D {calc}", shell=True, check=True)
            free_e = TestDriver(calc)
            inputs = crystal_input_from_test_generator_line(line, calc)
            for input in inputs:
                try:
                    results += free_e(**input)
                    with open("runs.json", "w") as f:
                        dump(results, f, indent=2)
                except Exception as e:
                    print(f"Got exception {repr(e)}")
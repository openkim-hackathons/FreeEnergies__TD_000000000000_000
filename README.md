# Free Energy Crystal Genome Test Driver

Gibbs free energy per conventional unit cell, with rhombohedral crystals using the hexagonal cell, of a crystal at constant temperature and pressure using the Frenkel-Ladd Hamiltonian integration algorithm (Freitas, R., Asta, M. & de Koning, M. Nonequilibrium free-energy calculation of solids using LAMMPS. Computational Materials Science 112, 333–341 (2016)). The free energy is computed through one equilibrium NPT simulation ('preFL') and one NONequilibrium NVT simulation ('FL').

## Example of usage 
```python 
model_name = "LJ_Shifted_Bernardes_1958MedCutoff_Ar__MO_126566794224_004"
subprocess.run(f"kimitems install {model_name}", shell=True, check=True)
subprocess.run("mkdir -p output", shell=True, check=True)
test_driver = TestDriver(model_name)
test_driver(
    bulk("Ar", "fcc", a=5.248),
    size=(3, 3, 3),
    temperature=20.0,
    pressure=0.0,
)
```

## Notes
Frenkel-Ladd method is more accurate for larger systems (~10k atoms) and longer nonequilibrium simulation time ($\geq 10\text{k}$ integration steps)

## Citation
Citation for the algorithm implented by this test driver:

```
@article{freitas_nonequilibrium_2016,
title = {Nonequilibrium free-energy calculation of solids using LAMMPS},
journal = {Computational Materials Science},
volume = {112},
pages = {333-341},
year = {2016},
issn = {0927-0256},
doi = {https://doi.org/10.1016/j.commatsci.2015.10.050},
author = {Rodrigo Freitas and Mark Asta and Maurice {de Koning}},
}
```

## Contact
If any troubles, feel free to contact us at ``ksheriff at mit dot edu``, and ``jogbebor at mit dot edu``


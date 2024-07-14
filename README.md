# FreeEnergies__TD_000000000000_000

Gibbs free energy per unit cell of a crystal at constant temperature and pressure using Frenkel-Ladd Hamiltonian integration algorithm. Computed through one equilibrium NPT simulation ('preFL') and one NONequilibrium NVT simulation ('FL').

## Example of usage 
```python 
model_name = "LJ_Shifted_Bernardes_1958MedCutoff_Ar__MO_126566794224_004"
subprocess.run(f"kimitems install {model_name}", shell=True, check=True)
subprocess.run("mkdir -p output", shell=True, check=True)
test_driver = FrenkelLaddFreeEnergies(model_name)
test_driver(
    bulk("Ar", "fcc", a=5.248),
    size=(3, 3, 3),
    temperature=20.0,
    pressure=0.0,
)
```

## Notes
Frenkel-Ladd method is more accurate for larger systems (~10k atoms) and longer nonequilibrium simulation time ($\geq$ 10k integration steps)

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


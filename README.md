# FreeEnergies__TD_000000000000_000

Gibbs free energy of a crystal at constant temperature and pressure using Frenkel-Ladd Hamiltonian integration algorithm. Computed through one equilibrium NPT simulation ('preFL') and one NONequilibrium NVT simulation ('FL').

## Example of usage 
```python 
model_name = "LJ_Shifted_Bernardes_1958MedCutoff_Ar__MO_126566794224_004"
subprocess.run(f"kimitems install {model_name}", shell=True, check=True)
test_driver = FrenkelLaddFreeEnergies(model_name)
test_driver(
    bulk("Ar", "fcc", a=5.248),
    size=(3, 3, 3),
    temperatures=100.0,
    pressures=1.0,
)
```

## Caution
- "fix ti/spring" is part of the EXTRA-FIX package, which isn't enabled in the lammps build that the test driver is using. 
- preFL computes equilibrium lattice parameters, yet right now it would only work for cubic systems. Do we have a test driver specifically for this, that works for arbitrary crystal structures?

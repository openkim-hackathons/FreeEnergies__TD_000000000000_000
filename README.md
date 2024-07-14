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

# Notes/Caution
Frenkel-Ladd method is most accurate for large systems (~10k atoms) and long nonequilibrium simulation time (>=)

# Contact
If any troubles, feel free to contact us at ``ksheriff at mit dot edu``, and ``jogbebor at mit dot edu``
# FreeEnergies__TD_000000000000_000

Gibbs free energy of a crystal at constant temperature and pressure using Frenkel-Ladd Hamiltonian integration algorithm. Computed through one equilibrium NPT simulation ('preFL') and one NONequilibrium NVT simulation ('FL').
  
## Some notes
- Free energy calculations broken down into three steps: preFL, FL, and RS.
- preFL computes the lattice parameter and spring constants as functions of pressure at the starting temperature (T = 100K).
- FL computes the free energy as a function of pressure at the starting temperature.
- RS computes the free energy at each pressure over a mesh of temperatures.

Each step has a lammps template written in the style of the example temp driver (can be merged later).

The runner has been expanded in the same style as the example driver, with additional functions for each step. All free energies are computed in runner.py via integration of data obtained from lammps.

Some parameters (switching times, equilibration times, starting temperature) are currently hard coded, but could be rewritten as input variables.

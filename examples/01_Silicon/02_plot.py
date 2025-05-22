import lovelyplots
import matplotlib.pyplot as plt
from cte import FREE_ENERGIES, TEMPERATURES
from kim_edn import load

plt.style.use('paper')

# Read in the computed free-energies from the edn files
computed_free_energies = []

for temperature in TEMPERATURES:
    data = load(f"output/results/results_{temperature}K.edn")
    free_energy = data[-1]["gibbs_free_energy_per_atom"]["source-value"]
    computed_free_energies.append(free_energy)

# Plot the computed free-energies against the reference free-energies
fig, ax = plt.subplots()

ax.plot(TEMPERATURES, FREE_ENERGIES, label='Freitas et al.')
ax.plot(TEMPERATURES, computed_free_energies, label='This test driver')

ax.set_xlabel('Temperature [K]')
ax.set_ylabel('Free energy [eV/atom]')

ax.legend(loc='best', frameon=True)
fig.savefig('free_energies.pdf')






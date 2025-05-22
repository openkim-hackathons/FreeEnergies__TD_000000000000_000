import lovelyplots
import matplotlib.pyplot as plt
from cte import FREE_ENERGIES, TEMPERATURES

plt.style.use('paper')






# Read in the computed free-energies from the edn files
computed_free_energies = []







# Plot the computed free-energies against the reference free-energies
fig, ax = plt.subplots()

ax.plot(TEMPERATURES, FREE_ENERGIES, label='Freitas et al.')
ax.plot(TEMPERATURES, computed_free_energies, label='This test driver')

ax.set_xlabel('Temperature [K]')
ax.set_ylabel('Free energy [eV/atom]')

ax.legend(loc='best', frameon=True)
fig.savefig('free_energies.pdf')
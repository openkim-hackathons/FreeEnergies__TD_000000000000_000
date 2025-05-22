import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

timestep = 0.001 # ps

data = np.loadtxt("output/equilibration.dat", skiprows=1)

step = data[:,0] * timestep
vol = data[:,1]
temp = data[:,2]
lx = data[:,3]
ly = data[:,4]
lz = data[:,5]
xy = data[:,6]
xz = data[:,7]
yz = data[:,8]

#=============================================================#

fig = plt.figure()

plt.plot(step, temp, color='red', marker='none')

plt.xlabel('Time [ps]')
plt.ylabel('T [K]')

plt.savefig('./output/figures/T.pdf')
plt.close(fig)

#=============================================================#

fig = plt.figure()

plt.plot(step, vol, color='red', marker='none')

plt.xlabel('Time [ps]')
plt.ylabel(r'V [Ang$^3$]')

plt.savefig('./output/figures/V.pdf')
plt.close(fig)

#=============================================================#

fig = plt.figure()

plt.plot(step, lx, color='red', marker='none', label='lx')
plt.plot(step, ly, color='blue', marker='none', label='ly')
plt.plot(step, lz, color='green', marker='none', label='lz')

plt.xlabel('Time [ps]')
plt.ylabel(r'[Ang]')

plt.legend(loc='best')

plt.savefig('./output/figures/lxlylz.pdf')
plt.close(fig)

#=============================================================#

fig = plt.figure()

plt.plot(step, xy, color='red', marker='none', label='xy')
plt.plot(step, xz, color='blue', marker='none', label='xz')
plt.plot(step, yz, color='green', marker='none', label='yz')

plt.xlabel('Time [ps]')
plt.ylabel(r'[Ang]')

plt.legend(loc='best')

plt.savefig('./output/figures/tilt_factors.pdf')
plt.close(fig)

#=============================================================#

print("Done")
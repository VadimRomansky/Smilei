from matplotlib import animation
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np
import h5py
def plot_smilei_concentration1d(ntot, file_name, sampling, dx, prefix=''):
    print("plot concentration 1d")
    plt.rcParams.update({'font.size': 40})
    plt.rcParams['text.usetex'] = True
    plt.rcParams['axes.linewidth'] = 4
    f1 = plt.figure(figsize=[10, 10])
    ax = f1.add_subplot(111)
    ax.tick_params(axis='x', size=10, width=4)
    ax.tick_params(axis='y', size=10, width=4)
    ax.minorticks_on()

    file = h5py.File(file_name, 'r')
    print("Keys: %s" % file.keys())
    l = list(file.keys())

    V = np.array(file.get(l[ntot-1])).T
    Nx = V.shape[0]

    x = np.zeros([Nx])
    for i in range(Nx):
        x[i] = i*sampling*dx

    ax.plot(x, V, linewidth=4)
    ax.set_xlabel(r'$X \omega_e/c$', fontsize=40,fontweight='bold')
    ax.set_ylabel(r'$n$', fontsize=40,fontweight='bold')
    ax.minorticks_on()
    plt.savefig('smilei_concentration1d'+prefix +'.png', bbox_inches='tight')



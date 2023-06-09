from matplotlib import animation
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np
import h5py
def plot_smilei_spectrum(ntot, file_name, prefix, minE, maxE, xmin, xmax):
    print("plot spectrum")
    plt.rcParams.update({'font.size': 40})
    plt.rcParams['text.usetex'] = True
    plt.rcParams['axes.linewidth'] = 4
    f1 = plt.figure(figsize=[10,10])
    ax = f1.add_subplot(111)
    ax.tick_params(axis='x', size=10, width=4)
    ax.tick_params(axis='y', size=10, width=4)
    ax.minorticks_on()

    file = h5py.File(file_name, 'r')
    print("Keys: %s" % file.keys())
    l = list(file.keys())

    V = np.array(file.get(l[ntot]))


    Nx = V.shape[0]
    Np = V.shape[1]

    factor = (maxE / minE) ** (1.0 / (Np - 1));

    energy = np.zeros([Np])
    de = np.zeros([Np])
    energy[0] = minE
    for i in range(1,Np):
        energy[i] = energy[i-1]*factor
    de[0] = energy[1]-energy[0]
    for i in range(1,Np):
        de[i] = energy[i]-energy[i-1]

    f = np.zeros([Np])
    for i in range(1,Np):
        for j in range(xmin,xmax):
            f[i] = f[i] + V[j][i]/de[i]

    minF = np.amin(f)
    maxF = np.amax(f)
    minF = maxF / 1E14

    ax.plot(energy, f, linewidth=4)  # plotting fluid data.
    ax.set_xlabel(r'$E_{kin}/m_e c^2$', fontsize=40,fontweight='bold')
    ax.set_ylabel(r'$F(E_{kin})$', fontsize=40,fontweight='bold')
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.savefig('smilei_distribution' + prefix + '.png', bbox_inches='tight')


from matplotlib import animation
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np
import h5py
def plot_smilei_momentum_spectrum(ntot, file_name, prefix, minE, maxE, xmin, xmax, mass):
    print("plot momentum spectrum")
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

    V = np.array(file.get(l[ntot-1]))


    Nx = V.shape[0]
    Np = V.shape[1]

    factor = (maxE / minE) ** (1.0 / (Np - 1));

    energy = np.zeros([Np])
    momentum = np.zeros([Np])
    de = np.zeros([Np])
    energy[0] = minE
    momentum[0] = np.sqrt((energy[0]+mass)*(energy[0]+mass) - mass*mass)
    for i in range(1,Np):
        energy[i] = energy[i-1]*factor
        momentum[i] = np.sqrt((energy[i] + mass) * (energy[i] + mass) - mass * mass)
    de[0] = energy[1]-energy[0]
    for i in range(1,Np):
        de[i] = energy[i]-energy[i-1]

    f = np.zeros([Np])
    for i in range(1,Np):
        for j in range(xmin,xmax):
            f[i] = f[i] + (V[j][i]/de[i])*momentum[i]/(energy[i]+mass)

    minF = np.amin(f)
    maxF = 2*np.amax(f)
    minF = maxF / 1E8

    ax.plot(momentum, f, linewidth=4)  # plotting fluid data.
    ax.set_ylim([minF, maxF])
    ax.set_xlabel(r'$p/m_e c$', fontsize=40,fontweight='bold')
    ax.set_ylabel(r'$F(p)$', fontsize=40,fontweight='bold')
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.savefig('smilei_momentum_distribution' + prefix + '.png', bbox_inches='tight')
    plt.close()


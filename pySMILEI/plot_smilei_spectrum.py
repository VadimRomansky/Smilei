from matplotlib import animation
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np
import h5py
def plot_smilei_spectrum(ntot, file_name, prefix, xmin, xmax):
    print("plot spectrum")
    plt.rcParams.update({'font.size': 40})
    plt.rcParams['text.usetex'] = True
    plt.rcParams['axes.linewidth'] = 4
    f1 = plt.figure(figsize=[10,8])
    ax = f1.add_subplot(111)


    file = h5py.File(file_name, 'r')
    print("Keys: %s" % file.keys())
    l = list(file.keys())

    V = np.array(file.get(l[ntot]))


    Nx = V.shape[0]
    Np = V.shape[1]

    minEe = 0.001;
    maxEe = 5000;
    minEp = 0.1;
    maxEp = 5000;
    minE = minEe;
    maxE = maxEe;

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

    ax.plot(energy, f)  # plotting fluid data.
    ax.set_xlabel(r'$E_{kin}/me c^2$', fontsize=18)
    ax.set_ylabel(r'$F(E_{kin})$', fontsize=18)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.minorticks_on()
    plt.savefig('smilei_distribution' + prefix + '.png', bbox_inches='tight')


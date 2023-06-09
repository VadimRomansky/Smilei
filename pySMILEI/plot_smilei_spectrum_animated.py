from matplotlib import animation
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np
import h5py
def plot_smilei_spectrum_animated(ntot, file_name, prefix, mass, minE, maxE, xmin, xmax):
    print("plot spectrum animated")
    plt.rcParams.update({'font.size': 40})
    plt.rcParams['text.usetex'] = True
    plt.rcParams['axes.linewidth'] = 4
    f1 = plt.figure(figsize=[10, 8])
    ax = f1.add_subplot(111)
    ax.tick_params(axis='x', size=10, width=4)
    ax.tick_params(axis='y', size=10, width=4)
    ax.minorticks_on()



    file = h5py.File(file_name, 'r')
    print("Keys: %s" % file.keys())
    l = list(file.keys())

    V = np.array(file.get(l[ntot])).T
    Nx = V.shape[1]
    Np = V.shape[0]

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
    norm = 0
    for i in range(Np):
        for j in range(xmin,xmax):
            #f[i] = f[i] + V[i][j]*(energy[i] + mass) * (energy[i] + mass)/de[i]
            f[i] = f[i] + V[i][j]/de[i]
            norm = norm+V[i][j]

    for i in range(Np):
        f[i] = f[i]/norm

    minF = np.amin(f)
    maxF = np.amax(f)

    for i in range(ntot):
        f = np.zeros([Np])
        norm = 0
        V = np.array(file.get(l[i])).T
        for k in range(Np):
            for j in range(xmin, xmax):
                #f[k] = f[k] + V[k][j] *(energy[k] + mass) * (energy[k] + mass) / de[k]
                f[k] = f[k] + V[k][j]/ de[k]
                norm = norm + V[k][j]

        for k in range(Np):
            f[k] = f[k]/norm

        if (np.amin(f) < minF):
            minF = np.amin(f)
        if (np.amax(f) > maxF):
            maxF = np.amax(f)

    maxF = 2*maxF
    minF = maxF / 1E6

    ax.plot(energy, f)  # plotting fluid data.
    ax.set_xlabel(r'$E_{kin}/m_e c^2$', fontsize=40,fontweight='bold')
    ax.set_ylabel(r'$F(E_{kin}$', fontsize=40,fontweight='bold')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.minorticks_on()

    def update(frame_number):
        print(frame_number)
        ax.clear()
        ax.set_ylim([minF, maxF])
        ax.set_xlabel(r'$E_{kin}/m_e c^2$', fontsize=40, fontweight='bold')
        ax.set_ylabel(r'$F(E_{kin}$', fontsize=40, fontweight='bold')
        ax.set_xscale('log')
        ax.set_yscale('log')
        V = np.array(file.get(l[frame_number])).T
        f = np.zeros([Np])
        norm = 0
        for i in range(Np):
            for j in range(xmin, xmax):
                #f[i] = f[i] + V[i][j] *(energy[i] + mass) * (energy[i] + mass)/ de[i]
                f[i] = f[i] + V[i][j]/ de[i]
                norm = norm +V[i][j]

        for i in range(Np):
            f[i] = f[i]/norm

        im2 = ax.plot(energy, f, linewidth=4)
        return im2

    anim = FuncAnimation(f1, update, interval=10, frames=ntot)

    # plt.show()

    f = r"smilei_distribution" + prefix + ".gif"
    writergif = animation.PillowWriter(fps=4)
    anim.save(f, writer=writergif)

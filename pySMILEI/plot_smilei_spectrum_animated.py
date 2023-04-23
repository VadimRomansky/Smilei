from matplotlib import animation
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np
import h5py
def plot_smilei_spectrum_animated(ntot, file_name, prefix, xmin, xmax):
    f1 = plt.figure(figsize=[10,8])
    ax = f1.add_subplot(111)


    file = h5py.File(file_name, 'r')
    print("Keys: %s" % file.keys())
    l = list(file.keys())

    V = np.array(file.get(l[ntot])).T
    Nx = V.shape[1]
    Np = V.shape[0]

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
    for i in range(Np):
        for j in range(xmin,xmax):
            f[i] = f[i] + V[i][j]/de[i]
    minF = np.amin(f)
    maxF = np.amax(f)

    for i in range(ntot):
        V = np.array(file.get(l[i])).T
        for k in range(Np):
            for j in range(xmin, xmax):
                f[k] = f[k] + V[k][j] / de[k]
        if (np.amin(f) < minF):
            minF = np.amin(f)
        if (np.amax(f) > maxF):
            maxF = np.amax(f)

    minF = maxF / 1E14

    ax.plot(energy, f)  # plotting fluid data.
    ax.set_xlabel(r'Ekin/me c^2', fontsize=18)
    ax.set_ylabel(r'F', fontsize=18)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.minorticks_on()

    def update(frame_number):
        print(frame_number)
        ax.clear()
        ax.set_ylim([minF, maxF])
        ax.set_xlabel(r'Ekin/me c^2', fontsize=18)
        ax.set_ylabel(r'F', fontsize=18)
        ax.set_xscale('log')
        ax.set_yscale('log')
        V = np.array(file.get(l[frame_number])).T
        f = np.zeros([Np])
        for i in range(Np):
            for j in range(xmin, xmax):
                f[i] = f[i] + V[i][j] / de[i]
        im2 = ax.plot(energy, f)
        return im2

    anim = FuncAnimation(f1, update, interval=10, frames=ntot)

    # plt.show()

    f = r"smilei_distribution" + prefix + ".gif"
    writergif = animation.PillowWriter(fps=4)
    anim.save(f, writer=writergif)

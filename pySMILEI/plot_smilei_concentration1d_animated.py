from matplotlib import animation
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np
import h5py
def plot_smilei_concentration1d_animated(ntot, file_name, sampling, dx, prefix=''):
    print("plot concentration 1d animated")
    plt.rcParams.update({'font.size': 40})
    plt.rcParams['text.usetex'] = True
    plt.rcParams['axes.linewidth'] = 4
    f1 = plt.figure(figsize=[10, 10])
    ax = f1.add_subplot(111)
    ax.tick_params(axis='x', size=10, width=4)
    ax.tick_params(axis='y', size=10, width=4)
    ax.minorticks_on()

    #todo
    V0 = 8.0

    file = h5py.File(file_name, 'r')
    print("Keys: %s" % file.keys())
    l = list(file.keys())

    V = np.array(file.get(l[ntot-1])).T
    Nx = V.shape[0]
    minV = np.min(V)/V0;
    maxV = np.max(V)/V0;
    minV=0
    for i in range(ntot):
        V = np.array(file.get(l[i])).T/V0
        if (np.amin(V) < minV):
            minV = np.amin(V)
        if (np.amax(V) > maxV):
            maxV = np.amax(V)

    maxV = 1.5*maxV
    maxV = 10
    x = np.zeros([Nx])
    for i in range(Nx):
        x[i] = i*sampling*dx

    ax.plot(x, V, linewidth=4)
    ax.set_xlabel(r'$X \omega_e/c$', fontsize=40,fontweight='bold')
    ax.set_ylabel(r'$n/n_0$', fontsize=40,fontweight='bold')
    ax.set_ylim([minV, maxV])
    ax.minorticks_on()

    def update(frame_number):
        print(frame_number)
        ax.clear()
        ax.set_ylim([minV, maxV])
        V = np.array(file.get(l[frame_number])).T/V0
        im2 = ax.plot(x, V, linewidth=4)
        ax.set_xlabel(r'$X \omega_e/c$', fontsize=40, fontweight='bold')
        ax.set_ylabel(r'$n/n_0$', fontsize=40, fontweight='bold')
        return im2

    anim = FuncAnimation(f1, update, interval=10, frames=ntot)

    # plt.show()

    f = r"smilei_concentration1d" + prefix + ".gif"
    writergif = animation.PillowWriter(fps=4)
    anim.save(f, writer=writergif)
    plt.close()


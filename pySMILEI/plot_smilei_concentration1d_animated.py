from matplotlib import animation
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np
import h5py
def plot_smilei_concentration1d_animated(ntot, file_name):
    f1 = plt.figure(figsize=[10,8])
    ax = f1.add_subplot(111)

    file = h5py.File(file_name, 'r')
    print("Keys: %s" % file.keys())
    l = list(file.keys())

    V = np.array(file.get(l[ntot])).T
    Nx = V.shape[0]
    minV = np.min(V);
    maxV = np.max(V);
    for i in range(ntot):
        V = np.array(file.get(l[i])).T
        if (np.amin(V) < minV):
            minV = np.amin(V)
        if (np.amax(V) > maxV):
            maxV = np.amax(V)

    x = np.zeros([Nx])
    for i in range(Nx):
        x[i] = i

    ax.plot(x, V)
    ax.set_xlabel(r'X', fontsize=18)
    ax.set_ylabel(r'n', fontsize=18)
    ax.minorticks_on()

    def update(frame_number):
        print(frame_number)
        ax.clear()
        ax.set_ylim([minV, maxV])
        V = np.array(file.get(l[frame_number])).T
        im2 = ax.plot(x, V)
        return im2

    anim = FuncAnimation(f1, update, interval=10, frames=ntot)

    # plt.show()

    f = r"smilei_concentration1d.gif"
    writergif = animation.PillowWriter(fps=4)
    anim.save(f, writer=writergif)


from matplotlib import animation
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np
import h5py
def plot_smilei_concentration1d(ntot, file_name, prefix=''):
    f1 = plt.figure(figsize=[10,8])
    ax = f1.add_subplot(111)

    file = h5py.File(file_name, 'r')
    print("Keys: %s" % file.keys())
    l = list(file.keys())

    V = np.array(file.get(l[ntot])).T
    Nx = V.shape[0]

    x = np.zeros([Nx])
    for i in range(Nx):
        x[i] = i

    ax.plot(x, V)
    ax.set_xlabel(r'X', fontsize=18)
    ax.set_ylabel(r'n', fontsize=18)
    ax.minorticks_on()
    plt.savefig('smilei_concentration1d'+prefix +'.png')


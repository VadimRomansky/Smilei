from matplotlib import animation
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np
import h5py
def plot_smilei_fields1d_animated(ntot, file_name, parameter_name, Nmin, Nmax):
    print("plot fields 1d")
    f1 = plt.figure(figsize=[10,8])
    ax = f1.add_subplot(111)

    file = h5py.File(file_name, 'r')
    print("Keys: %s" % file.keys())
    l = list(file.keys())
    a_group_key = list(file.keys())[0]
    data = list(file[a_group_key])
    str = '/data/'+data[ntot] +'/' + parameter_name
    V = np.array(file[str][:])
    
    Nx = V.shape[0]
    Ny = V.shape[1]

    
    V1 = V[:,int(Ny/2)]
    
    minV = np.min(V1);
    maxV = np.max(V1);
    for i in range(ntot):
        str = '/data/' + data[i] +'/' + parameter_name
        V = np.array(file[str][:])
    
        V1 = V[:,int(Ny/2)]
        if (np.amin(V1) < minV):
            minV = np.amin(V1)
        if (np.amax(V1) > maxV):
            maxV = np.amax(V1)

    x = np.zeros([Nx])
    for i in range(Nx):
        x[i] = i

    ax.plot(x, V1)
    ax.set_xlabel(r'X', fontsize=18)
    ax.set_ylabel(r'n', fontsize=18)
    ax.minorticks_on()

    def update(frame_number):
        print(frame_number)
        ax.clear()
        ax.set_ylim([minV, maxV])
        str = '/data/' + data[frame_number] +'/' + parameter_name
        V = np.array(file[str][:])
        V1 = V[:,int(Ny/2)]
        im2 = ax.plot(x[Nmin:Nmax], V1[Nmin:Nmax])
        return im2

    anim = FuncAnimation(f1, update, interval=10, frames=ntot)

    # plt.show()

    f = r"smilei_1d_" + parameter_name + ".gif"
    writergif = animation.PillowWriter(fps=4)
    anim.save(f, writer=writergif)


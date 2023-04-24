from matplotlib import animation
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np
import h5py
def plot_smilei_fields1d(ntot, file_name, parameter_name, Nmin, Nmax):
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

    x = np.zeros([Nx])
    for i in range(Nx):
        x[i] = i

    ax.plot(x[Nmin:Nmax], V1[Nmin:Nmax])
    ax.set_xlabel(r'X', fontsize=18)
    ax.set_ylabel(r'B', fontsize=18)
    ax.minorticks_on()
    plt.savefig('smilei_1d_'+ parameter_name +'.png')


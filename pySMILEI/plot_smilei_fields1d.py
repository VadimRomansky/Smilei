from matplotlib import animation
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np
import h5py
def plot_smilei_fields1d(ntot, file_name, parameter_name, Nmin, Nmax, sampling, dx):
    print("plot fields 1d")
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
    a_group_key = list(file.keys())[0]
    data = list(file[a_group_key])
    str = '/data/'+data[ntot] +'/' + parameter_name
    V = np.array(file[str][:])
    
    Nx = V.shape[0]
    Ny = V.shape[1]
    
    V1 = V[:,int(Ny/2)]

    x = np.zeros([Nx])
    for i in range(Nx):
        x[i] = i*sampling*dx

    ax.plot(x[Nmin:Nmax], V1[Nmin:Nmax], linewidth=4)
    ax.set_xlabel(r'$X \omega_e/c$', fontsize=40,fontweight='bold')
    ax.set_ylabel(r'$B$', fontsize=40,fontweight='bold')
    ax.minorticks_on()
    plt.savefig('smilei_1d_'+ parameter_name +'.png', bbox_inches='tight')


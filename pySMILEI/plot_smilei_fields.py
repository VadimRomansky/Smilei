from matplotlib import animation
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np
import h5py
def plot_smilei_fields(ntot, file_name, parameter_name, xmin, xmax, ymin, ymax):
    f1 = plt.figure(figsize=[10,8])
    ax = f1.add_subplot(111)


    #cax1 = f1.add_axes([0.91, 0.12, 0.03, 0.75])
    cax2 = f1.add_axes([0.125, 0.92, 0.75, 0.03])

    file = h5py.File(file_name, 'r')
    print("Keys: %s" % file.keys())
    l = list(file.keys())
    a_group_key = list(file.keys())[0]
    data = list(file[a_group_key])
    str = '/data/'+data[ntot] +'/' + parameter_name
    V = np.array(file[str][:]).T
    minV = np.min(V);
    maxV = np.max(V);


    im2 = ax.imshow(V, origin='lower', aspect='auto',
                    extent=[xmin, xmax, ymin, ymax])  # plotting fluid data.
    im2.set_clim(minV, maxV)
    plt.colorbar(im2, cax=cax2, orientation='horizontal')  # vertical colorbar for fluid data.
    ax.set_xlabel(r'R', fontsize=18)
    ax.set_ylabel(r'Z', fontsize=18)
    ax.minorticks_on()
    plt.savefig('smilei_'+parameter_name+'.png')


from matplotlib import animation
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np
import h5py
def plot_smilei_fields_animated(ntot, file_name, parameter_name, xmin, xmax, ymin, ymax):
    print("plot fields animated")
    f1 = plt.figure(figsize=[10,8])
    ax = f1.add_subplot(111)


    #cax1 = f1.add_axes([0.91, 0.12, 0.03, 0.75])
    cax2 = f1.add_axes([0.125, 0.92, 0.75, 0.03])

    file = h5py.File(file_name, 'r')
    print("Keys: %s" % file.keys())
    l = list(file.keys())
    a_group_key = list(file.keys())[0]
    data = list(file[a_group_key])
    str = '/data/'+data[0] +'/' + parameter_name
    V = np.array(file[str][:]).T
    minV = V[0][0]
    maxV = V[0][0]
    for i in range(ntot):
        str = '/data/' + data[i] +'/' + parameter_name
        V = np.array(file[str][:]).T
        if(np.amin(V) < minV):
            minV = np.amin(V)
        if(np.amax(V) > maxV):
            maxV = np.amax(V)

    im2 = ax.imshow(V, origin='lower', aspect='auto',
                    extent=[xmin, xmax, ymin, ymax])  # plotting fluid data.
    im2.set_clim(minV, maxV)
    plt.colorbar(im2, cax=cax2, orientation='horizontal')  # vertical colorbar for fluid data.
    ax.set_xlabel(r'R', fontsize=18)
    ax.set_ylabel(r'Z', fontsize=18)
    ax.minorticks_on()

    def update(frame_number):
        print(frame_number)
        ax.clear()
        ax.set_xlim([xmin, xmax])
        ax.set_ylim([ymin, ymax])
        str = '/data/' + data[frame_number] +'/' + parameter_name
        V = np.array(file[str][:]).T
        im2 = ax.imshow(V, origin='upper', aspect = 'auto',
                        extent=[xmin, xmax, ymin, ymax])  # plotting fluid data.
        im2.set_clim(minV, maxV)
        plt.colorbar(im2, cax=cax2, orientation='horizontal')  # vertical colorbar for fluid data.
        #time.sleep(1)
        return im2

    anim = FuncAnimation(f1, update, interval=10, frames = ntot)

    #plt.show()

    f = r"smilei_"+parameter_name+".gif"
    writergif = animation.PillowWriter(fps=4)
    anim.save(f, writer=writergif)
    plt.close()

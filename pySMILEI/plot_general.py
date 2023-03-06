from matplotlib import animation
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np
import h5py

def plot_general(file_name, Nspecies):
    f1 = plt.figure(figsize=[10, 8])
    ax = f1.add_subplot(111)

    file = open(file_name, 'r')
    lines_list = file.readlines()
    Ntot = len(lines_list)
    Nhead = 0;
    for line in lines_list:
        if(line.split()[0]=="#"):
            Nhead = Nhead + 1

    Ntail = Ntot - Nhead
    energy = np.zeros([Ntail, 3 + Nspecies])
    time = np.zeros([Ntail])
    labels = ['Utot', 'Ukin', 'Uelm']
    for i in range(Nhead, Ntot):
        split_line = lines_list[i].split()
        time[i-Nhead] = float(split_line[0])
        energy[i-Nhead,0] = float(split_line[3]) #Utot
        energy[i-Nhead,1] = float(split_line[5]) #Ukin
        energy[i-Nhead,2] = float(split_line[8]) #Uelm
        for j in range(Nspecies):
            labels.append('species' + str(j))
            energy[i-Nhead,3+j] = float(split_line[19 + 5*j])

    ax.plot(time, energy)
    ax.set_xlabel(r't', fontsize=18)
    ax.set_ylabel(r'E', fontsize=18)
    ax.minorticks_on()
    ax.legend(labels)
    plt.savefig('smilei_general.png')

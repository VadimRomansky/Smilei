from matplotlib import animation
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np
import h5py
import matplotlib.ticker as mticker

def plot_general(file_name, prefix, Nspecies):
    print("plot general")
    plt.rcParams.update({'font.size': 40})
    plt.rcParams['text.usetex'] = True
    plt.rcParams['axes.linewidth'] = 4
    f1 = plt.figure(figsize=[10, 10])
    ax = f1.add_subplot(111)
    ax.tick_params(axis='x', size=10, width=4)
    ax.tick_params(axis='y', size=10, width=4)
    ax.minorticks_on()

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
    labels = [r'$U_{tot}$', r'$U_{kin}$', r'$U_{elm}$']
    for i in range(Nhead, Ntot):
        split_line = lines_list[i].split()
        time[i-Nhead] = float(split_line[0])
        energy[i-Nhead,0] = float(split_line[3]) #Utot
        energy[i-Nhead,1] = float(split_line[5]) #Ukin
        energy[i-Nhead,2] = float(split_line[8]) #Uelm
        for j in range(Nspecies):
            labels.append('$species_' + str(j)+"$")
            energy[i-Nhead,3+j] = float(split_line[19 + 5*j])


    U0 = energy[0,0]
    energy = energy/U0
    ax.plot(time, energy, linewidth=4)
    # after plotting the data, format the labels
    current_values = plt.gca().get_xticks()
    # using format string '{:.0f}' here but you can choose others
    ax.xaxis.set_major_locator(mticker.FixedLocator(current_values))
    new_labels = ['{:g}'.format(x) for x in current_values]
    plt.gca().set_xticklabels(new_labels)

    #ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: '{:g}'.format(x)))

    #ax.ticklabel_format(axis = 'x',style = 'sci')
    ax.set_xlabel(r'$t/\omega_e$', fontsize=40,fontweight='bold')
    ax.set_ylabel(r'$E/E_0$', fontsize=40,fontweight='bold')
    ax.minorticks_on()
    ax.legend(labels, fontsize=30)
    plt.savefig('smilei_general' + prefix + '.png', bbox_inches='tight')
    plt.close()

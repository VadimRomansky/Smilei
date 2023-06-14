import numpy as np
import h5py
import matplotlib.pyplot as plt
import numpy as np
import h5py
import matplotlib.ticker as mticker


def plot_shock_wave(ntot, file_name, treshold, dx, sampling, dt, Nt):
    file = h5py.File(file_name, 'r')
    l = list(file.keys())

    n = np.array(file.get(l[ntot - 1])).T
    Nx = n.shape[0]

    n0 = (np.array(file.get(l[0])).T)[1]
    time = np.zeros(ntot)
    for i in range(ntot):
        time[i] = i*dt*Nt
    swpoints = np.zeros([ntot])
    for j in range(ntot):
        swpoints[j] = Nx-1;
        n = np.array(file.get(l[j])).T
        for i in range(Nx-1,-1,-1):
            if n[i] > treshold*n0:
                swpoints[j] = i*dx*sampling
                break
    swpoints[0] = 0
    print("plot shock wave")
    plt.rcParams.update({'font.size': 40})
    plt.rcParams['text.usetex'] = True
    plt.rcParams['axes.linewidth'] = 4

    f1 = plt.figure(figsize=[10, 10])
    ax = f1.add_subplot(111)
    ax.tick_params(axis='x', size=10, width=4)
    ax.tick_params(axis='y', size=10, width=4)
    ax.minorticks_on()
    ax.plot(time, swpoints, linewidth=4)
    ax.set_xlabel(r'$t/\omega_e$', fontsize=40, fontweight='bold')
    ax.set_ylabel(r'$x_{sh} \omega_e/c$', fontsize=40, fontweight='bold')
    ax.minorticks_on()
    plt.savefig('smilei_shock_wave.png', bbox_inches='tight')
    plt.close()

    v = np.zeros(ntot-1)
    for i in range(ntot-1):
        v[i] = (swpoints[i+1] - swpoints[i])/(dt*Nt)
    f1 = plt.figure(figsize=[10, 10])
    ax = f1.add_subplot(111)
    ax.tick_params(axis='x', size=10, width=4)
    ax.tick_params(axis='y', size=10, width=4)
    ax.minorticks_on()
    ax.plot(time[0:ntot-1], v, linewidth=4)
    ax.set_xlabel(r'$t/\omega_e$', fontsize=40, fontweight='bold')
    ax.set_ylabel(r'$V_{sh}/c$', fontsize=40, fontweight='bold')
    ax.minorticks_on()
    plt.savefig('smilei_shock_velocity.png', bbox_inches='tight')
    plt.close()
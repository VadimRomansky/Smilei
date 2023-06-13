import numpy as np
import h5py


def get_sw_point(ntot, file_name):
    file = h5py.File(file_name, 'r')
    l = list(file.keys())

    n = np.array(file.get(l[ntot - 1])).T
    Nx = n.shape[0]

    n0 = (np.array(file.get(l[0])).T)[1]

    swpoint = Nx-1;
    for i in range(Nx-1,-1,-1):
        if n[i] > 1.5*n0:
            return i
    print('shock wave not found')
    return 100
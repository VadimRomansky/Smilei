from get_smilei_number import get_smilei_number
from plot_general import plot_general
from plot_smilei_concentration1d import plot_smilei_concentration1d
from plot_smilei_concentration1d_animated import plot_smilei_concentration1d_animated
from plot_smilei_concentration2d import plot_smilei_concentration2d
from plot_smilei_concentration2d_animated import plot_smilei_concentration2d_animated
from plot_smilei_fields_animated import plot_smilei_fields_animated
from plot_smilei_fields import plot_smilei_fields
from plot_smilei_spectrum import plot_smilei_spectrum
from plot_smilei_spectrum_animated import plot_smilei_spectrum_animated

dx = [0.1,0.1]
Npatches = [64, 32]
Ncells = [200, 50]
grid_length = [0,0]
samplingFactor = 20
for i in range(2):
    grid_length[i] = dx[i]*Npatches[i]*Ncells[i]


ntot = get_smilei_number("../ParticleBinning0.h5")
plot_general("../scalars.txt",4)
plot_smilei_fields_animated(ntot, "../Fields0.h5", "By", 0, grid_length[0], 0, grid_length[1])
plot_smilei_fields(ntot, "../Fields0.h5", "By",  0, grid_length[0], 0, grid_length[1])
plot_smilei_spectrum(ntot,"../ParticleBinning6.h5", 0, int(Npatches[0]*Ncells[0]/samplingFactor))
plot_smilei_spectrum_animated(ntot,"../ParticleBinning6.h5",0,int(Npatches[0]*Ncells[0]/samplingFactor))
plot_smilei_concentration2d(ntot,"../ParticleBinning3.h5", 0, grid_length[0], 0, grid_length[1])
plot_smilei_concentration2d_animated(ntot,"../ParticleBinning5.h5", 0, grid_length[0], 0, grid_length[1])
plot_smilei_concentration1d(ntot,"../ParticleBinning0.h5")
plot_smilei_concentration1d_animated(ntot,"../ParticleBinning2.h5")
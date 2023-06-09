from get_smilei_number import get_smilei_number
from plot_general import plot_general
from plot_smilei_concentration1d import plot_smilei_concentration1d
from plot_smilei_concentration1d_animated import plot_smilei_concentration1d_animated
from plot_smilei_concentration2d import plot_smilei_concentration2d
from plot_smilei_concentration2d_animated import plot_smilei_concentration2d_animated
from plot_smilei_fields_animated import plot_smilei_fields_animated
from plot_smilei_fields import plot_smilei_fields
from plot_smilei_fields1d_animated import plot_smilei_fields1d_animated
from plot_smilei_fields1d import plot_smilei_fields1d
from plot_smilei_spectrum import plot_smilei_spectrum
from plot_smilei_spectrum_animated import plot_smilei_spectrum_animated

dx = [0.2,0.2]
Npatches = [256, 16]
Ncells = [100, 100]
grid_length = [0,0]
samplingFactor = 20
for i in range(2):
    grid_length[i] = dx[i]*Npatches[i]*Ncells[i]

#path = "../output_shear_jet_gamma1.5_nj0.01_n1_sigma6/"
path = "../output/"
ntot = get_smilei_number(path + "ParticleBinning0.h5")-1
#ntot = 0
print("ntot = ",ntot)
plot_general(path + "scalars.txt","",2)

plot_smilei_fields_animated(ntot, path + "Fields0.h5", "By", 0, grid_length[0], 0, grid_length[1])
plot_smilei_fields_animated(ntot, path + "Fields0.h5", "Bz", 0, grid_length[0], 0, grid_length[1])
plot_smilei_fields(ntot, path + "Fields0.h5", "By",  0, grid_length[0], 0, grid_length[1])
plot_smilei_fields(ntot, path + "Fields0.h5", "Bz",  0, grid_length[0], 0, grid_length[1])

plot_smilei_fields1d_animated(ntot, path + "Fields0.h5", "By", int(0.0*Npatches[0]*Ncells[0]/samplingFactor), int(1.0*Npatches[0]*Ncells[0]/samplingFactor))
plot_smilei_fields1d_animated(ntot, path + "Fields0.h5", "Bz", int(0.0*Npatches[0]*Ncells[0]/samplingFactor), int(1.0*Npatches[0]*Ncells[0]/samplingFactor))
plot_smilei_fields1d(ntot, path + "Fields0.h5", "By", int(0.0*Npatches[0]*Ncells[0]/samplingFactor), int(1.0*Npatches[0]*Ncells[0]/samplingFactor))
plot_smilei_fields1d(ntot, path + "Fields0.h5", "Bz", int(0.0*Npatches[0]*Ncells[0]/samplingFactor), int(1.0*Npatches[0]*Ncells[0]/samplingFactor))

plot_smilei_spectrum(ntot,path + "ParticleBinning7.h5","protons", 0, int(0.25*Npatches[0]*Ncells[0]/samplingFactor))
plot_smilei_spectrum_animated(ntot,path + "ParticleBinning7.h5","protons", 100,0,int(0.25*Npatches[0]*Ncells[0]/samplingFactor))
plot_smilei_spectrum(ntot,path + "ParticleBinning6.h5","electrons", 0, int(0.25*Npatches[0]*Ncells[0]/samplingFactor))
plot_smilei_spectrum_animated(ntot,path + "ParticleBinning6.h5","electrons", 1,0,int(0.25*Npatches[0]*Ncells[0]/samplingFactor))

#plot_smilei_concentration2d(ntot,path + "ParticleBinning2.h5","electrons", 0, grid_length[0], 0, grid_length[1])
#plot_smilei_concentration2d_animated(ntot,path + "ParticleBinning2.h5","electrons", 0, grid_length[0], 0, grid_length[1])
#plot_smilei_concentration2d(ntot,path + "ParticleBinning5.h5","positrons", 0, grid_length[0], 0, grid_length[1])
#plot_smilei_concentration2d_animated(ntot,path + "ParticleBinning5.h5","positrons", 0, grid_length[0], 0, grid_length[1])
#plot_smilei_concentration2d(ntot,path + "ParticleBinning3.h4","protons", 0, grid_length[0], 0, grid_length[1])
#plot_smilei_concentration2d_animated(ntot,path + "ParticleBinning3.h5","protons", 0, grid_length[0], 0, grid_length[1])

plot_smilei_concentration1d(ntot,path + "ParticleBinning0.h5","electrons")
plot_smilei_concentration1d_animated(ntot,path + "ParticleBinning0.h5","electrons")
plot_smilei_concentration1d(ntot,path + "ParticleBinning1.h5","protons")
plot_smilei_concentration1d_animated(ntot,path + "ParticleBinning1.h5","protons")

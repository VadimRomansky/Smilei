from matplotlib import pyplot as plt

from get_smilei_number import get_smilei_number
from get_sw_point import get_sw_point
from plot_general import plot_general
from plot_shock_wave import plot_shock_wave
from plot_smilei_concentration1d import plot_smilei_concentration1d
from plot_smilei_concentration1d_animated import plot_smilei_concentration1d_animated
from plot_smilei_concentration2d import plot_smilei_concentration2d
from plot_smilei_concentration2d_animated import plot_smilei_concentration2d_animated
from plot_smilei_fields_animated import plot_smilei_fields_animated
from plot_smilei_fields import plot_smilei_fields
from plot_smilei_fields1d_animated import plot_smilei_fields1d_animated
from plot_smilei_fields1d import plot_smilei_fields1d
from plot_smilei_momentum_spectrum import plot_smilei_momentum_spectrum
from plot_smilei_momentum_spectrum_animated import plot_smilei_momentum_spectrum_animated
from plot_smilei_spectrum import plot_smilei_spectrum
from plot_smilei_spectrum_animated import plot_smilei_spectrum_animated

plt.set_cmap('jet')

dx = [0.2,0.2]
Npatches = [1024, 4]
Ncells = [250, 50]
grid_length = [0,0]
samplingPart = 20
samplingFields = 4
dt = 0.5*dx[0]
Nt = 20000

for i in range(2):
    grid_length[i] = dx[i]*Npatches[i]*Ncells[i]

#path = "../output_shear_jet_gamma1.5_nj0.01_n1_sigma6/"
path = "../output/"
#path = "../output_theta0-90_gamma1.5_sigma0.004/"
ntot = get_smilei_number(path + "ParticleBinning5.h5")
#ntot = 0
print("ntot = ",ntot)
#lot_general(path + "scalars.txt","",2)

#plot_smilei_fields_animated(ntot, path + "Fields03.h5", "Bx", 0, grid_length[0], 0, grid_length[1])
#plot_smilei_fields_animated(ntot, path + "Fields03.h5", "By", 0, grid_length[0], 0, grid_length[1])
#plot_smilei_fields_animated(ntot, path + "Fields03.h5", "Bz", 0, grid_length[0], 0, grid_length[1])
#plot_smilei_fields(ntot, path + "Fields03.h5", "Bx",  0, grid_length[0], 0, grid_length[1])
#plot_smilei_fields(ntot, path + "Fields03.h5", "By",  0, grid_length[0], 0, grid_length[1])
#plot_smilei_fields(ntot, path + "Fields03.h5", "Bz",  0, grid_length[0], 0, grid_length[1])

#plot_smilei_fields_animated(ntot, path + "Fields03.h5", "Ex", 0, grid_length[0], 0, grid_length[1])
#plot_smilei_fields_animated(ntot, path + "Fields03.h5", "Ey", 0, grid_length[0], 0, grid_length[1])
#plot_smilei_fields_animated(ntot, path + "Fields03.h5", "Ez", 0, grid_length[0], 0, grid_length[1])
#plot_smilei_fields(ntot, path + "Fields03.h5", "Ex",  0, grid_length[0], 0, grid_length[1])
#plot_smilei_fields(ntot, path + "Fields03.h5", "Ey",  0, grid_length[0], 0, grid_length[1])
#plot_smilei_fields(ntot, path + "Fields03.h5", "Ez",  0, grid_length[0], 0, grid_length[1])

#plot_smilei_fields1d_animated(ntot, path + "Fields0.h5", "Bx", int(0.0*Npatches[0]*Ncells[0]/samplingFields), int(1.0*Npatches[0]*Ncells[0]/samplingFields), samplingFields, dx[0])
#plot_smilei_fields1d_animated(ntot, path + "Fields0.h5", "By", int(0.0*Npatches[0]*Ncells[0]/samplingFields), int(1.0*Npatches[0]*Ncells[0]/samplingFields), samplingFields, dx[0])
#plot_smilei_fields1d_animated(ntot, path + "Fields0.h5", "Bz", int(0.0*Npatches[0]*Ncells[0]/samplingFields), int(1.0*Npatches[0]*Ncells[0]/samplingFields), samplingFields, dx[0])
#plot_smilei_fields1d(ntot, path + "Fields0.h5", "Bx", int(0.0*Npatches[0]*Ncells[0]/samplingFields), int(1.0*Npatches[0]*Ncells[0]/samplingFields), samplingFields, dx[0])
#plot_smilei_fields1d(ntot, path + "Fields0.h5", "By", int(0.0*Npatches[0]*Ncells[0]/samplingFields), int(1.0*Npatches[0]*Ncells[0]/samplingFields), samplingFields, dx[0])
#plot_smilei_fields1d(ntot, path + "Fields0.h5", "Bz", int(0.0*Npatches[0]*Ncells[0]/samplingFields), int(1.0*Npatches[0]*Ncells[0]/samplingFields), samplingFields, dx[0])

#plot_smilei_fields1d_animated(ntot, path + "Fields03.h5", "Ex", int(0.0*Npatches[0]*Ncells[0]/samplingFields), int(1.0*Npatches[0]*Ncells[0]/samplingFields), samplingFields, dx[0])
#plot_smilei_fields1d_animated(ntot, path + "Fields03.h5", "Ey", int(0.0*Npatches[0]*Ncells[0]/samplingFields), int(1.0*Npatches[0]*Ncells[0]/samplingFields), samplingFields, dx[0])
#plot_smilei_fields1d_animated(ntot, path + "Fields03.h5", "Ez", int(0.0*Npatches[0]*Ncells[0]/samplingFields), int(1.0*Npatches[0]*Ncells[0]/samplingFields), samplingFields, dx[0])
#plot_smilei_fields1d(ntot, path + "Fields03.h5", "Ex", int(0.0*Npatches[0]*Ncells[0]/samplingFields), int(1.0*Npatches[0]*Ncells[0]/samplingFields), samplingFields, dx[0])
#plot_smilei_fields1d(ntot, path + "Fields03.h5", "Ey", int(0.0*Npatches[0]*Ncells[0]/samplingFields), int(1.0*Npatches[0]*Ncells[0]/samplingFields), samplingFields, dx[0])
#plot_smilei_fields1d(ntot, path + "Fields03.h5", "Ez", int(0.0*Npatches[0]*Ncells[0]/samplingFields), int(1.0*Npatches[0]*Ncells[0]/samplingFields), samplingFields, dx[0])

minEe = 0.1;
maxEe = 5000;
minEp = 0.1;
maxEp = 5000;
minE = minEp;
maxE = maxEp;

#swpoint = get_sw_point(ntot, path + "ParticleBinning03.h5", 2.0)
#print('shock wave point = ', swpoint)
#plot_shock_wave(ntot, path + "ParticleBinning03.h5", 4.0, dx[0], 1, dt, Nt)
#spectrumStartX = int(0.5*swpoint/samplingPart)
#spectrumStartX = int(100/samplingPart)
#spectrumEndX = int(2*swpoint/samplingPart)
#spectrumEndX = int(0.25 * Npatches[0] * Ncells[0] / samplingPart)

#plot_smilei_spectrum(ntot,path + "ParticleBinning7.h5","protons", minEp, maxEp, spectrumStartX, spectrumEndX)
#plot_smilei_spectrum_animated(ntot,path + "ParticleBinning7.h5","protons", 100, minEp, maxEp, spectrumStartX, spectrumEndX)
#plot_smilei_spectrum(ntot, path + "ParticleBinning6.h5","electrons", minEe, maxEe, spectrumStartX, spectrumEndX)
#plot_smilei_spectrum_animated(ntot, path + "ParticleBinning6.h5","electrons", 1, minEe, maxEe, spectrumStartX, spectrumEndX)

#plot_smilei_momentum_spectrum(ntot,path + "ParticleBinning7.h5","protons", minEp, maxEp, spectrumStartX, spectrumEndX, 100)
#plot_smilei_momentum_spectrum_animated(ntot,path + "ParticleBinning7.h5","protons", minEp, maxEp, spectrumStartX, spectrumEndX, 100)
#plot_smilei_momentum_spectrum(ntot, path + "ParticleBinning6.h5","electrons", minEe, maxEe, spectrumStartX, spectrumEndX, 1)
#plot_smilei_momentum_spectrum_animated(ntot, path + "ParticleBinning6.h5","electrons", minEe, maxEe, spectrumStartX, spectrumEndX, 1)


#plot_smilei_concentration2d(ntot,path + "ParticleBinning3.h5","electrons", 0, grid_length[0], 0, grid_length[1])
plot_smilei_concentration2d_animated(ntot,path + "ParticleBinning2.h5","electrons", 0, grid_length[0], 0, grid_length[1])
#plot_smilei_concentration2d(ntot,path + "ParticleBinning5.h5","positrons", 0, grid_length[0], 0, grid_length[1])
#plot_smilei_concentration2d_animated(ntot,path + "ParticleBinning5.h5","positrons", 0, grid_length[0], 0, grid_length[1])
#plot_smilei_concentration2d(ntot,path + "ParticleBinning3.h5","protons", 0, grid_length[0], 0, grid_length[1])
#plot_smilei_concentration2d_animated(ntot,path + "ParticleBinning3.h5","protons", 0, grid_length[0], 0, grid_length[1])

plot_smilei_concentration1d(ntot,path + "ParticleBinning0.h5", 1, dx[0], "electrons")
#plot_smilei_concentration1d_animated(ntot,path + "ParticleBinning03.h5", 1, dx[0],"electrons")
#plot_smilei_concentration1d(ntot,path + "ParticleBinning13.h5", 1, dx[0], "protons")
#plot_smilei_concentration1d_animated(ntot,path + "ParticleBinning13.h5", 1, dx[0], "protons")
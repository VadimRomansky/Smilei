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

ntot = get_smilei_number("../output/Fields0.h5")
plot_general("../output/scalars.txt",4)
plot_smilei_fields_animated(ntot, "../output/Fields0.h5", "By", 0, 1000, 0, 1000)
plot_smilei_fields(ntot, "../output/Fields0.h5", "By", 0, 1000, 0, 1000)
plot_smilei_spectrum(ntot,"../output/ParticleBinning6.h5", 0, 160)
plot_smilei_spectrum_animated(ntot,"../output/ParticleBinning6.h5",0,160)
plot_smilei_concentration2d(ntot,"../output/ParticleBinning3.h5",0,1000,0,1000)
plot_smilei_concentration2d_animated(ntot,"../output/ParticleBinning5.h5",0,1000,0,1000)
plot_smilei_concentration1d(ntot,"../output/ParticleBinning0.h5")
plot_smilei_concentration1d_animated(ntot,"../output/ParticleBinning2.h5")
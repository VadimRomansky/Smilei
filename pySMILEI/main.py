from plot_smilei_concentration1d import plot_smilei_concentration1d
from plot_smilei_concentration1d_animated import plot_smilei_concentration1d_animated
from plot_smilei_concentration2d import plot_smilei_concentration2d
from plot_smilei_concentration2d_animated import plot_smilei_concentration2d_animated
from plot_smilei_fields_animated import plot_smilei_fields_animated
from plot_smilei_fields import plot_smilei_fields
from plot_smilei_spectrum import plot_smilei_spectrum
from plot_smilei_spectrum_animated import plot_smilei_spectrum_animated

plot_smilei_fields_animated(10, "../output/Fields0.h5", "By", 0, 1000, 0, 1000)
plot_smilei_fields(10, "../output/Fields0.h5", "By", 0, 1000, 0, 1000)
plot_smilei_spectrum(10,"../output/ParticleBinning6.h5", 0, 160)
plot_smilei_spectrum_animated(10,"../output/ParticleBinning6.h5",0,160)
plot_smilei_concentration2d(10,"../output/ParticleBinning3.h5",0,1000,0,1000)
plot_smilei_concentration2d_animated(10,"../output/ParticleBinning5.h5",0,1000,0,1000)
plot_smilei_concentration1d(10,"../output/ParticleBinning0.h5")
plot_smilei_concentration1d_animated(10,"../output/ParticleBinning2.h5")
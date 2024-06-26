# ----------------------------------------------------------------------------------------
# SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
#
# Magnetic shower in 3D with particle merging
#
# ----------------------------------------------------------------------------------------

import math as math
import numpy as np

# Physical parameters
c = 299792458
lambdar = 1e-6
wr = 2*math.pi*c/lambdar

# laser wavelength
l0 = 2.0*math.pi

# Schwinger electric field
Schwinger_E_field= 1.3E18
# Normalization electric field at lambda = 1
Enorm = 3.2E12

# Parameters for QED
chi = 10.
B_field_amplitude = 1000
gamma = chi* Schwinger_E_field/(Enorm*B_field_amplitude)
synchrotron_radius = math.sqrt(gamma**2 - 1.)/B_field_amplitude
velocity = math.sqrt(1.-1./gamma**2)

print("Gamma: {}".format(gamma))
print("Velocity: {}".format(velocity))

# Static B field properties
B_field_direction = np.array([0,0,1])
B_field_direction = B_field_direction / np.linalg.norm(B_field_direction)
B_field_vector = B_field_amplitude* B_field_direction

# Initial velocity vector for particles
mean_velocity = [0,velocity,0]
# Domain size
Lx = 1.*synchrotron_radius
Ly = 1.*synchrotron_radius
Lz = 1.*synchrotron_radius
# particle density
n0 = 1e-5
# nb of cells in one synchrotron radius
res = 32.
# Initial number of particles per cell
particles_per_cell = 2

dt_factor = 0.95
dx = synchrotron_radius/res                            # space step
dy = synchrotron_radius/res                            # space step
dz = synchrotron_radius/res                            # space step
dt  = 1./math.sqrt(1./(dx*dx) + 1./(dy*dy) + 1./(dz*dz)) # timestep (CFL)
dt *= dt_factor
# duration of the simulation
simulation_time = 50*dt

pusher = "vay"                         # dynamic type
radiation = "Monte-Carlo"

# Density profile for inital location of the particles
def n0_photon(x,y,z):
                return 0.


# Density profile for inital location of the particles
def n0_electron(x,y,z):
                return n0


# ______________________________________________________________________________
# Namelists

Main(
    geometry = "3Dcartesian",
    interpolation_order = 2 ,
    cell_length = [dx,dy,dz],
    grid_length  = [Lx,Ly,Lz],
    number_of_patches = [4,4,4],
    time_fields_frozen = simulation_time,
    timestep = dt,
    simulation_time = simulation_time,
    EM_boundary_conditions = [ ["periodic"] ],
    print_every = 10,
    reference_angular_frequency_SI = wr,
)

ExternalField(
    field = "Bx",
    profile = constant(B_field_vector[0])
)

ExternalField(
    field = "By",
    profile = constant(B_field_vector[1])
)

ExternalField(
    field = "Bz",
    profile = constant(B_field_vector[2])
)

Species(
    name = "electron",
    position_initialization = "random",
    momentum_initialization = "rectangular",
    particles_per_cell = particles_per_cell,
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    charge_density = n0_electron,
    mean_velocity = mean_velocity,
    temperature = [0.],
    pusher = pusher,
    radiation_model = "Monte-Carlo",
    radiation_photon_species = "photon",
    radiation_photon_gamma_threshold = 10,
    boundary_conditions = [
    	["periodic", "periodic"],
    	["periodic", "periodic"],
    	["periodic", "periodic"],
    ],
    
    # Merging parameters
    merging_method = "vranic_cartesian",
    merge_every = 2,
    merge_min_particles_per_cell = 16,
    merge_max_packet_size = 4,
    merge_min_packet_size = 4,
    merge_momentum_cell_size = [8,8,8],
    merge_accumulation_correction = False,
    merge_discretization_scale = "linear",
    # merge_min_momentum = 1e-5,
#    track_every = 1,
#    track_flush_every = 1000,
    is_test = False
)

Species(
    name = "positron",
    position_initialization = "random",
    momentum_initialization = "rectangular",
    particles_per_cell = particles_per_cell,
    c_part_max = 1.0,
    mass = 1.0,
    charge = 1.0,
    charge_density = n0_electron,
    mean_velocity = mean_velocity,
    temperature = [0.],
    pusher = pusher,
    radiation_model = "Monte-Carlo",
    radiation_photon_species = "photon",
    radiation_photon_gamma_threshold = 10,
    boundary_conditions = [
    	["periodic", "periodic"],
    	["periodic", "periodic"],
    	["periodic", "periodic"],
    ],

    # Merging parameters
    merging_method = "vranic_cartesian",
    merge_every = 2,
    merge_min_particles_per_cell = 16,
    merge_max_packet_size = 4,
    merge_min_packet_size = 4,
    merge_momentum_cell_size = [8,8,8],
    merge_accumulation_correction = False,
    merge_discretization_scale = "linear",
    # merge_min_momentum = 1e-5,
)

Species(
    name = "photon",
    position_initialization = "random",
    momentum_initialization = "cold",
    particles_per_cell = 0,
    c_part_max = 1.0,
    mass = 0,
    charge = 0.,
    number_density = n0_photon,
    mean_velocity = [0.0 ,0.0, 0.0],
    temperature = [0.],
    pusher = "norm",
    multiphoton_Breit_Wheeler = ["electron","positron"],
    multiphoton_Breit_Wheeler_sampling = [1,1],
    boundary_conditions = [
    	["periodic", "periodic"],
    	["periodic", "periodic"],
    	["periodic", "periodic"],
    ],
    
    # Merging parameters
    merging_method = "vranic_cartesian",
    merge_every = 2,
    merge_min_particles_per_cell = 8,
    merge_max_packet_size = 4,
    merge_min_packet_size = 4,
    merge_momentum_cell_size = [8,8,8],
    merge_accumulation_correction = False,
    merge_discretization_scale = "linear",
    # merge_min_momentum = 1e-5,
#    track_every = 1,
#    track_flush_every = 1000,
    is_test = False
)

Vectorization(
    mode = "on",
)

RadiationReaction(
    minimum_chi_discontinuous = 1e-2,
    # table_path = "./"
)

MultiphotonBreitWheeler(
    #table_path = "./"
)

DiagScalar(
    every = 1,
    vars=['Uelm','Ukin','Utot','Uexp','Ubal',
          'Urad',
          'UmBWpairs',
          'Ukin_electron',
          'Ukin_photon',
          'Ukin_positron',
          'Ntot_electron',
          'Ntot_photon',
          'Ntot_positron',
          'Dens_electron',
          'Dens_positron',
          'Dens_photon']
)

DiagPerformances(
    every = 1,
#    flush_every = 100,
#    patch_information = True,
)

species_list = ["electron","positron","photon"]

for species in species_list:

	DiagParticleBinning(
		deposited_quantity = "weight",
		every = 20,
		time_average = 1,
		species = [species],
		axes = [
		    ["gamma", 0, gamma, 128],
		]
	)

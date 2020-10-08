import math
import numpy as np

## ./smilei sim2d.py | tee "sim_log_$(date +'%F_%H-%M').txt"

# Mean velocity
mean_velocity = 0.75
# Electron temperature
T0 = 0.01
Te = T0
# Ion temperature
Ti = T0
# Ion charge
Zi = 1
# Density
n0 = 1
#Plasma parameter
beta = 0.25
#sigma
beta = 0.04
# Debye length
#Debye_length = 1. / np.sqrt( n0 / Te + Zi * n0 / Ti )
#ion mass
mp = 100.0
# Cell length
cell_length = [0.1, 0.1]
# Number of patches
number_of_patches = [1, 1]
# Cells per patches (patch shape)
cells_per_patch = [100., 100.]
# Grid length
grid_length = [10.,10.]
for i in range(2):
    grid_length[i] = number_of_patches[i] * cell_length[i] * cells_per_patch[i]
# Number of particles per cell
particles_per_cell = 4
# Position init
position_initialization = 'random'
# Timestep (Courant condition)
timestep = 0.45/np.sqrt(1./ cell_length[0]**2 + 1./ cell_length[1]**2)
# Total simulation timeremove
simulation_time = 100          # duration of the simulation
# Period of output for the diags
xdiag = 1000
vdiag = 1000
diag_step = 1000
#diag_step = 70 # after including the interrupt
diag_every = int(simulation_time / (diag_step*timestep))

Main(
    geometry = "2Dcartesian",
    interpolation_order = 2 ,
    cell_length = cell_length,
    grid_length  = grid_length,
    number_of_patches = number_of_patches,
    #cell_sorting = True,
    timestep = timestep,
    simulation_time = simulation_time,
    EM_boundary_conditions = [
        ['periodic'],
        ['periodic'],
    ],
    random_seed = smilei_mpi_rank,
)

# Initial plasma shape
fm = constant(n0)

Species(
	name = 'pon1',
	position_initialization = position_initialization,
	momentum_initialization = 'maxwell-juettner',
	ionization_model = 'none',
	particles_per_cell = particles_per_cell,
	#c_part_max = 1.0,
	pusher = 'vay',
	mass = mp,
	charge = 1.0,
	number_density = fm,
	mean_velocity = [-mean_velocity,0.,0.],
	temperature = [Ti],
	time_frozen = 0.0,
	boundary_conditions = [
		["periodic", "periodic"],
		["periodic", "periodic"],
	],
)

Species(
	name = 'eon1',
	position_initialization = position_initialization,
	momentum_initialization = 'maxwell-juettner',
	ionization_model = 'none',
	particles_per_cell = particles_per_cell,
	#c_part_max = 1.0,
	pusher = 'vay',
	mass = 1.0,
	charge = -1.0,
	number_density = fm,
	mean_velocity = [-mean_velocity,0.,0.],
	temperature = [Te],
	time_frozen = 0.0,
	boundary_conditions = [
		["periodic", "periodic"],
		["periodic", "periodic"],
	],
)

ExternalField(
   field = "By",
   profile = np.sqrt(0.04*4*mp*n0*1.5*3.14)
)

ExternalField(
   field = "Ez",
   profile = np.sqrt(0.04*4*mp*n0*1.5*3.14)*mean_velocity
)

CurrentFilter(
    model = "binomial",
    passes = [8],
    kernelFIR = [0.25,0.5,0.25]
)

FieldFilter(
    model = "Friedman",
    theta = 0.,
)

DiagScalar(every=1)

DiagFields(
    #0,
    every = diag_every,
    time_average = 1,
    fields = ["Ex", "Ey", "Ez", "Bx", "By", "Bz"],
    #subgrid = None
)

#number density via x coordinate distribution
DiagParticleBinning( #0
    deposited_quantity = "weight",
    every = diag_every,
    time_average = 1,
    species = ["eon1"],
    axes = [
        ["x", 0., grid_length[0], xdiag]#number_of_patches[0] * cell_length[0]], #128 is the number of x points 
    ]
)

    
DiagParticleBinning( #1
    deposited_quantity = "weight",
    every = diag_every,
    time_average = 1,
    species = ["pon1"],
    axes = [
        ["x", 0., grid_length[0], xdiag]#number_of_patches[0] * cell_length[0]], 
    ]
)

DiagParticleBinning( #2
    deposited_quantity = "weight",
    every = diag_every,
    time_average = 1,
    species = ["eon1"],
    axes = [
        ["ekin", 0.01, 10, 1000, "logscale"],
    ]
)

DiagParticleBinning( #3
    deposited_quantity = "weight",
    every = diag_every,
    time_average = 1,
    species = ["pon1"],
    axes = [
        ["ekin", 1, 1000, 1000, "logscale"],
    ]
)
    
##number density via kinetic energy density distribution "Temperature profile"
DiagParticleBinning( #4
    deposited_quantity = "weight_ekin",
    every = diag_every,
    time_average = 1,
    species = ["eon1"],
    axes = [
        ["x", 0., grid_length[0], xdiag],
    ]
)
    
DiagParticleBinning( #5
    deposited_quantity = "weight_ekin",
    every = diag_every,
    time_average = 1,
    species = ["pon1"],
    axes = [
        ["x", 0., grid_length[0], xdiag],
    ]
)


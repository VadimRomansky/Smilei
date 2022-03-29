import math
import numpy as np
from numpy import s_
## ./smilei sim2d.py | tee "sim_log_$(date +'%F_%H-%M').txt"

#gamma
gamma = 1.5
# Mean velocity
mean_velocity = np.sqrt(1.0-1.0/(gamma*gamma))
# Electron temperature
T0 = 0.1
Te = T0
# Ion temperature
Ti = T0
# Ion charge
Zi = 1
# Density
n0 = 1
#sigma
sigma = 0.004
#theta
theta=0.0*np.pi/180.0
# Debye length
#Debye_length = 1. / np.sqrt( n0 / Te + Zi * n0 / Ti )
#ion mass
mp = 100.0
# Cell length
cell_length = [0.2, 0.2]
# Number of patches
number_of_patches = [4096, 4]
# Cells per patches (patch shape)
cells_per_patch = [100., 100.]
# Grid length
grid_length = [81920.,80.]
for i in range(2):
    grid_length[i] = number_of_patches[i] * cell_length[i] * cells_per_patch[i]
# Number of particles per cell
particles_per_cell = 4
# Position init
position_initialization = 'random'
# Timestep (Courant condition)
timestep = 0.5*cell_length[0]
# Total simulation timeremove
simulation_time = 180000          # duration of the simulation
# Period of output for the diags
xdiag = 40
vdiag = 40
diag_step = 80000
#diag_step = 70 # after including the interrupt
diag_every = 80000

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
        ['reflective','silver-muller'],
        ['periodic'],
    ],
    random_seed = smilei_mpi_rank,
    solve_poisson = True,
    poisson_max_iteration = 10000,
    poisson_max_error = 1E-13
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
		["reflective", "remove"],
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
		["reflective", "remove"],
		["periodic", "periodic"],
	],
)

#ParticleInjector(
#    name      = "injector1",
#    species   = "pon1",
#    box_side  = "xmax",
#)

#ParticleInjector(
#    name      = "injector2",
#    species   = "eon1",
#    box_side  = "xmax",
#    position_initialization = "injector1"
#)

ExternalField(
   field = "Bx",
   profile = np.sqrt(sigma*mp*n0*gamma)*np.cos(theta)
)

ExternalField(
   field = "Bz",
   profile = np.sqrt(sigma*mp*n0*gamma)*np.sin(theta)
)

ExternalField(
   field = "Ey",
   profile = -np.sqrt(sigma*mp*n0*gamma)*np.sin(theta)*mean_velocity
)

CurrentFilter(
    model = "binomial",
    passes = [8],
    kernelFIR = [0.25,0.5,0.25]
)

FieldFilter(
    model = "Friedman",
    theta = 0.1,
)


DiagScalar(every=100)

#number density via x coordinate distribution
DiagParticleBinning( #0
    deposited_quantity = "weight",
    every = diag_every,
    time_average = 1,
    species = ["eon1"],
    axes = [
        ["x", 0., grid_length[0], number_of_patches[0] * cells_per_patch[0]]#number_of_patches[0] * cell_length[0]], #128 is the number of x points 
    ]
)

    
DiagParticleBinning( #1
    deposited_quantity = "weight",
    every = diag_every,
    time_average = 1,
    species = ["pon1"],
    axes = [
        ["x", 0., grid_length[0], number_of_patches[0] * cells_per_patch[0]], 
    ]
)

DiagParticleBinning( #2
    deposited_quantity = "weight",
    every = diag_every,
    time_average = 1,
    species = ["eon1"],
    axes = [
        ["ekin", 0.001, 5000, 1000, "logscale"],
    ]
)

DiagParticleBinning( #3
    deposited_quantity = "weight",
    every = diag_every,
    time_average = 1,
    species = ["pon1"],
    axes = [
        ["ekin", 0.1, 5000, 1000, "logscale"],
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

DiagParticleBinning(#6
    deposited_quantity = "weight",
    every = diag_every,
    time_average = 1,
    species = ["eon1"],
    axes = [
        ["x", 0., grid_length[0], number_of_patches[0] * cells_per_patch[0]/20],
        ["ekin", 0.001, 5000, 200, "logscale", "edge_inclusive"]
    ]
)

DiagParticleBinning(#7
    deposited_quantity = "weight",
    every = diag_every,
    time_average = 1,
    species = ["pon1"],
    axes = [
        ["x", 0., grid_length[0], number_of_patches[0] * cells_per_patch[0]/20],
        ["ekin", 0.1, 5000, 200, "logscale", "edge_inclusive"]
    ]
)

DiagFields(
    #8
    #name = "fields",
    every = diag_every,
    time_average = 1,
    fields = ["Ex", "Ey", "Ez","Bx", "By", "Bz"],
    subgrid = s_[::4, ::4]
)	

DiagParticleBinning(
    #9,
    deposited_quantity = "weight",
    every = diag_every,
    time_average = 1,
    species = ["pon1"],
    axes = [
        ["x", 0., grid_length[0], number_of_patches[0] * cells_per_patch[0]/20],
        ["vx", -1.0, 1.0, 200, "edge_inclusive"]
    ]
)

DiagParticleBinning(
    #10,
    deposited_quantity = "weight",
    every = diag_every,
    time_average = 1,
    species = ["eon1"],
    axes = [
        ["x", 0., grid_length[0], number_of_patches[0] * cells_per_patch[0]/20],
        ["vx", -1.0, 1.0, 200, "edge_inclusive"]
    ]
)


import math
import numpy as np
from numpy import s_
## ./smilei sim2d.py | tee "sim_log_$(date +'%F_%H-%M').txt"


#gamma
gamma = 10.0
# Mean velocity
mean_velocity = np.sqrt(1.0-1.0/(gamma*gamma))
# Electron temperature
T0 = 0.00001
Te = T0
# Ion temperature
Ti = T0
# Ion charge
Zi = 1
# Density
n0 = 1
#sigma
sigma = 0.001
#B0
B0 = np.sqrt(2*sigma*1.0*n0*gamma)
#theta
theta=90.0*np.pi/180.0
# Debye length
#Debye_length = 1. / np.sqrt( n0 / Te + Zi * n0 / Ti )
#ion mass
mp = 100.0
# Cell length
cell_length = [0.2, 0.2]
# Number of patches
number_of_patches = [256, 8]
# Cells per patches (patch shape)
cells_per_patch = [150., 100.]
# Grid length
grid_length = [7680.,160.]
for i in range(2):
    grid_length[i] = number_of_patches[i] * cell_length[i] * cells_per_patch[i]
# Number of particles per cell
particles_per_cell = 4
# Position init
position_initialization = 'random'
# Timestep (Courant condition)
timestep = 0.45* cell_length[0]
# Total simulation timeremove
simulation_time = 180000          # duration of the simulation
# Period of output for the diags
xdiag = 40
vdiag = 40
diag_step = 2000
#diag_step = 70 # after including the interrupt
diag_every = 2000

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
	mass = 1.0,
	charge = 1.0,
	number_density = fm,
	mean_velocity = [-mean_velocity,0.,0.],
	temperature = [Te],
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
   profile = B0*np.cos(theta)
)

ExternalField(
   field = "Bz",
   profile = B0*np.sin(theta)
)

ExternalField(
   field = "Ey",
   profile = -B0*np.sin(theta)*mean_velocity
)

CurrentFilter(
    model = "binomial",
    passes = [3],
    kernelFIR = [0.25,0.5,0.25]
)

FieldFilter(
    model = "Friedman",
    theta = 0.3,
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
        ["ekin", 0.01, 1000, 1000, "logscale"],
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

def ekin_moving(particles):
    vx = mean_velocity
    gam = 1.0/np.sqrt(1.0 - vx*vx)
    E0 = np.sqrt(1.0 + particles.px*particles.px + particles.py*particles.py + particles.pz*particles.pz)
    energy1 = E0*gam - particles.px*vx*gam - 1.0
    return energy1
	

DiagParticleBinning(#9
    deposited_quantity = ekin_moving,
    every = diag_every,
    time_average = 1,
    species = ["eon1"],
    axes = [
        ["x", 0., grid_length[0], number_of_patches[0] * cells_per_patch[0]/20],
	[ekin_moving, 0.001, 5000, 200, "logscale", "edge_inclusive"]

    ]
)
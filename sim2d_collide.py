import math
import numpy as np
from numpy import s_
## ./smilei sim2d.py | tee "sim_log_$(date +'%F_%H-%M').txt"

#gamma
gamma = 10.0
# Mean velocity
mean_velocity = np.sqrt(1.0-1.0/(gamma*gamma))
# Ion charge
Zi = 1
#right density
n2 = 1
# left density
n1 = 2.2
#sigma
sigma = 0.1
#theta
theta=60.0*np.pi/180.0
# Debye length
#Debye_length = 1. / np.sqrt( n1 / T1 + Zi * n2 / T2 )
#relativistic plasma frequency
omega_p = 1.0/np.sqrt(gamma)
# left temperature
T1 = 0.01
# right temperature
T2 = 0.001
#ion mass
me = 1.0
mp = 100.0
# Cell length
cell_length = [0.1, 0.1]
# Number of patches
number_of_patches = [1024, 1]
# Cells per patches (patch shape)
cells_per_patch = [200., 50.]
# Grid length
grid_length = [0.,0.]
for i in range(2):
    grid_length[i] = number_of_patches[i] * cell_length[i] * cells_per_patch[i]
# Number of particles per cell
particles_per_cell = 4
# Position init
position_initialization = 'random'
# Timestep (Courant condition)
timestep = 0.5*cell_length[0]
# Total simulation timeremove
simulation_time = 50000          # duration of the simulation
# Period of output for the diags
xdiag = 40
vdiag = 40
diag_step = 10000
#diag_step = 70 # after including the interrupt
diag_every = 10000

Bn = np.sqrt(sigma*me*n1*gamma)*np.cos(theta)
Bt1 = np.sqrt(sigma*me*n1*gamma)*np.sin(theta)
Et = np.sqrt(sigma*me*n1*gamma)*np.sin(theta)*mean_velocity

ratio = 0.5

def leftProfile(x,y):
    if(x < ratio*grid_length[0]):
        return n1
    else:
        return 0

def rightProfile(x,y):
    if(x > ratio*grid_length[0]):
        return n2
    else:
        return 0


def vel_Profile(x, y):
    if (x < ratio * grid_length[0]):
        return mean_velocity
    else:
        return 0

def By(x,y):
    if(x < ratio*grid_length[0]):
        return Bt1
    else:
        return -5*Bt1

def Ez(x,y):
    if(x < ratio*grid_length[0]):
        return -Bt1*mean_velocity
    else:
        return 0

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
        ['silver-muller','silver-muller'],
        ['periodic'],
    ],
    random_seed = smilei_mpi_rank,
    solve_poisson = True,
    poisson_max_iteration = 10000,
    poisson_max_error = 1E-13
)


Species(
	name = 'protons1',
	position_initialization = position_initialization,
	momentum_initialization = 'maxwell-juettner',
	ionization_model = 'none',
        particles_per_cell = particles_per_cell,
	#c_part_max = 1.0,
	pusher = 'vay',
	mass = mp,
	charge = 1.0,
	number_density = rightProfile,
	mean_velocity = [0.,0.,0.],
	temperature = [T2],
	time_frozen = 0.0,
	boundary_conditions = [
		["remove", "remove"],
		["periodic", "periodic"],
	],
)

Species(
	name = 'positrons1',
	position_initialization = position_initialization,
	momentum_initialization = 'maxwell-juettner',
	ionization_model = 'none',
        particles_per_cell = particles_per_cell,
	#c_part_max = 1.0,
	pusher = 'vay',
	mass = mp,
	charge = 1.0,
        number_density = leftProfile,
	mean_velocity = [mean_velocity,0.,0.],
	temperature = [T1],
	time_frozen = 0.0,
	boundary_conditions = [
		["remove", "remove"],
		["periodic", "periodic"],
	],
)

Species(
	name = 'electrons1',
	position_initialization = position_initialization,
	momentum_initialization = 'maxwell-juettner',
	ionization_model = 'none',
        particles_per_cell = particles_per_cell,
	#c_part_max = 1.0,
	pusher = 'vay',
	mass = 1.0,
	charge = -1.0,
	number_density = leftProfile,
	mean_velocity = [mean_velocity,0.,0.],
	temperature = [T1],
	time_frozen = 0.0,
	boundary_conditions = [
		["remove", "remove"],
		["periodic", "periodic"],
	],
)

Species(
	name = 'electrons2',
	position_initialization = position_initialization,
	momentum_initialization = 'maxwell-juettner',
	ionization_model = 'none',
        particles_per_cell = particles_per_cell,
	#c_part_max = 1.0,
	pusher = 'vay',
	mass = 1.0,
	charge = -1.0,
        number_density = rightProfile,
	mean_velocity = [0,0.,0.],
	temperature = [T2],
	time_frozen = 0.0,
	boundary_conditions = [
		["remove", "remove"],
		["periodic", "periodic"],
	],
)


ExternalField(
   field = "Bx",
   profile = Bn
)

ExternalField(
   field = "By",
   profile = By
)

ExternalField(
   field = "Ez",
   profile = Ez
)

CurrentFilter(
    model = "binomial",
    passes = [8],
    kernelFIR = [0.25,0.5,0.25]
)

FieldFilter(
    model = "Friedman",
    theta = 0.5,
)


DiagScalar(every=100)

#number density via x coordinate distribution
DiagParticleBinning( #0
    deposited_quantity = "weight",
    every = diag_every,
    time_average = 1,
    species = ["electrons1", "electrons2"],
    axes = [
        ["x", 0., grid_length[0], number_of_patches[0] * cells_per_patch[0]]#number_of_patches[0] * cell_length[0]], #128 is the number of x points 
    ]
)

DiagParticleBinning( #0
    deposited_quantity = "weight",
    every = diag_every,
    time_average = 1,
    species = ["protons1"],
    axes = [
        ["x", 0., grid_length[0], number_of_patches[0] * cells_per_patch[0]]#number_of_patches[0] * cell_length[0]], #128 is the number of x points 
    ]
)


DiagParticleBinning( #3
    deposited_quantity = "weight",
    every = diag_every,
    time_average = 1,
    species = ["positrons1"],
    axes = [
        ["x", 0., grid_length[0], number_of_patches[0] * cells_per_patch[0]], 
    ]
)
    
DiagParticleBinning(#6
    deposited_quantity = "weight",
    every = diag_every,
    time_average = 1,
    species = ["electrons1", "electrons2"],
    axes = [
        ["x", 0., grid_length[0], number_of_patches[0] * cells_per_patch[0]/20],
        ["ekin", 0.001, 50000, 200, "logscale", "edge_inclusive"]
    ]
)

DiagParticleBinning(#7
    deposited_quantity = "weight",
    every = diag_every,
    time_average = 1,
    species = ["protons1"],
    axes = [
        ["x", 0., grid_length[0], number_of_patches[0] * cells_per_patch[0]/20],
        ["ekin", 0.1, 50000, 200, "logscale", "edge_inclusive"]
    ]
)

DiagParticleBinning(#8
    deposited_quantity = "weight",
    every = diag_every,
    time_average = 1,
    species = ["positrons1"],
    axes = [
        ["x", 0., grid_length[0], number_of_patches[0] * cells_per_patch[0]/20],
        ["ekin", 0.001, 50000, 200, "logscale", "edge_inclusive"]
    ]
)


DiagFields(
    #8
    #name = "fields",
    every = diag_every,
    time_average = 1,
    fields = ["Ex", "Ey", "Ez","Bx", "By", "Bz"],
    subgrid = s_[::10, ::10]
)
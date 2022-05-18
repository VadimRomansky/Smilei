import math
import numpy as np
from numpy import s_
## ./smilei sim2d.py | tee "sim_log_$(date +'%F_%H-%M').txt"

# Mean velocity
mean_velocity = 0.75
#gamma
gamma = 1.0/np.sqrt(1.0 - mean_velocity*mean_velocity)
# Electron temperature
T1 = 0.001
Te = T1
# Ion temperature
Ti = T1
# Ion charge
Zi = 1
# Density
n0 = 1
#sigma
sigma = 0.000004
#theta
theta=90.0*np.pi/180.0
# Debye length
#Debye_length = 1. / np.sqrt( n0 / Te + Zi * n0 / Ti )
#ion mass
me = 1.0
mp = 100.0
# Cell length
cell_length = [0.2, 0.2]
# Number of patches
number_of_patches = [128, 4]
# Cells per patches (patch shape)
cells_per_patch = [200., 50.]
# Grid length
grid_length = [5120.,40.]
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

Npart = int(particles_per_cell*number_of_patches[0]*number_of_patches[1]*cells_per_patch[0]*cells_per_patch[1])
NpartLeft = int(Npart*0.8)
coords1 = np.zeros([2+1, Npart - NpartLeft])
coords2 = np.zeros([2+1, NpartLeft])
weight = n0*cell_length[0]*cell_length[1]/particles_per_cell
for i in range(NpartLeft):
    coords2[0][i] = np.random.rand()*grid_length[0]*0.8
    coords2[1][i] = np.random.rand()*grid_length[1]
    coords2[2][i] = weight
for i in range(Npart-NpartLeft):
    coords1[0][i] = (1.0 + 0.25*np.random.rand())*grid_length[0]*0.8
    coords1[1][i] = np.random.rand()*grid_length[1]
    coords1[2][i] = weight

Bn = np.sqrt(sigma*mp*n0*gamma)*np.cos(theta)
Bt1 = np.sqrt(sigma*mp*n0*gamma)*np.sin(theta)
Et = np.sqrt(sigma*mp*n0*gamma)*np.sin(theta)*mean_velocity

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

# Initial plasma shape
fm = constant(n0)

Species(
	name = 'protons1',
	position_initialization = coords1,
	momentum_initialization = 'maxwell-juettner',
	ionization_model = 'none',
	#c_part_max = 1.0,
	pusher = 'vay',
	mass = mp,
	charge = 1.0,
	mean_velocity = [0.,0.,0.],
	temperature = [T1],
	time_frozen = 0.0,
	boundary_conditions = [
		["remove", "remove"],
		["periodic", "periodic"],
	],
)

Species(
	name = 'positrons',
	position_initialization = coords2,
	momentum_initialization = 'maxwell-juettner',
	ionization_model = 'none',
	#c_part_max = 1.0,
	pusher = 'vay',
	mass = 1.0,
	charge = 1.0,
	mean_velocity = [mean_velocity,0.,0.],
	temperature = [0.2],
	time_frozen = 0.0,
	boundary_conditions = [
		["remove", "remove"],
		["periodic", "periodic"],
	],
)

Species(
	name = 'electrons1',
	position_initialization = coords1,
	momentum_initialization = 'maxwell-juettner',
	ionization_model = 'none',
	#c_part_max = 1.0,
	pusher = 'vay',
	mass = 1.0,
	charge = -1.0,
	mean_velocity = [0.,0.,0.],
	temperature = [T1],
	time_frozen = 0.0,
	boundary_conditions = [
		["remove", "remove"],
		["periodic", "periodic"],
	],
)

Species(
	name = 'electrons2',
	position_initialization = coords2,
	momentum_initialization = 'maxwell-juettner',
	ionization_model = 'none',
	#c_part_max = 1.0,
	pusher = 'vay',
	mass = 1.0,
	charge = -1.0,
	mean_velocity = [mean_velocity,0.,0.],
	temperature = [0.2],
	time_frozen = 0.0,
	boundary_conditions = [
		["remove", "remove"],
		["periodic", "periodic"],
	],
)

ParticleInjector(
    name      = "injector1",
    species   = "positrons",
    box_side  = "xmin",
)

ParticleInjector(
    name      = "injector2",
    species   = "electrons2",
    box_side  = "xmin",
    position_initialization = "injector1"
)

def By(x,y):
    if(x < grid_length[0]/2):
        return Bt2
    else:
        return Bt1

def Ey(x,y):
    if(x < grid_length[0]/2):
        return Bt1*mean_velocity
    else:
        return 0

ExternalField(
   field = "Bx",
   profile = Bn
)

ExternalField(
   field = "Bz",
   profile = Bt1
)

ExternalField(
   field = "Ey",
   profile = Ey
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
    species = ["electrons1", "electrons2"],
    axes = [
        ["x", 0., grid_length[0], number_of_patches[0] * cells_per_patch[0]]#number_of_patches[0] * cell_length[0]], #128 is the number of x points 
    ]
)

    
DiagParticleBinning( #1
    deposited_quantity = "weight",
    every = diag_every,
    time_average = 1,
    species = ["protons1"],
    axes = [
        ["x", 0., grid_length[0], number_of_patches[0] * cells_per_patch[0]], 
    ]
)

DiagParticleBinning( #2
    deposited_quantity = "weight",
    every = diag_every,
    time_average = 1,
    species = ["electrons1", "electrons2"],
    axes = [
        ["ekin", 0.01, 1000, 1000, "logscale"],
    ]
)

DiagParticleBinning( #3
    deposited_quantity = "weight",
    every = diag_every,
    time_average = 1,
    species = ["positrons"],
    axes = [
        ["x", 0., grid_length[0], number_of_patches[0] * cells_per_patch[0]], 
    ]
)
    
##number density via kinetic energy density distribution "Temperature profile"
DiagParticleBinning( #4
    deposited_quantity = "weight_ekin",
    every = diag_every,
    time_average = 1,
    species = ["electrons1", "electrons2"],
    axes = [
        ["x", 0., grid_length[0], xdiag],
    ]
)
    
DiagParticleBinning( #5
    deposited_quantity = "weight_ekin",
    every = diag_every,
    time_average = 1,
    species = ["protons1"],
    axes = [
        ["x", 0., grid_length[0], xdiag],
    ]
)

DiagParticleBinning(#6
    deposited_quantity = "weight",
    every = diag_every,
    time_average = 1,
    species = ["electrons1", "electrons2"],
    axes = [
        ["x", 0., grid_length[0], number_of_patches[0] * cells_per_patch[0]/20],
        ["ekin", 0.001, 1000, 200, "logscale", "edge_inclusive"]
    ]
)

DiagParticleBinning(#7
    deposited_quantity = "weight",
    every = diag_every,
    time_average = 1,
    species = ["protons1"],
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

DiagParticleBinning(#9
    deposited_quantity = "weight",
    every = diag_every,
    time_average = 1,
    species = ["positrons"],
    axes = [
        ["x", 0., grid_length[0], number_of_patches[0] * cells_per_patch[0]/20],
        ["ekin", 0.001, 1000, 200, "logscale", "edge_inclusive"]
    ]
)
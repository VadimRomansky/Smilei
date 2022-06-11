import math
import numpy as np
from numpy import s_
## ./smilei sim2d.py | tee "sim_log_$(date +'%F_%H-%M').txt"

#gamma
gamma = 100.0
# Mean velocity
mean_velocity = np.sqrt(1.0-1.0/(gamma*gamma))
# Ion charge
Zi = 1
#right density
n2 = 1.0/2.2
# left density
n1 = 1.0
#sigma
sigma = 10
#theta
theta=90.0*np.pi/180.0
# Debye length
#Debye_length = 1. / np.sqrt( n1 / T1 + Zi * n2 / T2 )
#relativistic plasma frequency
omega_p = 1.0/np.sqrt(gamma)
#eta
eta = 3
#left hot
T3 = sigma/(2.0*eta)
# left cold
T1 = 0.01
# right temperature
T2 = 0.01
#lambda
lambd = 640*np.sqrt(gamma)
#alpha
alpha = 0.1
#Delta
Delta = 1.0*np.sqrt(gamma)
#delta
delta = 2*np.pi*Delta/lambd
#ion mass
me = 1.0
mp = 100.0
# Cell length
cell_length = [0.05, 0.05]
# Number of patches
number_of_patches = [1024, 1]
# Cells per patches (patch shape)
cells_per_patch = [300., 20.]
# Grid length
grid_length = [30720.,2.]
for i in range(2):
    grid_length[i] = number_of_patches[i] * cell_length[i] * cells_per_patch[i]
# Number of particles per cell
particles_per_cell = 4
# Position init
position_initialization = 'regular'
# Timestep (Courant condition)
timestep = 0.5* cell_length[0]
# Total simulation timeremove
simulation_time = 30000          # duration of the simulation
# Period of output for the diags
xdiag = 40
vdiag = 40
diag_step = 20000
diag_every = 20000

Bn = np.sqrt(sigma*me*n1*gamma)*np.cos(theta)
Bt1 = np.sqrt(sigma*me*n1*gamma)*np.sin(theta)
Et = np.sqrt(sigma*me*n1*gamma)*np.sin(theta)*mean_velocity

betah = (np.sqrt(sigma)/(eta*gamma))/(omega_p*Delta)

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
fm = constant(n1)
ratio = 0.5

def dzeta(x):
    return (alpha + np.cos(2*np.pi*x/lambd))/delta

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

def hotProfile(x,y):
    if (x < ratio * grid_length[0]):
        return eta*n1/np.square(np.cosh(dzeta(x)))
    else:
        return 0

def velz_e_Profile(x,y):
    if (x < ratio * grid_length[0]):
        return (betah/gamma)*np.sign(np.sin(2*np.pi*x/lambd))
    else:
        return 0

def velz_p_Profile(x, y):
    if (x < ratio * grid_length[0]):
        return -(betah/gamma) * np.sign(np.sin(2 * np.pi * x / lambd))
    else:
        return 0

def By(x,y):
    if(x < ratio*grid_length[0]):
        return Bt1*np.tanh(dzeta(x))
    else:
        return Bt1

def Ez(x,y):
    if(x < ratio*grid_length[0]):
        return -Bt1*np.tanh(dzeta(x))*mean_velocity
    else:
        return 0

Npart = int(particles_per_cell*number_of_patches[0]*number_of_patches[1]*cells_per_patch[0]*cells_per_patch[1])
NpartLeft = int(Npart*ratio)
#coordse1 = np.zeros([2+1, Npart - NpartLeft])
#coordsp = np.zeros([2+1, Npart - NpartLeft])

#coordse2 = np.zeros([2+1, NpartLeft])
#coordspos = np.zeros([2+1, NpartLeft])
#weight = n0*cell_length[0]*cell_length[1]/particles_per_cell
#for i in range(NpartLeft):
#    coordse2[0][i] = np.random.rand()*grid_length[0]*0.5
#    coordse2[1][i] = np.random.rand()*grid_length[1]
#    coordse2[2][i] = weight
#    coordspos[0][i] = coordse2[0][i]
#    coordspos[1][i] = coordse2[1][i]
#    coordspos[2][i] = weight

#for i in range(Npart-NpartLeft):
#    coordse1[0][i] = (1.0 + np.random.rand())*grid_length[0]*0.5
#    coordse1[1][i] = np.random.rand()*grid_length[1]
#    coordse1[2][i] = weight
#    coordsp[0][i] = coordse1[0][i]
#    coordsp[1][i] = coordse1[1][i]
#    coordsp[2][i] = weight


Species(
	name = 'protons1',
	#position_initialization = coordsp,
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
	name = 'positrons_cold',
	#position_initialization = coordspos,
        position_initialization = position_initialization,
	momentum_initialization = 'maxwell-juettner',
	ionization_model = 'none',
	particles_per_cell = particles_per_cell,
	#c_part_max = 1.0,
	pusher = 'vay',
	mass = 1.0,
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
	name = 'positrons_hot',
	#position_initialization = coordspos,
        position_initialization = position_initialization,
	momentum_initialization = 'maxwell-juettner',
	ionization_model = 'none',
	particles_per_cell = particles_per_cell,
	#c_part_max = 1.0,
	pusher = 'vay',
	mass = 1.0,
	charge = 1.0,
	number_density = hotProfile,
	mean_velocity = [mean_velocity,0.,velz_p_Profile],
	temperature = [T3],
	time_frozen = 0.0,
	boundary_conditions = [
		["remove", "remove"],
		["periodic", "periodic"],
	],
)

Species(
	name = 'electrons_right',
	#position_initialization = coordse1,
        position_initialization = 'protons1',
	momentum_initialization = 'maxwell-juettner',
	ionization_model = 'none',
	particles_per_cell = particles_per_cell,
	#c_part_max = 1.0,
	pusher = 'vay',
	mass = 1.0,
	charge = -1.0,
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
	name = 'electrons_hot',
	#position_initialization = coordse2,
        position_initialization = 'positrons_hot',
	momentum_initialization = 'maxwell-juettner',
	ionization_model = 'none',
	particles_per_cell = particles_per_cell,
	#c_part_max = 1.0,
	pusher = 'vay',
	mass = 1.0,
	charge = -1.0,
	number_density = hotProfile,
	mean_velocity = [mean_velocity,0.,velz_e_Profile],
	temperature = [T3],
	time_frozen = 0.0,
	boundary_conditions = [
		["remove", "remove"],
		["periodic", "periodic"],
	],
)

Species(
	name = 'electrons_cold',
	#position_initialization = coordse2,
        position_initialization = 'positrons_cold',
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
    theta = 0.1,
)


DiagScalar(every=100)

#number density via x coordinate distribution
DiagParticleBinning( #0
    deposited_quantity = "weight",
    every = diag_every,
    time_average = 1,
    species = ["electrons_right", "electrons_hot","electrons_cold"],
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
    species = ["electrons_right", "electrons_hot","electrons_cold"],
    axes = [
        ["ekin", 0.001, 50000, 1000, "logscale"],
    ]
)

DiagParticleBinning( #3
    deposited_quantity = "weight",
    every = diag_every,
    time_average = 1,
    species = ["positrons_hot","positrons_cold"],
    axes = [
        ["x", 0., grid_length[0], number_of_patches[0] * cells_per_patch[0]], 
    ]
)
    
##number density via kinetic energy density distribution "Temperature profile"
DiagParticleBinning( #4
    deposited_quantity = "weight_ekin",
    every = diag_every,
    time_average = 1,
    species = ["electrons_right", "electrons_hot","electrons_cold"],
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
    species = ["electrons_right", "electrons_hot","electrons_cold"],
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
    species = ["positrons_hot","positrons_cold"],
    axes = [
        ["x", 0., grid_length[0], number_of_patches[0] * cells_per_patch[0]/20],
        ["ekin", 0.001, 50000, 200, "logscale", "edge_inclusive"]
    ]
)

DiagParticleBinning(
    #10,
    deposited_quantity = "weight",
    every = diag_every,
    time_average = 1,
    species = ["protons1"],
    axes = [
        ["x", 0., grid_length[0], number_of_patches[0] * cells_per_patch[0]/20],
        ["vx", -1.0, 1.0, 200, "edge_inclusive"]
    ]
)

DiagParticleBinning(
    #11,
    deposited_quantity = "weight",
    every = diag_every,
    time_average = 1,
    species = ["electrons_right", "electrons_hot","electrons_cold"],
    axes = [
        ["x", 0., grid_length[0], number_of_patches[0] * cells_per_patch[0]/20],
        ["vx", -1.0, 1.0, 200, "edge_inclusive"]
    ]
)
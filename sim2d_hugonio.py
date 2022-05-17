import math
import numpy as np
from numpy import s_
## ./smilei sim2d.py | tee "sim_log_$(date +'%F_%H-%M').txt"

# Mean velocity
mean_velocity = 0.3
#gamma
gamma = 1.0/np.sqrt(1.0 - mean_velocity*mean_velocity)
# Electron temperature
T1 = 0.00001
Te = T1
# Ion temperature
Ti = T1
# Ion charge
Zi = 1
# Density
n0 = 1
#sigma
sigma = 0.0002
#theta
theta=10.0*np.pi/180.0
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
cells_per_patch = [100., 50.]
# Grid length
grid_length = [2560.,40.]
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

ratio = 3.9
Npart = int(particles_per_cell*number_of_patches[0]*number_of_patches[1]*cells_per_patch[0]*cells_per_patch[1]*(1 + ratio)*0.5)
NpartLeft = int(Npart*ratio/(1.0 + ratio))
coords1 = np.zeros([2+1, Npart - NpartLeft])
coords2 = np.zeros([2+1, NpartLeft])
weight = n0*cell_length[0]*cell_length[1]/particles_per_cell
for i in range(NpartLeft):
    coords2[0][i] = np.random.rand()*grid_length[0]*0.5
    coords2[1][i] = np.random.rand()*grid_length[1]
    coords2[2][i] = weight
for i in range(Npart-NpartLeft):
    coords1[0][i] = (1.0 + np.random.rand())*grid_length[0]*0.5
    coords1[1][i] = np.random.rand()*grid_length[1]
    coords1[2][i] = weight

Bn = np.sqrt(sigma*mp*n0*gamma)*np.cos(theta)
Bt1 = np.sqrt(sigma*mp*n0*gamma)*np.sin(theta)
Et = np.sqrt(sigma*mp*n0*gamma)*np.sin(theta)*mean_velocity

densflux = (mp + me)*n0*mean_velocity
Bt2 = Bt1*(1 - Bn*Bn*n0*(mp+me)/(4*np.pi*densflux*densflux))/(1.0/ratio - Bn*Bn*n0*(mp+me)/(4*np.pi*densflux*densflux))

print(Bt1, Bt2)

adiab = 5.0/3.0
P1 = 2*n0*T1*me
P2 = ((1.0/(adiab - 1) - 0.5*(1.0/ratio - 1.0))*P1 - (1.0/ratio - 1.0)*(Bt2 - Bt1)*(Bt2 - Bt1)/(16*np.pi))/(1.0/(ratio*(adiab - 1.0)) + 0.5*(1.0/ratio - 1.0))
print(P1,P2)
T2 = T1*(P2/P1)/ratio	

print(T1, T2)

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
	mean_velocity = [-mean_velocity,0.,0.],
	temperature = [T1],
	time_frozen = 0.0,
	boundary_conditions = [
		["remove", "remove"],
		["periodic", "periodic"],
	],
)

Species(
	name = 'protons2',
	position_initialization = coords2,
	momentum_initialization = 'maxwell-juettner',
	ionization_model = 'none',
	#c_part_max = 1.0,
	pusher = 'vay',
	mass = mp,
	charge = 1.0,
	mean_velocity = [-mean_velocity/ratio,0.,0.],
	temperature = [T2],
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
	mean_velocity = [-mean_velocity,0.,0.],
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
	mean_velocity = [-mean_velocity/ratio,0.,0.],
	temperature = [T2],
	time_frozen = 0.0,
	boundary_conditions = [
		["remove", "remove"],
		["periodic", "periodic"],
	],
)

def By(x,y):
    if(x < grid_length[0]/2):
        return Bt2
    else:
        return Bt1

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
   profile = Et
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

    
DiagParticleBinning( #1
    deposited_quantity = "weight",
    every = diag_every,
    time_average = 1,
    species = ["protons1","protons2"],
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
    species = ["protons1","protons2"],
    axes = [
        ["ekin", 1, 5000, 1000, "logscale"],
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
    species = ["protons1","protons2"],
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
    species = ["protons1","protons2"],
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
    vx = -mean_velocity
    gam = 1.0/np.sqrt(1.0 - vx*vx)
    E0 = np.sqrt(1.0 + particles.px*particles.px + particles.py*particles.py + particles.pz*particles.pz)
    energy1 = E0*gam - particles.px*vx*gam - 1.0
    return energy1
	

DiagParticleBinning(#9
    deposited_quantity = ekin_moving,
    every = diag_every,
    time_average = 1,
    species = ["electrons1", "electrons2"],
    axes = [
        ["x", 0., grid_length[0], number_of_patches[0] * cells_per_patch[0]/20],
	[ekin_moving, 0.001, 1000, 200, "logscale", "edge_inclusive"]

    ]
)
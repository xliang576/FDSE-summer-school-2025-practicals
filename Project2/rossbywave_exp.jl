# A script to explore Rossby wave propagation and nonlinear behavior using Oceananigans

# Load some packages that we will need
using Oceananigans
# using Oceanostics
using Printf
using Oceananigans.Units
using JLD2

# Set the domain size in dimensional coordinates
Lx = 1000kilometers  
Ly = 1000kilometers

# Set the grid size
Nx = 128   # number of grid points in x-direction
Ny = 128   # number of grid points in y-direction

# Construct a rectilinear grid that is periodic in x-direction and bounded in y-direction
grid = RectilinearGrid(size = (Nx, Ny),
                       x = (0, Lx), y = (0, Ly),
                       topology = (Bounded, Bounded, Flat)
)

# Set up a model for Rossby waves
model = NonhydrostaticModel(; grid,
              advection = UpwindBiased(),   # Specify the advection scheme.  Another good choice is WENO() which is more accurate but slower
            timestepper = :RungeKutta3,   # Set the timestepping scheme, here 3rd order Runge-Kutta
                tracers = :c,
                coriolis = BetaPlane(rotation_rate = 7.292115e-5, latitude = 20, radius = 6371e3)   # set Coriolis parameter using the Beta-plane approximation 
)

# Set wavenumbers associated with the initial condition
k = 2 * pi / 200kilometers 
l = 2 * pi / 200kilometers

# Define functions for the initial conditions
u₀ = 0.01   # units: m/s
U = 0  # add a zonal mean background flow
uᵢ(x, y) = u₀ * sin(k * x) * sin(l * y) + U
vᵢ(x, y) = u₀ * (k / l) * cos(k * x) * cos(l * y)
wᵢ(x, y) = 0
cᵢ(x, y) = x # Here, we set the function for c so that it is proportional to the streamfunction associated with (u,v)

# Send the initial conditions to the model to initialize the variables
set!(model, u = uᵢ, v = vᵢ, w = wᵢ, c = cᵢ)
ϵ = KineticEnergyDissipationRate(model)
# Create a 'simulation' to run the model for a specified length of time
simulation = Simulation(model, Δt = 10hours, stop_iteration = 6000)

# Add callback that prints progress message during simulation
progress(sim) = @info string("Iter: ", iteration(sim),
                             ", time: ", prettytime(sim))

simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

# Save output from the simulation
filename = "rossbywave_bounded"

u, v, w = model.velocities
c = model.tracers.c
u_prime = u - Field(Average(u, dims = 1))
v_prime = v - Field(Average(v, dims = 1))
uv_prime = Field(Average(u_prime*v_prime,dims = 1)) # has some problems here

simulation.output_writers[:jld2] = JLD2Writer(model, (; u, v, w, c, uv_prime),
                                                    schedule = IterationInterval(50),
                                                    filename = filename * ".jld2",
                                                    overwrite_existing = true,
                                                    with_halos = false
)

# Run the simulation                                                  
run!(simulation)

# Make a plot of u at y=Ly/2 and save a movie
include("plot_rossbywave_exp.jl")

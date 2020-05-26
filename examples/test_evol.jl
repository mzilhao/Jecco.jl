
using OrdinaryDiffEq

using Jecco
using Jecco.AdS5_3_1

grid = SpecCartGrid3D(
    x_min            = -5.0,
    x_max            =  5.0,
    x_nodes          =  128,
    y_min            = -5.0,
    y_max            =  5.0,
    y_nodes          =  128,
    u_outer_min      =  0.1,
    u_outer_max      =  1.0,
    # u_outer_domains  =  1,
    # u_outer_nodes    =  64,
    u_outer_domains  =  2,
    u_outer_nodes    =  16,
    u_inner_nodes    =  12,
)

id   = BlackBrane()

# evoleq = EvolTest0()

evoleq = AffineNull(
    phi0       = 0.0,
    potential  = ZeroPotential(),
)

# atlas of grid configuration and respective SystemPartition
atlas     = Atlas(grid)
systems   = SystemPartition(grid)

# allocate variables
boundary       = Boundary(grid)
gauge          = Gauge(grid)
bulkevols      = BulkEvolvedPartition(grid)
bulkconstrains = BulkConstrainedPartition(grid)

# initial conditions
init_data!(bulkevols, boundary, gauge, systems, id)

# full state vector
evolvars  = EvolVars(boundary, gauge, bulkevols)

# function that updates the state vector
rhs! = AdS5_3_1.setup_rhs(bulkconstrains, systems)


dt0   = 0.001
tspan = (0.0, 0.01)

prob  = ODEProblem(rhs!, evolvars, tspan, evoleq)

# integrator = init(prob, RK4(), save_everystep=false, dt=dt0, adaptive=false)
integrator = init(prob, AB3(), save_everystep=false, dt=dt0, adaptive=false)

for (f,t) in tuples(integrator)
    B1  = AdS5_3_1.getB1(f, 1)
    a4  = AdS5_3_1.geta4(f)
    @show t, B1[1], a4[1]
end

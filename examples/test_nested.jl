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

# id = AdS5_3_1.IDTest0(
#     b14_0      = 0.01,
#     b24_0      = 0.02,
#     phi0       = 0.0,
#     phi2_0     = 0.01,
#     a4_0       = -1.0,
#     fx2_0      = 0.02,
#     fy2_0      = 0.1,
#     potential  = ZeroPotential(),
# )

id = BlackBrane()

evoleq = AffineNull(
    phi0       = 0.0,
    potential  = ZeroPotential(),
)

# atlas of grid configuration and respective SystemPartition
atlas     = Atlas(grid)
systems   = SystemPartition(grid)

# allocate variables
boundary  = Boundary(grid)
gauge     = Gauge(grid)

bulkevols      = BulkEvolvedPartition(grid)
bulkconstrains = BulkConstrainedPartition(grid)

# initial conditions
init_data!(bulkevols, boundary, gauge, systems, id)

# function to solve the nested system, given the initial data
solve_nested! = nested_solver(systems)

# solve nested system for the constrained variables
solve_nested!(bulkconstrains, bulkevols, boundary, gauge, evoleq)

# analyze data

i = 2

bulk = Bulk(bulkevols[i], bulkconstrains[i])

chart = atlas.charts[i]
uu    = chart.coords[1][:]
xx    = chart.coords[2][:]
yy    = chart.coords[3][:]

bulk.A[:,1,1]

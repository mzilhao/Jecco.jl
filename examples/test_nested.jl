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

# ibvp = AdS5_3_1.IDTest0(
#     b14_0      = 0.01,
#     b24_0      = 0.02,
#     phi0       = 0.0,
#     phi2_0     = 0.01,
#     a4_0       = -1.0,
#     fx2_0      = 0.02,
#     fy2_0      = 0.1,
#     potential  = ZeroPotential(),
# )

ibvp = BlackBrane()

# atlas of grid configuration and respective SystemPartition
atlas     = Atlas(grid)
systems   = SystemPartition(grid)

# evolved variables
bulkevols = BulkEvols(grid)
boundary  = Boundary(grid)
gauge     = Gauge(grid)

# and their initial conditions
init_data!(bulkevols, systems, ibvp)
init_data!(boundary, systems[1],   ibvp)
init_data!(gauge,    systems[end], ibvp)

# function to solve the nested system, given the initial data
solve_nested = nested_solver(systems, ibvp)

# initialize all bulk variables
bulks = Bulk.(bulkevols)

# solve nested system
solve_nested(bulks, boundary, gauge)

# analyze data

i = 2

chart = atlas.charts[i]
uu    = chart.coords[1][:]
xx    = chart.coords[2][:]
yy    = chart.coords[3][:]

bulk  = bulks[i]

bulk.A[:,1,1]

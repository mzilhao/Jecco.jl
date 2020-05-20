using Jecco
using Jecco.AdS5_3_1

import Base.Threads.@threads
import Base.Threads.@spawn
# using LinearAlgebra

par_grid = Grid3D(
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

# par_evol = ParamEvol(
#     ODE_method = "AB3",
#     # ODE_method = "RK4",
#     dt      = 0.008,
#     tmax    = 1.0,
# )


ibvp = AdS5_3_1.IDTest0(
    b14_0      = 0.01,
    b24_0      = 0.02,
    phi0       = 0.0,
    phi2_0     = 0.01,
    a4_0       = -1.0,
    fx2_0      = 0.02,
    fy2_0      = 0.1,
    potential  = ZeroPotential(),
)

systems = SystemPartition(par_grid)

bulkevols = BulkEvols(par_grid)

boundary  = Boundary(par_grid)
gauge     = Gauge(par_grid)

init_data!(bulkevols, systems, ibvp)
init_data!(boundary, systems[1],   ibvp)
init_data!(gauge,    systems[end], ibvp)


bulks = Bulk.(bulkevols)

solve_nested = nested_solver(systems, ibvp)
solve_nested(bulks, boundary, gauge)

i = 2

u = systems[i].ucoord[:]

bulk = bulks[i]
# BC   = BCs[i]
# dBC  = dBCs[i]
# nested = nesteds[i]

bulk.A[:,1,1]

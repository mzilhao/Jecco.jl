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
    u_outer_max      =  1.003,
    # u_outer_domains  =  1,
    # u_outer_nodes    =  64,
    u_outer_domains  =  3,
    u_outer_nodes    =  24,
    u_inner_nodes    =  12,
)

id   = BlackBrane()

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
gauge_t   = similar(gauge)

bulkevols      = BulkEvolvedPartition(grid)
bulkconstrains = BulkConstrainedPartition(grid)
bulkderivs     = BulkDerivPartition(grid)

# initial conditions
init_data!(bulkevols, boundary, gauge, systems, id)

# function to solve the nested system, given the initial data
nested = Nested(systems, bulkconstrains, bulkderivs)

# solve nested system for the constrained variables
nested(bulkevols, boundary, gauge, evoleq)


cache  = AdS5_3_1.HorizonCache(systems[end], 2)

AdS5_3_1.compute_xi_t!(gauge_t, bulkconstrains[end], bulkevols[end], bulkderivs[end],
                       gauge, cache, systems[end], evoleq.gaugecondition)

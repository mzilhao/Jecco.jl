
using Jecco
using Jecco.AdS5_3_1

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


systems   = SystemPartition(par_grid)

bulkevols = BulkEvols(par_grid)

boundary  = Boundary(par_grid)
gauge     = Gauge(par_grid)


# bulks = Bulk.(bulkevols)


ibvp = BlackBrane()

init_data!(bulkevols, systems, ibvp)
init_data!(boundary, systems[1],   ibvp)
init_data!(gauge,    systems[end], ibvp)


evolpartition0 = EvolPartition([boundary.a4, boundary.fx2, boundary.fy2, gauge.xi])

evolpartition1 = EvolPartition(boundary, gauge, bulkevols)

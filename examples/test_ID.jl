
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

# ibvp = BlackBrane()
ibvp = AdS5_3_1.IDTest0(
    b14_0 = 0.01,
    b24_0 = 0.02,
    g4_0  = 0.0,

    a4_0  = -1.0,
    fx2_0 = 0.001,
    fy2_0 = 0.002,

    xi_0  = 0.0,
)

init_data!(bulkevols, systems, ibvp)
init_data!(boundary, systems[1],   ibvp)
init_data!(gauge,    systems[end], ibvp)

evolpartition = EvolPartition(boundary, gauge, bulkevols)

boundary_  = AdS5_3_1.getboundary(evolpartition)
gauge_     = AdS5_3_1.getgauge(evolpartition)
bulkevol1_ = AdS5_3_1.getbulkevol(evolpartition, 1)
bulkevol2_ = AdS5_3_1.getbulkevol(evolpartition, 2)
bulkevol3_ = AdS5_3_1.getbulkevol(evolpartition, 3)

@show boundary_  === boundary
@show gauge_     === gauge
@show bulkevol1_ === bulkevols[1]
@show bulkevol2_ === bulkevols[2]
@show bulkevol3_ === bulkevols[3]

using Jecco
using Jecco.AdS5_3_1

grid = SpecCartGrid3D(
    x_min            = -5.0,
    x_max            =  5.0,
    x_nodes          =  32,
    y_min            = -5.0,
    y_max            =  5.0,
    y_nodes          =  48,
    u_outer_min      =  0.1,
    u_outer_max      =  1.005,
    u_outer_domains  =  5,
    u_outer_nodes    =  32,
    u_inner_nodes    =  12,
)

id = BlackBranePert(
    # B1_amp  = 0.005,
    # B1_amp = 0.0,
    # B1_amp = 0.0001,
    B1_nx   = 1,
    B1_ny   = 1,
    a4_ampx = 0.01,
    a4_ampy = 0.002,
    a4_kx   = 1,
    a4_ky   = 1,
    xmax    = grid.x_max,
    xmin    = grid.x_min,
    ymax    = grid.y_max,
    ymin    = grid.y_min,
    AH_pos  = 1.0,
)

evoleq = AffineNull(
    phi0           = 0.0,
    potential      = ZeroPotential(),
    gaugecondition = ConstantAH(u_AH = id.AH_pos),
)

ahf = AHF()

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

# AH_fd_order  = 4
AH_fd_order  = 2

horizoncache   = AdS5_3_1.HorizonCache(systems[end], AH_fd_order)

# initial conditions
init_data!(bulkconstrains, bulkevols, boundary, gauge, systems, evoleq, id)

# guess
# uAH = id.AH_pos
uAH = 0.9
# uAH = 1.9
# uAH = 2.0

sigma = similar(gauge.xi)
fill!(sigma, 1/uAH)

res = 0 * sigma

AdS5_3_1.find_AH!(sigma, bulkconstrains[end], bulkevols[end], bulkderivs[end], gauge,
                  horizoncache, systems[end], ahf)

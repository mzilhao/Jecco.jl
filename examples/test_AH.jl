using Jecco
using Jecco.AdS5_3_1

grid = SpecCartGrid3D(
    x_min            = -5.0,
    x_max            =  5.0,
    x_nodes          =  32,
    # x_nodes          =  12,
    y_min            = -5.0,
    y_max            =  5.0,
    # y_nodes          =  32,
    y_nodes          =  12,
    u_outer_min      =  0.1,
    # u_outer_max      =  2.005,
    u_outer_max      =  1.005,
    # u_outer_domains  =  3,
    u_outer_domains  =  5,
    u_outer_nodes    =  32,
    # u_outer_nodes    =  48,
    # u_inner_nodes    =  12,
    u_inner_nodes    =  24,
)

id = BlackBranePert(
    # B1_amp  = 0.005,
    # B1_amp = 0.0,
    # B1_amp = 0.0001,
    B1_nx  = 1,
    B1_ny  = 1,
    a4_amp = 0.00001,
    # a4_amp = 0.0,
    a4_k   = 1,
    xmax   = grid.x_max,
    xmin   = grid.x_min,
    ymax   = grid.y_max,
    ymin   = grid.y_min,
    # AH_pos = 2.0,
    AH_pos = 1.0,
)

evoleq = AffineNull(
    phi0           = 0.0,
    potential      = ZeroPotential(),
    gaugecondition = ConstantAH(u_AH = id.AH_pos),
)

ahf = AHF(itmax=10)
# ahf = AHF(itmax=0)
# ahf = AHF()

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
# uAH = 1.0
# uAH = 1.9
# uAH = 2.0

sigma = similar(gauge.xi)
fill!(sigma, uAH)

res = 0 * sigma


AdS5_3_1.compute_coeffs_AH!(sigma, gauge, horizoncache, systems[end])

AdS5_3_1.find_AH!(sigma, bulkconstrains[end], bulkevols[end], bulkderivs[end], gauge,
                  horizoncache, systems[end], ahf)

Du   = systems[end].Du
Duu  = systems[end].Duu
Dx   = systems[end].Dx
Dxx  = systems[end].Dxx

sigma_x  = Dx * sigma
sigma_xx = Dxx * sigma
S        = horizoncache.bulkhorizon.S_uAH
Sd       = horizoncache.bulkhorizon.Sd_uAH
Sp       = -sigma.^2 .* horizoncache.bulkhorizon.Du_S_uAH
Sdp      = -sigma.^2 .* horizoncache.bulkhorizon.Du_Sd_uAH
Spp      =  2 .* sigma.^3 .* horizoncache.bulkhorizon.Du_S_uAH +
    sigma.^4 .* horizoncache.bulkhorizon.Duu_S_uAH

# should be zero, but let's check
S_x     = Dx * horizoncache.bulkhorizon.S_uAH - horizoncache.bulkhorizon.Du_S_uAH .* sigma_x

tmp = Sp .* sigma_x .- S_x
bxx = tmp[1,:,:]

axx = -S[1,:,:]

tmp = 6 * S .* Sd .* Sp .+ 3 * S.^2 .* Sdp + 0.5 * Spp .* sigma_x.^2 .- Sp .* sigma_xx
cc  = tmp[1,:,:]

axx_ = reshape(horizoncache.axx, (grid.x_nodes,grid.y_nodes))
bxx_ = reshape(horizoncache.bx,  (grid.x_nodes,grid.y_nodes))
cc_  = reshape(horizoncache.cc,  (grid.x_nodes,grid.y_nodes))

diffa = abs.(axx-axx_)
diffb = abs.(bxx-bxx_)
diffc = abs.(cc - cc_)

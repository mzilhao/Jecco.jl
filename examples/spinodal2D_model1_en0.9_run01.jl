
using Jecco.AdS5_3_1

grid = SpecCartGrid3D(
    x_min            = -5.0,
    x_max            =  5.0,
    x_nodes          =  32,
    y_min            = -5.0,
    y_max            =  5.0,
    y_nodes          =  32,
    u_outer_min      =  0.1,
    u_outer_max      =  1.005,
    u_outer_domains  =  3,
    u_outer_nodes    =  32,
    u_inner_nodes    =  12,
    fd_order         =  6,
    sigma_diss       =  0.2,
)

potential = Phi8Potential(
    oophiM2 = -1.0,
    oophiQ  = 0.1,
)

id = BlackBranePert(
    energy_dens = 0.9,
    phi0        = 1.0,
    phi2        = 0.29819,
    oophiM2     = potential.oophiM2,
    a4_ampx     = -1.e-2,
    a4_ampy     = -1.e-2,
    a4_kx       = 1,
    a4_ky       = 1,
    AH_pos      = 1.0,
    xmax        = grid.x_max,
    xmin        = grid.x_min,
    ymin        = grid.y_min,
    ymax        = grid.y_max,
)

evoleq = AffineNull(
    phi0           = 1.0,
    potential      = potential,
    gaugecondition = ConstantAH(u_AH = 1.0),
)

diag = DiagAH(
    find_AH_every_t    = 1.0,
)

io = InOut(
    out_boundary_every_t        = 0.1,
    out_bulk_every_t            = 0.5,
    out_gauge_every_t           = 0.1,
    out_bulkconstrained_every_t = 0.5,
    checkpoint_every_walltime_hours = 1,
)

integration = Integration(
    tmax            = 200.0,
    ODE_method      = AdS5_3_1.VCABM3(),
    adaptive        = true,
    filter_poststep = true,
)

run_model(grid, id, evoleq, diag, integration, io)

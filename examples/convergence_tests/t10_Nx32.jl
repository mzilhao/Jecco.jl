
using Jecco.AdS5_3_1

grid = SpecCartGrid3D(
    x_min            = -5.0,
    x_max            =  5.0,
    x_nodes          =  32,
    y_min            = -5.0,
    y_max            =  5.0,
    y_nodes          =  32,
    u_outer_min      =  0.1,
    u_outer_max      =  1.2,
    u_outer_domains  =  3,
    u_outer_nodes    =  28,
    u_inner_nodes    =  12,
    fd_order         =  2,
    sigma_diss       =  0.01,
)

id = BlackBrane_xi1(
    a40    = -1.0,
    AH_pos = 1.0,
    xi_0   = 0.0,
    xi_nx  = 1,
    xi_Ax  = 0.1,
    xi_ny  = 1,
    xi_Ay  = 0.0,
    xmax = grid.x_max,
    xmin = grid.x_min,
    ymax = grid.y_max,
    ymin = grid.y_min,
)

evoleq = AffineNull(
    phi0           = 0.0,
    potential      = ZeroPotential(),
    gaugecondition = Advect_xi(xi_vx = 1.0, xi_vy = 0.0),
)

diag = DiagAH(
    find_AH_every_t    = 2.0,
)

io = InOut(
    out_boundary_every_t  = 0.1,
    out_bulk_every_t      = 0.1,
    out_gauge_every_t     = 0.1,
    out_bulkconstrained_every_t = 0.1,
    checkpoint_every_walltime_hours = 1,
)

integration = Integration(
    #dt              = 0.0025,
    tmax            = 10.,
    ODE_method      = AdS5_3_1.VCABM3(),
    filter_poststep = true,
)

@time run_model(grid, id, evoleq, diag, integration, io)

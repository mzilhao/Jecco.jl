
using Jecco.AdS5_3_1

grid = SpecCartGrid3D(
    x_min            = -10.0,
    x_max            =  10.0,
    x_nodes          =  128,
    y_min            = -5.0,
    y_max            =  5.0,
    y_nodes          =  6,
    u_outer_min      =  0.1,
    u_outer_max      =  1.01,
    u_outer_domains  =  1,
    u_outer_nodes    =  48,
    u_inner_nodes    =  12,
    fd_order         =  4,
    sigma_diss       =  0.2,
)

id = BlackBraneB1_a4_pert(
    energy_dens = 1.0,
    AH_pos = 1.0,
)

evoleq = AffineNull(
    phi0           = 0.0,
    potential      = ZeroPotential(),
    gaugecondition = ConstantAH(u_AH = 1.0),
)

io = InOut(
    out_boundary_every  = 200,
    out_bulk_every      = 1000,
    out_gauge_every     = 200,
    folder              = "",
)

integration = Integration(
    dt              = 0.001,
    tmax            = 2.0,
    ODE_method      = AdS5_3_1.AB4(),
    filter_poststep = true,
)

run_model(grid, id, evoleq, integration, io)
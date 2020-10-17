
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
    gaugecondition = ConstantAH(u_AH=0.99),
)

io = InOut(
    out_boundary_every  = 50,
    out_bulk_every      = 200,
    out_gauge_every     = 10,
    folder              = "./data",
    remove_existing     = true,
)

integration = Integration(
    # dt   = 0.001,
    dt   = 0.002,
    tmax = 0.01,
    # tmax = 3.0,
)

run_model(grid, id, evoleq, integration, io)

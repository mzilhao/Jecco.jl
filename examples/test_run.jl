
# using Jecco
using Jecco.AdS5_3_1

grid = SpecCartGrid3D(
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

id   = BlackBrane()

evoleq = AffineNull(
    phi0       = 0.0,
    potential  = ZeroPotential(),
)

io = InOut(
    out_boundary_every  = 1,
    out_bulk_every      = 2,
    out_gauge_every     = 1,
    folder              = "./data",
    overwrite           = true,
)

integration = Integration(
    dt   = 0.001,
    tmax = 0.01,
)

run(grid, id, evoleq, integration, io)

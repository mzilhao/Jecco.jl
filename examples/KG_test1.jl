using Jecco.KG_3_1

grid = SpecCartGrid3D(
    x_min            = -5.0,
    x_max            =  5.0,
    x_nodes          =  64,
    y_min            = -5.0,
    y_max            =  5.0,
    y_nodes          =  64,
    u_max            =  1.001,
    u_domains        =  1,
    u_nodes          =  48,
)

potential = SquarePotential()

evoleq = AffineNull(
    potential      = potential,
)

id = Sine2D(
    kx   = 1,
    ky   = 2,
    Lx   = grid.x_max - grid.x_min,
    Ly   = grid.y_max - grid.y_min,
)

io = InOut(
    out_bulk_every_t            = 0.04,
    checkpoint_every_walltime_hours = 1,
    remove_existing   = true,
)

integration = Integration(
    dt              = 0.002,
    ODE_method      = KG_3_1.RK4(),
    adaptive        = false,
    tmax            = 4.0,
)

run_model(grid, id, evoleq, integration, io)

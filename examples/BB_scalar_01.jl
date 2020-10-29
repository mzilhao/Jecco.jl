
using Jecco.AdS5_3_1

grid = SpecCartGrid3D(
    x_min            = -5.0,
    x_max            =  5.0,
    x_nodes          =  8,
    y_min            = -5.0,
    y_max            =  5.0,
    y_nodes          =  8,
    u_outer_min      =  0.1,
    u_outer_max      =  1.003,
    u_outer_domains  =  3,
    u_outer_nodes    =  24,
    u_inner_nodes    =  12,
    fd_order         =  4,
    sigma_diss       =  0.2,
)

potential = Phi8Potential(
    oophiM2 = 1/10^2,
    oophiQ  = 0.0,
)

id = PhiGaussian_u(
    energy_dens   = 1.0,
    AH_pos        = 1.0,
    phi0          = 1.0,
    phi2          = 0.1,
    amp           = 0.005,
    u0            = 0.5,
    sigma         = 0.15,
    oophiM2       = potential.oophiM2,
)

evoleq = AffineNull(
    phi0           = 1.0,
    potential      = potential,
    gaugecondition = ConstantAH(u_AH = 1.0),
)

io = InOut(
    out_boundary_every  = 10,
    out_bulk_every      = 200,
    out_gauge_every     = 10,

    checkpoint_every_walltime_hours = 1,
)

integration = Integration(
    dt              = 0.002,
    tmax            = 5.0,
    ODE_method      = AdS5_3_1.AB3(),
    filter_poststep = true,
)

run_model(grid, id, evoleq, integration, io)

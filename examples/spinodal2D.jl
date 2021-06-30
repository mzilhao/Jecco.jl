
using Jecco.AdS5_3_1

grid = SpecCartGrid3D(
    x_min            = -50.,
    x_max            =  50.,
    x_nodes          =  125,
    y_min            = -50.,
    y_max            =  50.,
    y_nodes          =  125,
    u_outer_min      =  0.1,
    u_outer_max      =  1.005,
    u_outer_domains  =  1,
    u_outer_nodes    =  48,
    u_inner_nodes    =  12,
    fd_order         =  4,
    sigma_diss       =  0.2,
)


phiM2 = -(0.85)^2
phiQ  = 10.

potential = AdS5_3_1.Phi8Potential(
    #alpha   = -0.01,
    #beta    = 8,
    #gamma   = 0.1,
    oophiM2 = 1/phiM2,
    oophiQ  = 1/phiQ,
)

id = BlackBranePert(
    energy_dens = .209,
    phi0        = 1.0,
    phi2        = 0.9,
    oophiM2     = potential.oophiM2,
    #phi5        = 1.1,
    #a4_ampx     = -0.05,
    #a4_kx       = 1,
    #a4_ampy     = -0.05,
    #a4_ky       = 1,
    xi0         = 0.02,
    AH_pos      = 0.9,
    xmax        = grid.x_max,
    xmin        = grid.x_min,
    ymin        = grid.y_min,
    ymax        = grid.y_max,
)

evoleq = AffineNull(
    phi0           = id.phi0,
    potential      = potential,
    gaugecondition = ConstantAH(u_AH = id.AH_pos),
)

diag = DiagAH(
    find_AH_every_t    = 0.1,
)

outdir = "/home/mikel/Documents/Jecco.jl/data/test/"

io = InOut(
    out_boundary_every_t        = 1.,
    out_bulk_every_t            = 1.,
    out_gauge_every_t           = 1.,
    #out_bulkconstrained_every_t = 5.0,
    checkpoint_every_walltime_hours = 1,
    out_dir                     = outdir,
    recover                     = :yes,
    recover_dir                 = "/home/mikel/Documents/Jecco.jl/data/new_data/",
    checkpoint_dir              = outdir,
    remove_existing             = true,
)

integration = Integration(
    #dt              = 0.0002,
    tmax            = 60.,
    ODE_method      = AdS5_3_1.VCABM3(),
    #ODE_method      = AdS5_3_1.AB3(),
    adaptive        = true,
    filter_poststep = true,
)

run_model(grid, id, evoleq, diag, integration, io)

convert_to_mathematica(io.out_dir)

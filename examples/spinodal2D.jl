using Jecco.AdS5_3_1

grid = SpecCartGrid3D(
    x_min            = -31.8,
    x_max            =  31.8,
    x_nodes          =  50,
    y_min            = -31.8,
    y_max            =  31.8,
    y_nodes          =  50,
    u_outer_min      =  0.1,
    u_outer_max      =  1.005,
    u_outer_domains  =  3,
    u_outer_nodes    =  24,
    u_inner_nodes    =  12,
    fd_order         =  4,
    sigma_diss       =  0.2,
)


phiM2 = -(1.0)^2
phiQ  = 10.

potential = AdS5_3_1.Phi8Potential(
    #alpha   = -0.01,
    #beta    = 8,
    #gamma   = 0.1,
    oophiM2 = 1/phiM2,
    oophiQ  = 1/phiQ,
)

id = AdS5_3_1.BlackBraneGaussPert(
    energy_dens = 0.95,
    phi0        = 1.0,
    phi2        = 0.3,
    oophiM2     = potential.oophiM2,
    #phi5        = 1.1,
    #a4_ampx     = -0.005,
    #a4_kx       = 1,
    #a4_ampy     = -0.005,
    #a4_ky       = 1,
    #xi0         = 0.02,
    sigma       = 1e-3,
    AH_pos      = 1.0,
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
    find_AH_every_t    = 1.,
)

outdir = "/home/mikel/Documents/Jecco.jl/data/test/"

io = InOut(
    out_boundary_every_t        = 5.,
    out_bulk_every_t            = 10.,
    out_gauge_every_t           = 5.,
    #out_bulkconstrained_every_t = 5.0,
    checkpoint_every_walltime_hours = 1,
    out_dir                     = outdir,
    recover                     = :no,
    recover_dir                 = "/home/mikel/Documents/Jecco.jl/data/new_data/",
    checkpoint_dir              = outdir,
    remove_existing             = true,
)

integration = Integration(
    #dt              = 0.0002,
    tmax            = 200.,
    ODE_method      = AdS5_3_1.VCABM3(),
    #ODE_method      = AdS5_3_1.AB3(),
    adaptive        = true,
    filter_poststep = true,
)

run_model(grid, id, evoleq, diag, integration, io)

convert_to_mathematica(io.out_dir)

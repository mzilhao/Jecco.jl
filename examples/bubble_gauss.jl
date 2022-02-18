using Jecco.AdS5_3_1

grid = SpecCartGrid3D(
    x_min            = -40.,
    x_max            =  40.,
    x_nodes          =  320,
    y_min            = -.5,
    y_max            =  .5,
    y_nodes          =  4,
    u_outer_min      =  0.1,
    u_outer_max      =  1.005,
    u_outer_domains  =  3,
    u_outer_nodes    =  24,
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

id = AdS5_3_1.Gaussian_a4(
    energy_dens = 1.318,
    phi0        = 1.0,
    phi2        = 0.3,
    oophiM2     = potential.oophiM2,
    xi0         = 0.0,
    amp         = -4.e-1,
    sigma       = 10.,
    AH_pos      = 0.95,
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

outdir = "/Users/apple/Documents/Jecco.jl/data/test_2_new_resol/"

io = InOut(
    out_boundary_every_t        = .5,
    out_bulk_every_t            = .5,
    out_gauge_every_t           = .5,
    #out_bulkconstrained_every_t = 5.0,
    checkpoint_every_walltime_hours = 1,
    out_dir                     = outdir,
    recover                     = :no,
    recover_dir                 = outdir,
    checkpoint_dir              = outdir,
    remove_existing             = false,
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

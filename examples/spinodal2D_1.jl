
using Jecco.AdS5_3_1

grid = SpecCartGrid3D(
    x_min            = -10.,
    x_max            =  10.,
    x_nodes          =  80,
    y_min            = -10.0,
    y_max            =  10.0,
    y_nodes          =  80,
    u_outer_min      =  0.1,
    u_outer_max      =  1.005,
    u_outer_domains  =  1,
    u_outer_nodes    =  48,
    u_inner_nodes    =  12,
    fd_order         =  4,
    sigma_diss       =  0.2,
)

potential = Phi8Potential(
    oophiM2 = -1.38408,
    oophiQ  = 0.1,
)

id = BlackBranePert(
    energy_dens = 1.0,
    phi0        = 1.0,
    phi2        = 0.3,
    oophiM2     = potential.oophiM2,
    a4_ampx     = -0.05,
    a4_ampy     = -0.05,
    a4_kx       = 1,
    a4_ky       = 1,
   # xi_0        = 0.1,
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
    find_AH_every_t    = 1.0,
)

outdir = "/home/mikel/Documents/Jecco.jl/data/test_3/"

io = InOut(
    out_boundary_every_t        = 1,
    out_bulk_every_t            = 1,
    out_gauge_every_t           = 1,
    out_bulkconstrained_every_t = 1,
    checkpoint_every_walltime_hours = 1,
    out_dir                     = outdir,
    recover                     = :yes,
    recover_dir                 = "/home/mikel/Documents/Jecco.jl/data/bubbles/phiM_0.85_phiQ_10/phase_separated/e_1.8_L_20_AH_0.95/",
    checkpoint_dir              = outdir,
    remove_existing             = true,
)

integration = Integration(
    tmax            = 10.0,
    ODE_method      = AdS5_3_1.VCABM3(),
    adaptive        = true,
    filter_poststep = true,
)

run_model(grid, id, evoleq, diag, integration, io)
# convert_to_mathematica(io.out_dir)

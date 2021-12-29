
using Jecco.AdS5_3_1

grid = SpecCartGrid3D(
    x_min            = -0.5,
    x_max            =  0.5,
    x_nodes          =  4,
    y_min            = -100.0,
    y_max            =  100.0,
    y_nodes          =  200,
    u_outer_min      =  0.1,
    u_outer_max      =  1.005,
    u_outer_domains  =  3,
    u_outer_nodes    =  28,
    u_inner_nodes    =  12,
    fd_order         =  4,
    #fd_order         =  6,
    sigma_diss       =  0.2,
)

potential = Phi8Potential(
    oophiM2 = -1.0,
    oophiQ  = 0.1,
)

id = BlackBranePert(
    energy_dens = 1.3,
    phi0        = 1.0,
    phi2        = 0.28,
    oophiM2     = -1.0,
    a4_ampx     = 0.0,
    a4_ampy     = 0.002,
    a4_kx       = 1,
    a4_ky       = 1,
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
    find_AH_every_t    = 1.0,
)

outdir = "/gpfs/projects/ub48/ub48946/Yago/E_1.3_a_0.002_k_0.0314_x_200_3_28_new_y/"

io = InOut(
    out_boundary_every_t        = 0.2,
    out_bulk_every_t            = 2.0,
    out_gauge_every_t           = 0.2,
    #out_bulkconstrained_every_t = 2.0,
    checkpoint_every_walltime_hours = 1,
    out_dir                     = outdir,
    recover                     = :no,
    #recover_dir                 = "/gpfs/projects/ub48/ub48139/jecco_runs/E_1.3_a_0.002_k_0.0314_x_200_3_28_new_y",
    checkpoint_dir              = outdir,
    remove_existing             = true,
)

integration = Integration(
    tmax            = 300.0,
    #ODE_method      = AdS5_3_1.VCABM3(),
    dt              = 0.001,
    ODE_method      = AdS5_3_1.AB3(),
    adaptive        = false,
    #adaptive        = true,
    filter_poststep = true,
)

run_model(grid, id, evoleq, diag, integration, io)

convert_to_mathematica(io.out_dir)

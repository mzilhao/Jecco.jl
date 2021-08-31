
using Jecco.AdS5_3_1

grid = SpecCartGrid3D(
    x_min            = -3.5,
    x_max            =  3.5,
    x_nodes          =  42,
    y_min            = -0.5,
    y_max            =  0.5,
    y_nodes          =  4,
    u_outer_min      =  0.1,
    u_outer_max      =  1.005,
    u_outer_domains  =  1,
    u_outer_nodes    =  64,
    u_inner_nodes    =  12,
    fd_order         =  4,
    sigma_diss       =  0.2,
)

potential = Phi8Potential(
    oophiM2 = -1.0,
    oophiQ  = 0.1,
)

id = BlackBranePert(
    energy_dens = 0.85,
    phi0        = 1.0,
    phi2        = 0.29819,
    oophiM2     = -1.0,
    a4_ampx     = 1.e-2,
    a4_ampy     = 0.0,
    a4_kx       = 1,
    a4_ky       = 0,
    AH_pos      = 1.0,
    xmax        = grid.x_max,
    xmin        = grid.x_min,
    ymin        = grid.y_min,
    ymax        = grid.y_max,
)

evoleq = AffineNull(
    phi0           = 1.0,
    potential      = potential,
    gaugecondition = ConstantAH(u_AH = 1.0),
)

diag = DiagAH(
    find_AH_every_t    = 1.0,
)

io = InOut(
    out_boundary_every_t        = 10.0,
    out_bulk_every_t            = 10.0,
    out_gauge_every_t           = 10.0,
    out_bulkconstrained_every_t = 10.0,
    checkpoint_every_walltime_hours = 1,
    out_dir                     ="/gpfs/projects/ub48/ub48946/jecco_test/spinodal1D_e_0.9_L_30_42/",
    recover                     = :no,
   # recover_dir                = "/home/ub48/ub48946/Jecco.jl/initial_data/",
    checkpoint_dir              = "/gpfs/projects/ub48/ub48946/jecco_test/spinodal1D_e_0.9_L_30_42/",
    remove_existing             = true, 
)

integration = Integration(
    tmax            = 30.0,
    ODE_method      = AdS5_3_1.VCABM3(),
    adaptive        = true,
    filter_poststep = true,
)

run_model(grid, id, evoleq, diag, integration, io)

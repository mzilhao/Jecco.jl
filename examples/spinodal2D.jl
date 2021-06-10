
using Jecco.AdS5_3_1

grid = SpecCartGrid3D(
    x_min            = -10.,
    x_max            =  10.,
    x_nodes          =  80,
    y_min            = -10.,
    y_max            =  10.,
    y_nodes          =  80,
    u_outer_min      =  0.1,
    u_outer_max      =  1.005,
    u_outer_domains  =  1,
    u_outer_nodes    =  48,
    u_inner_nodes    =  12,
    fd_order         =  4,
    sigma_diss       =  0.2,
)

potential = AdS5_3_1.Phi8Potential(
    #alpha   = -0.01,
    #beta    = 8,
    #gamma   = 0.1,
    oophiM2 = -1.156203029251937,
    oophiQ  = 0.05,
)

id = BlackBranePert(
    energy_dens = 0.63,
    phi0        = 1.0,
    phi2        = 0.3,
    oophiM2     = potential.oophiM2,
    #phi5        = 1.1,
    a4_ampx     = -0.05,
    a4_kx       = 1,
    a4_ampy     = -0.05,
    a4_ky       = 1,
    xi0         = 0.0,
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

outdir = "/home/mikel/Documents/Jecco.jl/data/spinodal2D_e_0.63_phiM_0.93i_phiQ_20_L_20_N_80/"

io = InOut(
    out_boundary_every_t        = 5.,
    out_bulk_every_t            = 5.,
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

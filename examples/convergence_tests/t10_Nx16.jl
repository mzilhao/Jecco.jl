
using Jecco.AdS5_3_1

# the numerical setup: grid, finite difference operators accuracy order, KO dissipation etc.
grid = SpecCartGrid3D(
    x_min            = -5.0,
    x_max            =  5.0,
    x_nodes          =  16,
    y_min            = -5.0,
    y_max            =  5.0,
    y_nodes          =  16,
    u_outer_min      =  0.1,
    u_outer_max      =  1.2,
    u_outer_domains  =  3,
    u_outer_nodes    =  28,
    u_inner_nodes    =  12,
    fd_order         =  2,
    sigma_diss       =  0.01,
)

# initial data
id = BlackBrane_xi1(
    a40    = -1.0,
    AH_pos = 1.0,
    xi_0   = 0.0,
    xi_nx  = 1,
    xi_Ax  = 0.1,
    xi_ny  = 1,
    xi_Ay  = 0.0,
    xmax = grid.x_max,
    xmin = grid.x_min,
    ymax = grid.y_max,
    ymin = grid.y_min,
)

# gauge specification
evoleq = AffineNull(
    phi0           = 0.0,
    potential      = ZeroPotential(),
    gaugecondition = Advect_xi(xi_vx = 1.0, xi_vy = 0.0),
)

# diagnostics, i.e. apparent horizon finder
# here the AH finder runs every t=2 (in code units)
diag = DiagAH(
    find_AH_every_t    = 2.0,
)

# input and output
# here t is in code units, i.e. every_t=0.1 the relevant output is written
io = InOut(
    out_boundary_every_t  = 0.1,
    out_bulk_every_t      = 0.1,
    out_gauge_every_t     = 0.1,
    out_bulkconstrained_every_t = 0.1,
    checkpoint_every_walltime_hours = 1,
)

# integration methos
# here we choose adaptive step, so the explicit specification of dt is commented out
# if you wish to specify a fixed timestep, just include dt and give the value
# in such case, even if you let the ODE_methof be an adaptive step one,
# it will use the fixed step, so effectivel it is a fixed step one
integration = Integration(
    #dt              = 0.0025,
    tmax            = 10.,
    ODE_method      = AdS5_3_1.VCABM3(),
    filter_poststep = true,
)

# run the mode (solve the AdS5 equations) for the above specifications
@time run_model(grid, id, evoleq, diag, integration, io)

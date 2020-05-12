
using Jecco
using Jecco.AdS5_3_1

par_grid = Grid3D(
    x_min            = -5.0,
    x_max            =  5.0,
    x_nodes          =  128,
    y_min            = -5.0,
    y_max            =  5.0,
    y_nodes          =  128,
    u_outer_min      =  0.1,
    u_outer_max      =  1.0,
    # u_outer_domains  =  1,
    # u_outer_nodes    =  64,
    u_outer_domains  =  2,
    u_outer_nodes    =  16,
    u_inner_nodes    =  12,
)


potential   = ZeroPotential()
phi0        = 0.0

base  = BaseVars(potential, phi0)

systems = Jecco.AdS5_3_1.Systems(par_grid)

sys = systems[1]

f  = Jecco.AdS5_3_1.init(sys, BlackBrane())

fs = Jecco.AdS5_3_1.init(systems, BlackBrane())


using Jecco
using Jecco.KG_3_1

p = Param(
    A0x         = 1.0,
    A0y         = 1.0,

    tmax        = 1.0,
    out_every_t = 1,

    xmin        = -5.0,
    xmax        =  5.0,
    xnodes      =  128,
    ymin        = -5.0,
    ymax        =  5.0,
    ynodes      =  128,
    umin        =  0.0,
    umax        =  1.0,
    unodes      =  32,

    dtfac       = 0.5,

    folder      = "./data",
)


initial_data = Jecco.KG_3_1.sine2D

sys = System(p)

phi = initial_data(sys, p)

bulk = BulkVars(phi)

a4 = -Jecco.KG_3_1.ones2D(sys)

boundary = BoundaryVars(a4)

Jecco.KG_3_1.Vf(phi)  = -1.0 + 0.5 * phi*phi
Jecco.KG_3_1.Vfp(phi) = phi


# solve_nested_g1! = Jecco.KG_3_1.nested_g1(sys)
# solve_nested_g1!(bulk, boundary)

rhs! = Jecco.KG_3_1.setup_rhs(phi, sys)
dphidt = similar(phi)


rhs!(dphidt, phi, sys, 0.0)

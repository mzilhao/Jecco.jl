
using Jecco
using Jecco.KG_3_1

p = Param(
    A0x         = 1.0,
    A0y         = 1.0,

    tmax        = 20.0,
    out_every   = 1,

    xmin        = -1.0,
    xmax        =  1.0,
    xnodes      =    6,
    ymin        = -1.0,
    ymax        =  1.0,
    ynodes      =    6,
    umin        =  0.0,
    umax        =  1.0,
    unodes      =  32,

    # dtfac       = 0.5,

    # dt          = 0.04, # for RK4
    # dt          = 0.08, # for RK4
    dt          = 0.1, # for RK4

    folder      = "./data",
)

Jecco.KG_3_1.Vf(phi)  = -1.0 + 0.5 * phi*phi
Jecco.KG_3_1.Vfp(phi) = phi

Jecco.KG_3_1.ibvp(p)


using Jecco
using Jecco.KG_3_1

p = Param(
    A0x         = 1.0,
    A0y         = 1.0,

    tmax        = 2.0,
    out_every   = 1,

    xmin        = -5.0,
    xmax        =  5.0,
    xnodes      =  128,
    ymin        = -5.0,
    ymax        =  5.0,
    ynodes      =  128,
    umin        =  0.0,
    umax        =  1.0,
    unodes      =  32,

    # dtfac       = 0.5,

    dt          = 0.06, # for RK4
    # dt          = 0.01,   # for AB3

    folder      = "./data",
)

Jecco.KG_3_1.Vf(phi)  = -1.0 + 0.5 * phi*phi
Jecco.KG_3_1.Vfp(phi) = phi

Jecco.KG_3_1.ibvp(p)

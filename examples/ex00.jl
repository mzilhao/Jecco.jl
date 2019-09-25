
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

phif = initial_data(sys, p)

# FIXME
phifd = copy(phif)
Af    = -ones(size(phif))

dphif = similar(phif)

bulk = BulkVars(phif, phifd, Af)

Jecco.KG_3_1.rhs!(dphif, bulk, sys)

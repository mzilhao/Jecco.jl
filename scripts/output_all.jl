using Jecco, Jecco.AdS5_3_1
import Jecco.AdS5_3_1: geta4, getfx2, getfy2, getxi, getB1, getB2, getG, getphi

grid = SpecCartGrid3D(
    x_min            = -5.0,
    x_max            =  5.0,
    x_nodes          =  32,
    y_min            = -5.0,
    y_max            =  5.0,
    y_nodes          =  32,
    u_outer_min      =  0.1,
    u_outer_max      =  1.005,
    u_outer_domains  =  3,
    u_outer_nodes    =  32,
    u_inner_nodes    =  12,
    fd_order         =  6,
    sigma_diss       =  0.2,
)

potential = Phi8Potential(
    oophiM2 = -1.0,
    oophiQ  = 0.1,
)

evoleq = AffineNull(
    phi0           = 1.0,
    potential      = potential,
    gaugecondition = ConstantAH(u_AH = 1.0),
)

dirname = "/home/mzilhao/dev/julia/Jecco/examples/spinodal2D_model1_en0.9_run01/"

it = 0

ts_bulk  = OpenPMDTimeSeries(dirname, prefix="bulk_")
ts_bdry  = OpenPMDTimeSeries(dirname, prefix="boundary_")
ts_gauge = OpenPMDTimeSeries(dirname, prefix="gauge_")


# atlas of grid configuration and respective SystemPartition
atlas     = Atlas(grid)
systems   = SystemPartition(grid)
Nsys      = length(systems)

# allocate variables
boundary       = Boundary(grid)
boundary_t     = Boundary(grid)
gauge          = Gauge(grid)
gauge_t        = Gauge(grid)
bulkevols      = BulkEvolvedPartition(grid)
bulkevols_t    = BulkEvolvedPartition(grid)
bulkconstrains = BulkConstrainedPartition(grid)
bulkderivs     = BulkDerivPartition(grid)
cache          = AdS5_3_1.HorizonCache(systems[end], evoleq.gaugecondition.fd_order)

charts = atlas.charts

# for the boundary/xi grid
empty   = Cartesian{1}("u", systems[1].ucoord[1], systems[1].ucoord[1], 1)
chart2D = Chart(empty, systems[1].xcoord, systems[1].ycoord)


# restore data
AdS5_3_1.restore!(bulkevols, ts_bulk, it)
AdS5_3_1.restore!(boundary, ts_bdry, it)
AdS5_3_1.restore!(gauge, ts_gauge, it)

t0 = ts_bulk.current_t

tinfo  = Jecco.TimeInfo(it, t0, 0.0, 0.0)


# nested system function
nested = AdS5_3_1.Nested(systems, bulkconstrains, bulkderivs)

AdS5_3_1.compute_boundary_t!(boundary_t, bulkevols[1], boundary, gauge, systems[1], evoleq)

# solve nested system for the constrained variables
nested(bulkevols, boundary, gauge, evoleq)

AdS5_3_1.compute_xi_t!(gauge_t, bulkconstrains[Nsys], bulkevols[Nsys],
                       bulkderivs[Nsys], gauge, cache, systems[Nsys],
                       evoleq.gaugecondition)

for aa in 1:Nsys
    sys           = systems[aa]
    bulkevol_t    = bulkevols_t[aa]
    bulkevol      = bulkevols[aa]
    bulkconstrain = bulkconstrains[aa]

    AdS5_3_1.compute_bulkevolved_t!(bulkevol_t, bulkconstrain, gauge_t, bulkevol,
                                    boundary, gauge, sys, evoleq)
end
AdS5_3_1.sync_bulkevolved!(bulkevols_t, bulkconstrains, gauge_t, systems, evoleq)


# output structures
out_bdry_t  = Jecco.Output(dirname, "boundary_t_", tinfo)
out_gauge_t = Jecco.Output(dirname, "gauge_t_", tinfo)
out_bulk_t  = Jecco.Output(dirname, "bulk_t_", tinfo)

# output fields
boundary_fields = (
    Jecco.Field("a4_t",   geta4(boundary_t),  chart2D),
    Jecco.Field("fx2_t",  getfx2(boundary_t), chart2D),
    Jecco.Field("fy2_t",  getfy2(boundary_t), chart2D),
)

gauge_fields = Jecco.Field("xi_t", getxi(gauge_t), chart2D)

bulkevols_fields = ntuple(i -> (
    Jecco.Field("B1_t c=$i",  getB1(bulkevols_t[i]),  charts[i]),
    Jecco.Field("B2_t c=$i",  getB2(bulkevols_t[i]),  charts[i]),
    Jecco.Field("G_t c=$i",   getG(bulkevols_t[i]),   charts[i]),
    Jecco.Field("phi_t c=$i", getphi(bulkevols_t[i]), charts[i])
), Nsys)


out_bdry_t(boundary_fields)
out_gauge_t(gauge_fields)
out_bulk_t.(bulkevols_fields)

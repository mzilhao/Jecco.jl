
using OrdinaryDiffEq

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


systems   = SystemPartition(par_grid)

bulkevols = BulkEvols(par_grid)

boundary  = Boundary(par_grid)
gauge     = Gauge(par_grid)


# bulks = Bulk.(bulkevols)


ibvp = BlackBrane()

init_data!(bulkevols, systems, ibvp)
init_data!(boundary, systems[1],   ibvp)
init_data!(gauge,    systems[end], ibvp)

evol = EvolPartition(boundary, gauge, bulkevols)


abstract type AbstractEvolEq end

struct EvolTest0 <: AbstractEvolEq end


function get_f_t!(ff_t, ff, systems, evoleq::EvolTest0)
    Nsys = length(systems)

    boundary  = AdS5_3_1.getboundary(ff)
    gauge     = AdS5_3_1.getgauge(ff)
    bulkevols = AdS5_3_1.getbulkevols(ff)

    boundary_t  = AdS5_3_1.getboundary(ff_t)
    gauge_t     = AdS5_3_1.getgauge(ff_t)
    bulkevols_t = AdS5_3_1.getbulkevols(ff_t)

    fill!(boundary_t.a4, 0)
    fill!(boundary_t.fx2, 0)
    fill!(boundary_t.fy2, 0)
    fill!(gauge_t.xi, 0)

    for aa in 1:Nsys
        sys = systems[aa]
        bulkevol   = bulkevols[aa]
        bulkevol_t = bulkevols_t[aa]

        fill!(bulkevol_t.B1,  0)
        fill!(bulkevol_t.B2,  0)
        fill!(bulkevol_t.G,   0)
        fill!(bulkevol_t.phi, 0)
    end
    nothing
end

function rhs!(df, f, (systems, evoleq), t)
    get_f_t!(df, f, systems, evoleq)
end


evoleq = EvolTest0()

dt0 = 0.001
tspan = (0.0, 0.05)

prob  = ODEProblem(rhs!, evol, tspan, (systems, evoleq))

integrator = init(prob, RK4(), save_everystep=false, dt=dt0, adaptive=false)

for (f,t) in tuples(integrator)
    B1  = AdS5_3_1.getB1(f, 1)
    a4  = AdS5_3_1.geta4(f)
    @show t, B1[1], a4[1]
end

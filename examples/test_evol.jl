
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

ibvp   = BlackBrane()

evoleq = EvolTest0()


systems   = SystemPartition(par_grid)

bulkevols = BulkEvols(par_grid)

boundary  = Boundary(par_grid)
gauge     = Gauge(par_grid)


init_data!(bulkevols, systems, ibvp)
init_data!(boundary, systems[1],   ibvp)
init_data!(gauge,    systems[end], ibvp)

evol = EvolPartition(boundary, gauge, bulkevols)


function rhs!(df, f, (systems, evoleq), t)
    get_f_t!(df, f, systems, evoleq)
end


dt0 = 0.001
tspan = (0.0, 0.05)

prob  = ODEProblem(rhs!, evol, tspan, (systems, evoleq))

integrator = init(prob, RK4(), save_everystep=false, dt=dt0, adaptive=false)

for (f,t) in tuples(integrator)
    B1  = AdS5_3_1.getB1(f, 1)
    a4  = AdS5_3_1.geta4(f)
    @show t, B1[1], a4[1]
end

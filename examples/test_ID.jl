
using Jecco
using Jecco.AdS5_3_1

using OrdinaryDiffEq
# using RecursiveArrayTools

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


systems = System(par_grid)

ibvp  = BlackBrane()

# evols = init_data(systems, ibvp)

sys = systems[1]

Nu, Nx, Ny = size(sys)

bulkevols = BulkEvol(par_grid)
boundary  = Boundary(par_grid)
gauge     = Gauge(par_grid)



# init_data!(bulkevol, sys, ibvp)
# init_data!(boundary, systems[1],   ibvp)
# init_data!(gauge,    systems[end], ibvp)



# abstract type AbstractEvolEq end

# struct EvolTest0 <: AbstractEvolEq end


# function get_evol_t!(evol_t::EvolVars, evol::EvolVars, sys::System, ::EvolTest0)
#     evol_t.B1  .= 0
#     evol_t.B2  .= 0
#     evol_t.G   .= 0
#     evol_t.phi .= 0
#     evol_t.a4  .= 0
#     evol_t.fx2 .= 0
#     evol_t.fy2 .= 0
#     evol_t.xi  .= 0
#     nothing
# end

# function rhs!(df, f, (sys, evoleq), t)
#     get_evol_t!(df, f, sys, evoleq)
#     nothing
# end

# evol = evols[1]
# sys  = systems[1]

# evoleq = EvolTest0()

# evol_dt = similar(evol)

# dt0 = 0.001
# tspan = (0.0, 0.01)

# prob  = ODEProblem(rhs!, evol, tspan, (sys, evoleq))

# integrator = init(prob, RK4(), save_everystep=false, dt=dt0, adaptive=false)

# for (f,t) in tuples(integrator)
#     B1  = Jecco.AdS5_3_1.getB1(f)
#     a4  = Jecco.AdS5_3_1.geta4(f)
#     @show t, B1[1], a4[1]
# end

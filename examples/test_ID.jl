
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


systems = Jecco.AdS5_3_1.Systems(par_grid)

ibvp  = BlackBrane()

evols = init_data(systems, ibvp)


abstract type AbstractEvolEq end

struct EvolTest0 <: AbstractEvolEq end


function get_evol_t!(evol_t::EvolVars, evol::EvolVars, sys, ::EvolTest0)
    evol_t.B1  .= 0
    evol_t.B2  .= 0
    evol_t.G   .= 0
    evol_t.phi .= 0
    evol_t.a4  .= 0
    evol_t.fx2 .= 0
    evol_t.fy2 .= 0
    evol_t.xi  .= 0
    nothing
end

function rhs!(df, f, (sys, evoleq), t)
    get_evol_t!(df, f, sys, evoleq)
    nothing
end

evol = evols[1]
sys  = systems[1]

evoleq = EvolTest0()

evol_dt = similar(evol)

tspan = (0.0, 0.1)

prob  = ODEProblem(rhs!, evol, tspan, (sys, evoleq))

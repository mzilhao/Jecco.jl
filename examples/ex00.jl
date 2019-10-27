
using Vivi
using Jecco
using Jecco.KG_3_1

using DifferentialEquations

p = Param(
    A0x         = 1.0,
    A0y         = 1.0,

    tmax        = 2.0,
#    tmax        = 8.0,
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


initial_data = Jecco.KG_3_1.sine2D

ucoord  = Vivi.SpectralCoord("u", p.umin, p.umax, p.unodes)
xcoord  = Vivi.CartCoord("x", p.xmin, p.xmax, p.xnodes, endpoint=false)
ycoord  = Vivi.CartCoord("y", p.ymin, p.ymax, p.ynodes, endpoint=false)


sys = System(ucoord, xcoord, ycoord)

phi0 = initial_data(sys, p)

bulk = BulkVars(phi0)


Jecco.KG_3_1.Vf(phi)  = -1.0 + 0.5 * phi*phi
Jecco.KG_3_1.Vfp(phi) = phi

rhs! = Jecco.KG_3_1.setup_rhs(phi0, sys)
# timestep = Jecco.KG_3_1.timestep
# timestep = p.dt

# dphidt = similar(phi)
# rhs!(dphidt, phi, sys, 0.0)

# dt0 = timestep(sys, phi0)
dt0 = p.dt

tspan = (0.0, p.tmax)

prob  = ODEProblem(rhs!, phi0, tspan, sys)
# http://docs.juliadiffeq.org/latest/basics/integrator.html
integrator = init(prob, RK4(), save_everystep=false, dt=dt0, adaptive=false)
# integrator = init(prob, AB3(), save_everystep=false, dt=dt0, adaptive=false)


out    = Vivi.Output(p.folder, p.prefix, p.out_every)

it = 0
# write initial data
Jecco.out_info(it, 0, phi0, "phi", 1, 200)
Vivi.output(out, Dict("phi" => (phi0, sys.coords)), it, 0, 0)

for (u,t) in tuples(integrator)
    # it += 1
    global it += 1

    dt = integrator.dt

    Jecco.out_info(it, t, u, "phi", 1, 200)
    Vivi.output(out, Dict("phi" => (u, sys.coords)), it, t, dt)
end

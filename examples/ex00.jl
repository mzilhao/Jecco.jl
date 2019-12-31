
using Vivi
using Jecco
using Jecco.KG_3_1

using DifferentialEquations


par_base = ParamBase(
    which_potential = "square",
)

par_grid = ParamGrid(
    xmin        = -5.0,
    xmax        =  5.0,
    xnodes      =  128,
    ymin        = -5.0,
    ymax        =  5.0,
    ynodes      =  128,
    umin        =  0.0,
    umax        =  1.0,
    unodes      =  64,
)

par_id = ParamID(
    ID_type      = "sine2D",
    A0x          = 1.0,
    A0y          = 1.0,
    Lx           = par_grid.xmax - par_grid.xmin,
    Ly           = par_grid.ymax - par_grid.ymin,
)

par_evol = ParamEvol(
    dt      = 0.02, # for RK4
    tmax    = 4.0,
)

par_io = ParamIO(
    out_every   = 4,
    folder      = "./data",
)

# define potential
Jecco.KG_3_1.setup(par_base)

initial_data = Jecco.KG_3_1.sine2D

ucoord  = Vivi.SpectralCoord("u", par_grid.umin, par_grid.umax, par_grid.unodes)
xcoord  = Vivi.CartCoord("x", par_grid.xmin, par_grid.xmax, par_grid.xnodes, endpoint=false)
ycoord  = Vivi.CartCoord("y", par_grid.ymin, par_grid.ymax, par_grid.ynodes, endpoint=false)


sys = System(ucoord, xcoord, ycoord)

phi0 = initial_data(sys, par_id)


Jecco.KG_3_1.Vf(phi)  = -1.0 + 0.5 * phi*phi
Jecco.KG_3_1.Vfp(phi) = phi

rhs! = Jecco.KG_3_1.setup_rhs(phi0, sys)

# timestep = Jecco.KG_3_1.timestep

# dt0 = timestep(sys, phi0)
dt0 = par_evol.dt

tspan = (0.0, par_evol.tmax)

prob  = ODEProblem(rhs!, phi0, tspan, sys)
# http://docs.juliadiffeq.org/latest/basics/integrator.html
integrator = init(prob, RK4(), save_everystep=false, dt=dt0, adaptive=false)
# integrator = init(prob, AB3(), save_everystep=false, dt=dt0, adaptive=false)


tinfo  = Vivi.TimeInfo()
out    = Vivi.Output(par_io.folder, par_io.prefix, par_io.out_every, tinfo)

# write initial data
Jecco.out_info(tinfo.it, tinfo.t, phi0, "phi", 1, 200)
Vivi.output(out, "phi", phi0, sys.coords)

for (u,t) in tuples(integrator)
    tinfo.it += 1
    tinfo.dt  = integrator.dt
    tinfo.t   = t

    Jecco.out_info(tinfo.it, tinfo.t, u, "phi", 1, 200)
    Vivi.output(out, "phi", u, sys.coords)
end

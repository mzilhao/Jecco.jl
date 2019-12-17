
using Vivi
using Jecco
using Jecco.KG_3_1

using DifferentialEquations

function unpack_dom(ucoord)
    Nsys = length(ucoord)
    Nus  = [ucoord[i].nodes for i in 1:Nsys]

    Nu_lims = zeros(typeof(Nsys), Nsys + 1)
    for i in 1:Nsys
        Nu_lims[i+1] = Nu_lims[i] + Nus[i]
    end

    function (f)
        [f[Nu_lims[i]+1:Nu_lims[i+1],:,:] for i in 1:Nsys]
    end
end

function write_out(out, fieldnames, coordss)
    Nsys   = length(fieldnames)
    ucoord = [coordss[i][1] for i in 1:Nsys]
    unpack = unpack_dom(ucoord)

    function (u)
        phis = unpack(u)
        Vivi.output(out, fieldnames, phis, coordss)
        nothing
    end
end


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
    udomains    =  4,
    unodes      =  16,
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
    # tmax    = 4.0,
    tmax    = 0.4,
)

par_io = ParamIO(
    out_every   = 4,
    folder      = "./data00_multi",
)


initial_data = Jecco.KG_3_1.sine2D

Nsys = par_grid.udomains

delta_udom = (par_grid.umax - par_grid.umin) / Nsys

ucoord = [Vivi.SpectralCoord("u", par_grid.umin + (i-1)*delta_udom, par_grid.umin + i*delta_udom,
                             par_grid.unodes) for i in 1:Nsys]
xcoord  = Vivi.CartCoord("x", par_grid.xmin, par_grid.xmax, par_grid.xnodes, endpoint=false)
ycoord  = Vivi.CartCoord("y", par_grid.ymin, par_grid.ymax, par_grid.ynodes, endpoint=false)

systems = [System(ucoord[i], xcoord, ycoord) for i in 1:Nsys]

phi0s = initial_data(systems, par_id)
ID    = vcat(phi0s...)




unpack = unpack_dom(ucoord)

rhs! = Jecco.KG_3_1.setup_rhs(phi0s, systems, unpack)


# timestep = Jecco.KG_3_1.timestep
# timestep = p.dt


# dt0 = timestep(sys, phi0)
dt0 = par_evol.dt

tspan = (0.0, par_evol.tmax)


prob  = ODEProblem(rhs!, ID, tspan, systems)
# http://docs.juliadiffeq.org/latest/basics/integrator.html
integrator = init(prob, RK4(), save_everystep=false, dt=dt0, adaptive=false)
# integrator = init(prob, AB3(), save_everystep=false, dt=dt0, adaptive=false)

tinfo  = Vivi.TimeInfo()

# write initial data
Jecco.out_info(tinfo.it, tinfo.t, ID, "phi", 1, 200)

fieldnames = ["phi c=$i" for i in 1:Nsys]
fields     = phi0s
coordss    = [systems[i].coords for i in 1:Nsys]

out    = Vivi.Output(par_io.folder, par_io.prefix, par_io.out_every, tinfo)
output = write_out(out, fieldnames, coordss)
output(ID)

for (u,t) in tuples(integrator)
    tinfo.it += 1
    tinfo.dt  = integrator.dt
    tinfo.t   = t

    Jecco.out_info(tinfo.it, tinfo.t, u, "phi", 1, 200)
    output(u)
end

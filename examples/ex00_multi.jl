
using Vivi
using Jecco
using Jecco.KG_3_1

using DifferentialEquations


function write_out(out, fieldnames, coordss)
    Nsys = length(fieldnames)
    function (u)
        phis = [u.x[i] for i in 1:Nsys]
        Vivi.output(out, fieldnames, phis, coordss)
        nothing
    end
end


p = Param(
    A0x         = 1.0,
    A0y         = 1.0,

    tmax        = 4.0,
#    tmax        = 8.0,
    out_every   = 4,

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

    # dtfac       = 0.5,

    dt          = 0.02, # for RK4

    folder      = "./data00_multi",
)


initial_data = Jecco.KG_3_1.sine2D

Nsys = p.udomains

delta_udom = (p.umax - p.umin) / Nsys

ucoord = [Vivi.SpectralCoord("u", p.umin + (i-1)*delta_udom, p.umin + i*delta_udom,
                             p.unodes) for i in 1:Nsys]
xcoord  = Vivi.CartCoord("x", p.xmin, p.xmax, p.xnodes, endpoint=false)
ycoord  = Vivi.CartCoord("y", p.ymin, p.ymax, p.ynodes, endpoint=false)

systems = [System(ucoord[i], xcoord, ycoord) for i in 1:Nsys]

phi0s = initial_data(systems, p)

bulks = BulkVars(phi0s)

phi0s_slice  = [phi0[1,:,:] for phi0 in phi0s]

boundaries = BulkVars(phi0s_slice)


Jecco.KG_3_1.Vf(phi)  = -1.0 + 0.5 * phi*phi
Jecco.KG_3_1.Vfp(phi) = phi


rhs! = Jecco.KG_3_1.setup_rhs(phi0s, systems)

# timestep = Jecco.KG_3_1.timestep
# timestep = p.dt

# dphidt = similar(phi)
# rhs!(dphidt, phi, sys, 0.0)

# dt0 = timestep(sys, phi0)
dt0 = p.dt

tspan = (0.0, p.tmax)

ID = ArrayPartition(phi0s...)

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
vars_dict = Dict("phi c=1" => (phi0s[1], systems[1].coords) )

out    = Vivi.Output(p.folder, p.prefix, p.out_every, tinfo)
output = write_out(out, fieldnames, coordss)
output(ID)

for (u,t) in tuples(integrator)
    tinfo.it += 1
    tinfo.dt  = integrator.dt
    tinfo.t   = t

    # phis = [u.x[i] for i in 1:Nsys]

    # Jecco.out_info(tinfo.it, tinfo.t, phis[1], "phi c=1", 1, 200)
    Jecco.out_info(tinfo.it, tinfo.t, u, "phi", 1, 200)
    output(u)
end

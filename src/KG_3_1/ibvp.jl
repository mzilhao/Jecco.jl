
using DifferentialEquations
using Vivi

function ibvp(p::Param)

    # TODO: make parameter?
    initial_data = Jecco.KG_3_1.sine2D
    # initial_data = Jecco.KG_3_1.uniform2D

    ucoord  = Vivi.SpectralCoord("u", p.umin, p.umax, p.unodes)
    xcoord  = Vivi.CartCoord("x", p.xmin, p.xmax, p.xnodes, endpoint=false)
    ycoord  = Vivi.CartCoord("y", p.ymin, p.ymax, p.ynodes, endpoint=false)

    sys = System(ucoord, xcoord, ycoord)

    phi0 = initial_data(sys, p)
    bulk = BulkVars(phi0)

    a4 = -Jecco.KG_3_1.ones2D(sys)
    boundary = BoundaryVars(a4)

    rhs! = Jecco.KG_3_1.setup_rhs(phi0, sys)

    # timestep = Jecco.KG_3_1.timestep
    # dt0 = timestep(sys, phi0)

    dt0 = p.dt

    tspan = (0.0, p.tmax)

    prob  = ODEProblem(rhs!, phi0, tspan, sys)
    # http://docs.juliadiffeq.org/latest/basics/integrator.html
    integrator = init(prob, RK4(), save_everystep=false, dt=dt0, adaptive=false)
    # integrator = init(prob, AB3(), save_everystep=false, dt=dt0, adaptive=false)

    tinfo  = Vivi.TimeInfo()
    out    = Vivi.Output(p.folder, p.prefix, p.out_every, tinfo)

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

    nothing
end

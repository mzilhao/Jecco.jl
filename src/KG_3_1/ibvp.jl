
using DifferentialEquations
using Vivi

function ibvp(p::Param)

    # TODO: make parameter?
    initial_data = Jecco.KG_3_1.sine2D
    # initial_data = Jecco.KG_3_1.uniform2D

    sys = System(p)

    phi0 = initial_data(sys, p)
    bulk = BulkVars(phi0)

    a4 = -Jecco.KG_3_1.ones2D(sys)
    boundary = BoundaryVars(a4)

    rhs! = Jecco.KG_3_1.setup_rhs(phi0, sys)
    timestep = Jecco.KG_3_1.timestep

    dt0 = timestep(sys, phi0)

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
        it += 1

        dt = integrator.dt

        Jecco.out_info(it, t, u, "phi", 1, 200)
        Vivi.output(out, Dict("phi" => (u, sys.coords)), it, t, dt)
    end

    nothing
end

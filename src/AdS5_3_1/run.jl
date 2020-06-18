
function run_model(grid::SpecCartGrid3D, id::InitialData, evoleq::EvolutionEquations,
                   integration::Integration, io::InOut)
    Jecco.startup()

    # atlas of grid configuration and respective SystemPartition
    atlas     = Atlas(grid)
    systems   = SystemPartition(grid)
    Nsys      = length(systems)

    # allocate variables
    boundary       = Boundary(grid)
    gauge          = Gauge(grid)
    bulkevols      = BulkEvolvedPartition(grid)
    bulkconstrains = BulkConstrainedPartition(grid)
    bulkderivs     = BulkDerivPartition(grid)
    horizoncache   = HorizonCache(systems[end], evoleq.gaugecondition.fd_order)

    # initial conditions
    init_data!(bulkevols, boundary, gauge, systems, id)

    # full state vector
    evolvars  = EvolVars(boundary, gauge, bulkevols)

    # function that updates the state vector
    rhs! = setup_rhs(bulkconstrains, bulkderivs, horizoncache, systems, integration)

    dt0  = integration.dt
    tmax = integration.tmax

    tspan = (0.0, tmax)
    alg   = integration.ODE_method

    prob  = ODEProblem(rhs!, evolvars, tspan, evoleq)
    # http://docs.juliadiffeq.org/latest/basics/integrator.html
    integrator = init(prob, alg, save_everystep=false, dt=dt0, adaptive=false)

    tinfo  = Jecco.TimeInfo()

    # for the boundary/xi grid
    empty   = Cartesian{1}("u", 0.0, 0.0, 1)
    chart2D = Chart(empty, systems[1].xcoord, systems[1].ycoord)

    # prepare functions to write data
    output_evol = output_writer(evolvars, chart2D, atlas.charts, tinfo, io)

    if io.out_bulkconstrained_every > 0
        output_constrained = output_writer(bulkconstrains, atlas.charts, tinfo, io)
    else
        output_constrained = x -> nothing
    end

    # prepare checkpointing function
    if io.checkpoint_every_walltime_hours > 0
        checkpoint = checkpoint_writer(evolvars, chart2D, atlas.charts, tinfo, io)
    else
        checkpoint = x -> nothing
    end

    # write initial data
    output_evol(evolvars)
    output_constrained(bulkconstrains)

    # for stdout info
    Jecco.out_info(tinfo.it, tinfo.t, 0.0, gauge.xi, "ξ", 1, 200)

    tstart = time()
    t0     = tinfo.t
    # start integration
    for (u,t) in tuples(integrator)
        tinfo.it     += 1
        tinfo.dt      = integrator.dt
        tinfo.t       = t
        tinfo.runtime = time() - tstart

        # write data
        output_evol(u)
        output_constrained(bulkconstrains)

        # checkpoint
        checkpoint(u)

        gauge = getgauge(u)

        telapsed = tinfo.runtime / 3600
        deltat   = t - t0
        Jecco.out_info(tinfo.it, tinfo.t, deltat/telapsed, gauge.xi, "ξ", 1, 200)
    end

    println("-------------------------------------------------------------")
    println("Done.")
end

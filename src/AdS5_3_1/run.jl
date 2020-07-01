
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

    do_id = true
    it0   = 0
    t0    = 0.0

    # are we recovering from a previous checkpoint?
    if io.recover == :yes
        it0, t0 = recover(bulkevols, boundary, gauge, io.recover_dir)
        do_id   = false
    elseif io.recover == :auto
        try
            it0, t0 = recover(bulkevols, boundary, gauge, io.recover_dir)
            do_id   = false
        catch e
            if isa(e, ErrorException) && e.msg == "No files found."
                println("INFO: No checkpoint data found. Will run initial data.")
            else
                throw(e)
            end
        end
    elseif io.recover == :no
        # do nothing; run the initial data functions instead.
    else
        error("Unknown option for io.recover")
    end

    # if we are not recovering from a checkpoint, run initial conditions
    if do_id
        init_data!(bulkevols, boundary, gauge, systems, id)
    end

    # full state vector
    evolvars  = EvolVars(boundary, gauge, bulkevols)

    # function that updates the state vector
    rhs! = setup_rhs(bulkconstrains, bulkderivs, horizoncache, systems, integration)

    dt0  = integration.dt
    tmax = integration.tmax

    # decide in the evolution loop when to terminate the run, so set here an
    # impossibly large value for tstop
    tspan = (0.0, 1.e20)
    alg   = integration.ODE_method

    prob  = ODEProblem(rhs!, evolvars, tspan, evoleq)
    # http://docs.juliadiffeq.org/latest/basics/integrator.html
    integrator = init(prob, alg, save_everystep=false, dt=dt0, adaptive=false)

    tinfo  = Jecco.TimeInfo(it0, t0, 0.0, 0.0)

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
    last_checkpoint_walltime = 0.0

    # remove termination trigger file, if it exists
    if io.termination_from_file
        finish_him = abspath(io.out_dir, io.termination_file)
        rm(finish_him, force=true)
    end

    # write initial data
    output_evol(evolvars)
    output_constrained(bulkconstrains)

    # for stdout info
    Jecco.out_info(tinfo.it, tinfo.t, 0.0, gauge.xi, "ξ", 1, 1)

    tstart = time()
    # start integration
    for (u,t) in tuples(integrator)
        tinfo.it     += 1
        tinfo.dt      = integrator.dt
        tinfo.t       = t0 + t
        tinfo.runtime = time() - tstart
        telapsed      = tinfo.runtime / 3600

        # write info to stdout
        gauge  = getgauge(u)
        deltat = t
        Jecco.out_info(tinfo.it, tinfo.t, deltat/telapsed, gauge.xi, "ξ", 1, 200)

        # write data
        output_evol(u)
        output_constrained(bulkconstrains)

        # checkpoint
        if telapsed >= last_checkpoint_walltime + io.checkpoint_every_walltime_hours
            last_checkpoint_walltime = telapsed
            checkpoint(u)
        end

        # terminate run?
        if t >= tmax || telapsed >= io.max_walltime
            checkpoint(u)
            terminate!(integrator)
        end
        if io.termination_from_file && tinfo.it % io.check_file_every == 0
            if isfile(finish_him)
                println("INFO: Found termination file.")
                println("INFO: Triggering termination...")
                checkpoint(u)
                terminate!(integrator)
            end
        end

    end

    println("-------------------------------------------------------------")
    println("Done.")
end


function estimate_dtmax(chart::Chart)
    ucoord, xcoord, ycoord = chart.coords
    dx = Jecco.delta(xcoord)
    dy = Jecco.delta(ycoord)
    # spacing in u is not uniform, so let's compute the average spacing
    du_avg = (ucoord[end] - ucoord[1]) / (ucoord.nodes - 1)
    0.8 * min(dx, dy, du_avg)
end
function estimate_dtmax(atlas::Atlas)
    dtmaxs = estimate_dtmax.(atlas)
    minimum(dtmaxs)
end

function run_model(grid::SpecCartGrid3D, id::InitialData, evoleq::EvolutionEquations,
                   integration::Integration, io::InOut)
    Jecco.startup(io.out_dir; remove_existing=io.remove_existing)

    # atlas of grid configuration and respective SystemPartition
    atlas     = Atlas(grid)
    systems   = SystemPartition(grid)
    Nsys      = length(systems)

    # allocate variables
    boundary       = Boundary(grid)
    bulkevols      = BulkEvolvedPartition(grid)
    bulkconstrains = BulkConstrainedPartition(grid)
    bulkderivs     = BulkDerivPartition(grid)

    do_id = true
    it0   = 0
    t0    = 0.0

    # are we recovering from a previous checkpoint?
    if io.recover == :yes
        it0, t0 = recover(bulkevols, boundary, systems, id, io.recover_dir)
        do_id   = false
    elseif io.recover == :auto
        try
            it0, t0 = recover(bulkevols, boundary, systems, id, io.recover_dir)
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
        id(bulkevols, boundary, systems)
    end

    # full state vector
    evolvars  = EvolVars(bulkevols)

    # function that updates the state vector
    rhs! = setup_rhs(bulkconstrains, boundary, bulkderivs, systems, integration)

    #=
    limit the default integrator dtmax and qmax values. see:
      https://diffeq.sciml.ai/latest/extras/timestepping/
      https://diffeq.sciml.ai/latest/basics/common_solver_opts/
    =#
    dtmax = estimate_dtmax(atlas)
    qmax  = 1.2

    if isa(integration.dt, Number)
        dt0   = integration.dt
    elseif integration.dt == :auto
        dt0   = 0.5 * dtmax
    else
        error("Unknown dt value")
    end
    tmax  = integration.tmax

    # decide in the evolution loop when to terminate the run, so set here an
    # impossibly large value for tstop
    tspan = (0.0, 1.e20)
    alg   = integration.ODE_method

    prob  = ODEProblem(rhs!, evolvars, tspan, evoleq)
    # https://diffeq.sciml.ai/stable/basics/integrator/
    integrator = init(prob, alg, save_everystep=false, dt=dt0, dtmax=dtmax, qmax=qmax,
                      adaptive=integration.adaptive, reltol=integration.reltol,
                      calck=false)

    tinfo  = Jecco.TimeInfo(it0, t0, 0.0, 0.0)

    # for the boundary grid
    empty   = Cartesian{1}("u", systems[1].ucoord[1], systems[1].ucoord[1], 1)
    chart2D = Chart(empty, systems[1].xcoord, systems[1].ycoord)

    # prepare functions to write data
    output_evol = output_writer(evolvars, chart2D, atlas, tinfo, io,
                                evoleq.potential)

    # prepare checkpointing function
    if io.checkpoint_every_walltime_hours > 0
        checkpoint = checkpoint_writer(evolvars, chart2D, atlas, tinfo, io)
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
    if do_id
        output_evol(evolvars)
    end

    # for stdout info
    Jecco.out_info(tinfo.it, tinfo.t, 0.0, getphi(bulkevols[1]), "ϕ c=1", 1, 1)

    tstart = time()
    # start integration
    for (u,t) in tuples(integrator)
        tinfo.it     += 1
        tinfo.dt      = integrator.dt
        tinfo.t       = t0 + t
        tinfo.runtime = time() - tstart
        telapsed      = tinfo.runtime / 3600

        # write info to stdout
        bulkevols = getbulkevolvedpartition(u)
        deltat    = t
        Jecco.out_info(tinfo.it, tinfo.t, deltat/telapsed, getphi(bulkevols[1]), "ϕ c=1", 1, 200)

        # write data
        output_evol(u)

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

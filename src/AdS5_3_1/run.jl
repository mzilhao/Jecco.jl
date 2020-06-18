
Base.@kwdef struct InOut
    out_boundary_every :: Int
    out_bulk_every     :: Int
    out_gauge_every    :: Int
    folder             :: String  = "./data"
    # be very careful with this option! it will remove the whole folder contents
    # if set to true! use only for fast debugging runs
    remove_existing    :: Bool    = false
end

function output_writer(u::EvolVars, chart2D::Chart, charts, tinfo::Jecco.TimeInfo,
                       io::InOut)
    Nsys = length(charts)

    # output structures
    out_bdry  = Jecco.Output(io.folder, "boundary_", tinfo;
                             remove_existing=io.remove_existing)
    out_gauge = Jecco.Output(io.folder, "gauge_", tinfo;
                             remove_existing=io.remove_existing)
    out_bulk  = Jecco.Output(io.folder, "bulk_", tinfo;
                             remove_existing=io.remove_existing)

    boundary  = getboundary(u)
    gauge     = getgauge(u)
    bulkevols = getbulkevolvedpartition(u)

    #=
    for analysis, it's useful to have the boundary behaviour of the bulk fields,
    so let's output those as well in the boundary file. note that the slicing
    operation [1,:,:] converts 3D arrays into 2D arrays; since we want the same
    shape as in the remaining boundary fields -- 3D arrays of size (1,Nx,Ny) --
    we must reshape them first.
    =#
    b14  = reshape(bulkevols[1].B1[1,:,:],  size(chart2D))
    b24  = reshape(bulkevols[1].B2[1,:,:],  size(chart2D))
    g4   = reshape(bulkevols[1].G[1,:,:],   size(chart2D))
    phi3 = reshape(bulkevols[1].phi[1,:,:], size(chart2D))

    #= output fields

    note that the field "phi3" is just the leading boundary behaviour of the
    "phi" field in the inner grid. to compute the usual "phi2" quantity one
    needs to do (at the post-processing stage)

      phi2 = phi0^3 phi3 - phi0 xi^2

    at a later stage we may wish to output this quantity directly
    =#
    boundary_fields = (
        Jecco.Field("a4",   boundary.a4,  chart2D),
        Jecco.Field("fx2",  boundary.fx2, chart2D),
        Jecco.Field("fy2",  boundary.fy2, chart2D),
        Jecco.Field("b14",  b14,          chart2D),
        Jecco.Field("b24",  b24,          chart2D),
        Jecco.Field("g4",   g4,           chart2D),
        Jecco.Field("phi3", phi3,         chart2D),
    )
    gauge_fields = Jecco.Field("xi", gauge.xi, chart2D)
    bulkevols_fields = ntuple(i -> (
        Jecco.Field("B1 c=$i",  bulkevols[i].B1,  charts[i]),
        Jecco.Field("B2 c=$i",  bulkevols[i].B2,  charts[i]),
        Jecco.Field("G c=$i",   bulkevols[i].G,   charts[i]),
        Jecco.Field("phi c=$i", bulkevols[i].phi, charts[i])
    ), Nsys)

    function (u::EvolVars)
        boundary  = getboundary(u)
        gauge     = getgauge(u)
        bulkevols = getbulkevolvedpartition(u)

        it = tinfo.it

        if it % io.out_boundary_every == 0
            boundary_fields[1].data = boundary.a4
            boundary_fields[2].data = boundary.fx2
            boundary_fields[3].data = boundary.fy2
            @views copyto!(boundary_fields[4].data, bulkevols[1].B1[1,:,:])
            @views copyto!(boundary_fields[5].data, bulkevols[1].B2[1,:,:])
            @views copyto!(boundary_fields[6].data, bulkevols[1].G[1,:,:])
            @views copyto!(boundary_fields[7].data, bulkevols[1].phi[1,:,:])
            # write data
            out_bdry(boundary_fields)
        end

        if it % io.out_gauge_every == 0
            gauge_fields.data = gauge.xi
            # write data
            out_gauge(gauge_fields)
        end

        if it % io.out_bulk_every == 0
            @inbounds for i in 1:Nsys
                bulkevols_fields[i][1].data = bulkevols[i].B1
                bulkevols_fields[i][2].data = bulkevols[i].B2
                bulkevols_fields[i][3].data = bulkevols[i].G
                bulkevols_fields[i][4].data = bulkevols[i].phi
            end
            # write data
            out_bulk.(bulkevols_fields)
        end

        nothing
    end
end

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

    # prepare function to write data
    output = output_writer(evolvars, chart2D, atlas.charts, tinfo, io)

    # write initial data
    output(evolvars)

    # for stdout info
    Jecco.out_info(tinfo.it, tinfo.t, 0.0, gauge.xi, "ξ", 1, 200)

    tstart = time()
    t0     = tinfo.t
    # start integration
    for (u,t) in tuples(integrator)
        tinfo.it += 1
        tinfo.dt  = integrator.dt
        tinfo.t   = t

        # write data
        output(u)

        gauge = getgauge(u)

        telapsed = (time() - tstart) / 3600
        deltat   = t - t0
        Jecco.out_info(tinfo.it, tinfo.t, deltat/telapsed, gauge.xi, "ξ", 1, 200)
    end

    println("-------------------------------------------------------------")
    println("Done.")
end

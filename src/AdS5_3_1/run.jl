
Base.@kwdef struct Integration
    dt          :: Float64
    tmax        :: Float64
    # ODE_method  :: String   = "RK4"
end

Base.@kwdef struct InOut
    out_boundary_every :: Int
    out_bulk_every     :: Int
    out_gauge_every    :: Int
    folder             :: String  = "./data"
    overwrite          :: Bool    = false
end

function run(grid::SpecCartGrid3D, id::InitialData, evoleq::EvolutionEquations,
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

    # initial conditions
    init_data!(bulkevols, boundary, gauge, systems, id)

    # full state vector
    evolvars  = EvolVars(boundary, gauge, bulkevols)

    # function that updates the state vector
    rhs! = setup_rhs(bulkconstrains, systems)

    dt0  = integration.dt
    tmax = integration.tmax

    tspan = (0.0, tmax)
    # TODO: make parameter
    alg   = AB3()

    prob  = ODEProblem(rhs!, evolvars, tspan, evoleq)
    # http://docs.juliadiffeq.org/latest/basics/integrator.html
    integrator = init(prob, alg, save_everystep=false, dt=dt0, adaptive=false)

    tinfo  = Jecco.TimeInfo()

    # output structures
    out_bdry  = Jecco.Output(io.folder, "boundary_", io.out_boundary_every, tinfo;
                             overwrite=io.overwrite)
    out_bulk  = Jecco.Output(io.folder, "bulk_", io.out_bulk_every, tinfo;
                             overwrite=io.overwrite)
    out_gauge = Jecco.Output(io.folder, "gauge_", io.out_gauge_every, tinfo;
                             overwrite=io.overwrite)

    # output fields
    charts  = atlas.charts

    # for the boundary/xi grid
    empty   = Cartesian{1}("u", 0.0, 0.0, 1)
    chart2D = Chart(empty, systems[1].xcoord, systems[1].ycoord)

    boundary_fields = (
        Jecco.Field("a4",  boundary.a4,  chart2D),
        Jecco.Field("fx2", boundary.fx2, chart2D),
        Jecco.Field("fy2", boundary.fy2, chart2D),
    )
    gauge_fields = Jecco.Field("xi", gauge.xi, chart2D)

    tmp = [ [Jecco.Field("B1 c=$i",  bulkevols[i].B1,  charts[i]),
             Jecco.Field("B2 c=$i",  bulkevols[i].B2,  charts[i]),
             Jecco.Field("G c=$i",   bulkevols[i].G,   charts[i]),
             Jecco.Field("phi c=$i", bulkevols[i].phi, charts[i])] for i in 1:Nsys ]
    bulkevols_fields = [tmp...;]

    # write initial data
    out_bdry(boundary_fields)
    out_gauge(gauge_fields)
    out_bulk(bulkevols_fields)

    # for stdout info
    Aend = bulkconstrains[end].A

    Jecco.out_info(tinfo.it, tinfo.t, 0.0, Aend, "A c=$Nsys", 1, 200)

    tstart = time()
    t0     = tinfo.t
    # start integration
    for (u,t) in tuples(integrator)
        tinfo.it += 1
        tinfo.dt  = integrator.dt
        tinfo.t   = t

        boundary  = getboundary(u)
        gauge     = getgauge(u)
        bulkevols = getbulkevolvedpartition(u)

        # output
        # FIXME: this bulkevols is not pointing to the one inside the
        # bulkevols_fields define above, so the output information won't be
        # updated... need to see a better way...
        out_bdry(boundary_fields)
        out_gauge(gauge_fields)
        out_bulk(bulkevols_fields)

        telapsed = (time() - tstart) / 3600
        deltat   = t - t0
        Jecco.out_info(tinfo.it, tinfo.t, deltat/telapsed, Aend, "A c=$Nsys", 1, 200)
    end

    println("-------------------------------------------------------------")
    println("Done.")
end

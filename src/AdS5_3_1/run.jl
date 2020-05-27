
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

    # SystemPartition for this grid configuration
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


    # write initial data

    # phis = unpack(ID)
    # fieldnames = ["phi c=$i" for i in 1:Nsys]
    # charts     = [Chart([sys.ucoord, sys.xcoord, sys.ycoord]) for sys in systems]
    # fields     = [Jecco.Field(fieldnames[i], phis[i], charts[i]) for i in 1:Nsys]


    # output = write_out(out, fields)
    # output(ID)

    Aend = bulkconstrains[end].A

    Jecco.out_info(tinfo.it, tinfo.t, 0.0, Aend, "A c=$Nsys", 1, 200)

    tstart = time()
    t0     = tinfo.t
    for (u,t) in tuples(integrator)
        tinfo.it += 1
        tinfo.dt  = integrator.dt
        tinfo.t   = t

        # output(u)

        telapsed = (time() - tstart) / 3600
        deltat   = t - t0
        Jecco.out_info(tinfo.it, tinfo.t, deltat/telapsed, Aend, "A c=$Nsys", 1, 200)
    end

    println("-------------------------------------------------------------")
    println("Done.")
end


function output_writer(u::EvolVars, chart2D::Chart, charts, tinfo::Jecco.TimeInfo,
                       io::InOut)
    Nsys = length(charts)

    # output structures
    out_bdry  = Jecco.Output(io.out_dir, "boundary_", tinfo;
                             remove_existing=io.remove_existing)
    out_gauge = Jecco.Output(io.out_dir, "gauge_", tinfo;
                             remove_existing=io.remove_existing)
    out_bulk  = Jecco.Output(io.out_dir, "bulk_", tinfo;
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


function output_writer(bulkconstrains::BulkPartition{Nsys,BulkConstrained{T}}, charts,
                       tinfo::Jecco.TimeInfo, io::InOut) where {Nsys,T}
    @assert Nsys == length(charts)

    # output structure
    out  = Jecco.Output(io.out_dir, "constrained_", tinfo;
                        remove_existing=io.remove_existing)

    # output fields
    fields = ntuple(i -> (
        Jecco.Field("S c=$i",    bulkconstrains[i].S,    charts[i]),
        Jecco.Field("Fx c=$i",   bulkconstrains[i].Fx,   charts[i]),
        Jecco.Field("Fy c=$i",   bulkconstrains[i].Fy,   charts[i]),
        Jecco.Field("B1d c=$i",  bulkconstrains[i].B1d,  charts[i]),
        Jecco.Field("B2d c=$i",  bulkconstrains[i].B2d,  charts[i]),
        Jecco.Field("Gd c=$i",   bulkconstrains[i].Gd,   charts[i]),
        Jecco.Field("phid c=$i", bulkconstrains[i].phid, charts[i]),
        Jecco.Field("Sd c=$i",   bulkconstrains[i].Sd,   charts[i]),
        Jecco.Field("A c=$i",    bulkconstrains[i].A,    charts[i])
    ), Nsys)

    function (bulkconstrains)
        it = tinfo.it

        if it % io.out_bulkconstrained_every == 0
            @inbounds for i in 1:Nsys
                fields[i][1].data = bulkconstrains[i].S
                fields[i][2].data = bulkconstrains[i].Fx
                fields[i][3].data = bulkconstrains[i].Fy
                fields[i][4].data = bulkconstrains[i].B1d
                fields[i][5].data = bulkconstrains[i].B2d
                fields[i][6].data = bulkconstrains[i].Gd
                fields[i][7].data = bulkconstrains[i].phid
                fields[i][8].data = bulkconstrains[i].Sd
                fields[i][9].data = bulkconstrains[i].A
            end
            # write data
            out.(fields)
        end

        nothing
    end
end


function checkpoint_writer(u::EvolVars, chart2D::Chart, charts, tinfo::Jecco.TimeInfo,
                           io::InOut)
    Nsys = length(charts)

    # output structure
    out  = Jecco.Output(io.checkpoint_dir, "checkpoint_it", tinfo;
                        remove_existing=false)

    boundary  = getboundary(u)
    gauge     = getgauge(u)
    bulkevols = getbulkevolvedpartition(u)

    # output fields
    boundary_fields = (
        Jecco.Field("a4",   boundary.a4,  chart2D),
        Jecco.Field("fx2",  boundary.fx2, chart2D),
        Jecco.Field("fy2",  boundary.fy2, chart2D),
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

        boundary_fields[1].data = boundary.a4
        boundary_fields[2].data = boundary.fx2
        boundary_fields[3].data = boundary.fy2

        gauge_fields.data = gauge.xi

        @inbounds for i in 1:Nsys
            bulkevols_fields[i][1].data = bulkevols[i].B1
            bulkevols_fields[i][2].data = bulkevols[i].B2
            bulkevols_fields[i][3].data = bulkevols[i].G
            bulkevols_fields[i][4].data = bulkevols[i].phi
        end

        # write data
        out(boundary_fields)
        out(gauge_fields)
        out.(bulkevols_fields)

        nothing
    end
end


# restore all bulk evolved fields
function restore!(bulkevols::BulkPartition{Nsys}, ts::OpenPMDTimeSeries,
                  it::Int) where {Nsys}
    for i in 1:Nsys
        B1, chart = get_field(ts, it=it, field="B1 c=$i")
        @assert size(B1) == size(bulkevols[i].B1)
        copyto!(bulkevols[i].B1, B1)

        B2, chart = get_field(ts, it=it, field="B2 c=$i")
        @assert size(B2) == size(bulkevols[i].B2)
        copyto!(bulkevols[i].B2, B2)

        G, chart = get_field(ts, it=it, field="G c=$i")
        @assert size(G) == size(bulkevols[i].G)
        copyto!(bulkevols[i].G, G)

        phi, chart = get_field(ts, it=it, field="phi c=$i")
        @assert size(phi) == size(bulkevols[i].phi)
        copyto!(bulkevols[i].phi, phi)
    end

    nothing
end

# restore all boundary fields
function restore!(boundary::Boundary, ts::OpenPMDTimeSeries, it::Int)

    a4, chart = get_field(ts, it=it, field="a4")
    @assert size(a4) == size(boundary.a4)
    copyto!(boundary.a4, a4)

    fx2, chart = get_field(ts, it=it, field="fx2")
    @assert size(fx2) == size(boundary.fx2)
    copyto!(boundary.fx2, fx2)

    fy2, chart = get_field(ts, it=it, field="fy2")
    @assert size(fy2) == size(boundary.fy2)
    copyto!(boundary.fy2, fy2)

    nothing
end

# restore gauge field
function restore!(gauge::Gauge, ts::OpenPMDTimeSeries, it::Int)
    xi, chart = get_field(ts, it=it, field="xi")
    @assert size(xi) == size(gauge.xi)
    copyto!(gauge.xi, xi)
    nothing
end


# recover all evolved variables from a checkpoint file
function recover(bulkevols::BulkPartition, boundary::Boundary, gauge::Gauge,
                 recovery_dir::String)
    ts = try
        OpenPMDTimeSeries(recovery_dir; prefix="checkpoint_it")
    catch e
        throw(e)
    end
    # grab last iteration of the timeseries, which should be the latest checkpoint
    it = ts.iterations[end]

    restore!(bulkevols, ts, it)
    restore!(boundary, ts, it)
    restore!(gauge, ts, it)

    ts.current_iteration, ts.current_t
end

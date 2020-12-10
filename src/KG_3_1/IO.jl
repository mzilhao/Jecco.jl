
function output_writer(u::EvolVars, chart2D::Chart, charts, tinfo::Jecco.TimeInfo,
                       io::InOut, potential::Potential)
    Nsys = length(charts)

    # output structures
    out_bdry  = Jecco.Output(io.out_dir, "boundary_", tinfo)
    out_bulk  = Jecco.Output(io.out_dir, "bulk_", tinfo)

    boundary  = getboundary(u)
    bulkevols = getbulkevolvedpartition(u)

    #=
    for analysis, it's useful to have the boundary behaviour of the bulk fields,
    so let's output those as well in the boundary file. note that the slicing
    operation [1,:,:] converts 3D arrays into 2D arrays; since we want the same
    shape as in the remaining boundary fields -- 3D arrays of size (1,Nx,Ny) --
    we must reshape them first.
    =#
    phi2 = reshape(getphi(bulkevols[1])[1,:,:], size(chart2D))

    # output fields
    boundary_fields = (
        Jecco.Field("a4",   boundary.a4,  chart2D),
        Jecco.Field("phi2", phi2,         chart2D),
    )
    bulkevols_fields = ntuple(i -> (
        Jecco.Field("phi c=$i", bulkevols[i].phi, charts[i])
    ), Nsys)

    last_output_boundary_t = -io.out_boundary_every_t
    last_output_bulk_t     = -io.out_bulk_every_t

    function (u::EvolVars)
        boundary  = getboundary(u)
        bulkevols = getbulkevolvedpartition(u)

        it = tinfo.it
        tt = tinfo.t

        do_output_boundary = false
        do_output_bulk     = false

        if (io.out_boundary_every > 0 && it % io.out_boundary_every == 0)
            do_output_boundary = true
        end
        if (io.out_bulk_every > 0 && it % io.out_bulk_every == 0)
            do_output_bulk = true
        end

        if (io.out_boundary_every_t > 0 &&
            tt >= last_output_boundary_t + io.out_boundary_every_t - 1e-12)
            do_output_boundary     = true
            last_output_boundary_t = tt
        end
        if (io.out_bulk_every_t > 0 &&
            tt >= last_output_bulk_t + io.out_bulk_every_t - 1e-12)
            do_output_bulk     = true
            last_output_bulk_t = tt
        end

        if do_output_boundary
            boundary_fields[1].data = geta4(boundary)
            @views copyto!(phi2, getphi(bulkevols[1])[1,:,:])
            boundary_fields[2].data = phi2

            # write data
            out_bdry(boundary_fields)
        end

        if do_output_bulk
            @inbounds for i in 1:Nsys
                bulkevols_fields[i][1].data = getphi(bulkevols[i])
            end
            # write data
            out_bulk.(bulkevols_fields)
        end

        nothing
    end
end


function checkpoint_writer(u::EvolVars, chart2D::Chart, charts, tinfo::Jecco.TimeInfo,
                           io::InOut)
    Nsys = length(charts)

    # output structure
    out  = Jecco.Output(io.checkpoint_dir, "checkpoint_it", tinfo)

    boundary  = getboundary(u)
    bulkevols = getbulkevolvedpartition(u)

    # output fields
    boundary_fields = (
        Jecco.Field("a4",   boundary.a4,  chart2D),
    )
    bulkevols_fields = ntuple(i -> (
        Jecco.Field("phi c=$i", getphi(bulkevols[i]), charts[i])
    ), Nsys)

    function (u::EvolVars)
        boundary  = getboundary(u)
        bulkevols = getbulkevolvedpartition(u)

        boundary_fields[1].data = geta4(boundary)

        @inbounds for i in 1:Nsys
            bulkevols_fields[i][1].data = getphi(bulkevols[i])
        end

        # write data
        out(boundary_fields)
        out.(bulkevols_fields)

        nothing
    end
end


# restore all bulk evolved fields
function restore!(bulkevols::BulkPartition{Nsys}, ts::OpenPMDTimeSeries,
                  it::Int) where {Nsys}
    for i in 1:Nsys
        phi, chart = get_field(ts, it=it, field="phi c=$i")
        @assert size(phi) == size(bulkevols[i].phi)
        copyto!(bulkevols[i].phi, phi)
    end

    nothing
end

# restore all boundary fields
function restore!(boundary::Boundary, ts::OpenPMDTimeSeries, it::Int)
    a4GF = geta4(boundary)
    a4, chart = get_field(ts, it=it, field="a4")
    @assert size(a4) == size(a4GF)
    copyto!(a4GF, a4)

    nothing
end


# recover all evolved variables from a checkpoint file
function recover(bulkevols::BulkPartition, boundary::Boundary,
                 recovery_dir::String)
    ts = try
        OpenPMDTimeSeries(recovery_dir; prefix="checkpoint_it")
    catch e
        throw(e)
    end
    # grab last iteration of the timeseries, which should be the latest checkpoint
    it = ts.iterations[end]
    println("INFO: Recovering from checkpoint file $(ts.files[end])")

    restore!(bulkevols, ts, it)
    restore!(boundary, ts, it)

    ts.current_iteration, ts.current_t
end

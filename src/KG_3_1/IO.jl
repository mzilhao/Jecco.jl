
function output_writer(u::EvolVars, chart2D::Chart, atlas, tinfo::Jecco.TimeInfo,
                       io::InOut, potential::Potential)
    Nsys = length(atlas)

    # output structures
    out_bulk  = Jecco.Output(io.out_dir, "bulk_", tinfo)
    bulkevols = getbulkevolvedpartition(u)

    bulkevols_fields = ntuple(i -> (
        Jecco.Field("phi c=$i", getphi(bulkevols[i]), atlas[i])
    ), Nsys)

    last_output_bulk_t     = -io.out_bulk_every_t

    function (u::EvolVars)
        bulkevols = getbulkevolvedpartition(u)

        it = tinfo.it
        tt = tinfo.t

        do_output_bulk     = false

        if (io.out_bulk_every > 0 && it % io.out_bulk_every == 0)
            do_output_bulk = true
        end

        if (io.out_bulk_every_t > 0 &&
            tt >= last_output_bulk_t + io.out_bulk_every_t - 1e-12)
            do_output_bulk     = true
            last_output_bulk_t = tt
        end

        if do_output_bulk
            @inbounds for i in 1:Nsys
                bulkevols_fields[i].data = getphi(bulkevols[i])
            end
            # write data
            out_bulk.(bulkevols_fields)
        end

        nothing
    end
end


function checkpoint_writer(u::EvolVars, chart2D::Chart, atlas, tinfo::Jecco.TimeInfo,
                           io::InOut)
    Nsys = length(atlas)

    # output structure
    out  = Jecco.Output(io.checkpoint_dir, "checkpoint_it", tinfo)

    bulkevols = getbulkevolvedpartition(u)

    # output fields
    bulkevols_fields = ntuple(i -> (
        Jecco.Field("phi c=$i", getphi(bulkevols[i]), atlas[i])
    ), Nsys)

    function (u::EvolVars)
        bulkevols = getbulkevolvedpartition(u)

        @inbounds for i in 1:Nsys
            bulkevols_fields[i].data = getphi(bulkevols[i])
        end

        # write data
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

# recover all necessary values
function recover(bulkevols::BulkPartition, boundary::Boundary, systems::SystemPartition,
                 id::InitialData, recovery_dir::String)
    ts = try
        OpenPMDTimeSeries(recovery_dir; prefix="checkpoint_it")
    catch e
        throw(e)
    end
    # grab last iteration of the timeseries, which should be the latest checkpoint
    it = ts.iterations[end]
    println("INFO: Recovering from checkpoint file $(ts.files[end])")

    restore!(bulkevols, ts, it)
    # for the boundary data
    init_data!(boundary, systems[1], id)

    ts.current_iteration, ts.current_t
end

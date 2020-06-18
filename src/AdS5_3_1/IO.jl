
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

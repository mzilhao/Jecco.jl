
function output_writer(u::EvolVars, chart2D::Chart, charts, tinfo::Jecco.TimeInfo,
                       io::InOut, potential::Potential, phi0)
    Nsys = length(charts)

    # output structures
    out_bdry  = Jecco.Output(io.out_dir, "boundary_", tinfo)
    out_gauge = Jecco.Output(io.out_dir, "gauge_", tinfo)
    out_bulk  = Jecco.Output(io.out_dir, "bulk_", tinfo)

    boundary  = getboundary(u)
    gauge     = getgauge(u)
    bulkevols = getbulkevolvedpartition(u)

    a4  = geta4(boundary)
    fx2 = getfx2(boundary)
    fy2 = getfy2(boundary)
    xi  = getxi(gauge)

    # NamedTuple with potential parameters
    params = parameters(potential)
    params = merge((phi0=phi0,), params)

    #=
    for analysis, it's useful to have the boundary behaviour of the bulk fields,
    so let's output those as well in the boundary file. note that the slicing
    operation [1,:,:] converts 3D arrays into 2D arrays; since we want the same
    shape as in the remaining boundary fields -- 3D arrays of size (1,Nx,Ny) --
    we must reshape them first.
    =#
    B1_1  = getB1(bulkevols[1])
    B2_1  = getB2(bulkevols[1])
    G_1   = getG(bulkevols[1])
    phi_1 = getphi(bulkevols[1])
    b14   = reshape(B1_1[1,:,:],  size(chart2D))
    b24   = reshape(B2_1[1,:,:],  size(chart2D))
    g4    = reshape(G_1[1,:,:],   size(chart2D))
    phi3  = reshape(phi_1[1,:,:], size(chart2D))

    # phi2 = phi0^3 phi3 - phi0 xi^2
    phi2 = phi0*phi0*phi0 .* phi3 .- phi0 .* xi .* xi

    # output fields
    boundary_fields = (
        Jecco.Field("a4",   a4,  chart2D),
        Jecco.Field("fx2",  fx2, chart2D),
        Jecco.Field("fy2",  fy2, chart2D),
        Jecco.Field("b14",  b14, chart2D),
        Jecco.Field("b24",  b24, chart2D),
        Jecco.Field("g4",   g4,  chart2D),
        Jecco.Field("phi2", phi2,chart2D),
    )
    gauge_fields = Jecco.Field("xi", xi, chart2D)
    bulkevols_fields = ntuple(i -> (
        Jecco.Field("B1 c=$i",  getB1(bulkevols[i]),  charts[i]),
        Jecco.Field("B2 c=$i",  getB2(bulkevols[i]),  charts[i]),
        Jecco.Field("G c=$i",   getG(bulkevols[i]),   charts[i]),
        Jecco.Field("phi c=$i", getphi(bulkevols[i]), charts[i])
    ), Nsys)

    last_output_boundary_t = -io.out_boundary_every_t
    last_output_gauge_t    = -io.out_gauge_every_t
    last_output_bulk_t     = -io.out_bulk_every_t

    function (u::EvolVars)
        boundary  = getboundary(u)
        gauge     = getgauge(u)
        bulkevols = getbulkevolvedpartition(u)

        it = tinfo.it
        tt = tinfo.t

        do_output_boundary = false
        do_output_gauge    = false
        do_output_bulk     = false

        if (io.out_boundary_every > 0 && it % io.out_boundary_every == 0)
            do_output_boundary = true
        end
        if (io.out_gauge_every > 0 && it % io.out_gauge_every == 0)
            do_output_gauge = true
        end
        if (io.out_bulk_every > 0 && it % io.out_bulk_every == 0)
            do_output_bulk = true
        end

        if (io.out_boundary_every_t > 0 &&
            tt >= last_output_boundary_t + io.out_boundary_every_t - 1e-12)
            do_output_boundary     = true
            last_output_boundary_t = tt
        end
        if (io.out_gauge_every_t > 0 &&
            tt >= last_output_gauge_t + io.out_gauge_every_t - 1e-12)
            do_output_gauge     = true
            last_output_gauge_t = tt
        end
        if (io.out_bulk_every_t > 0 &&
            tt >= last_output_bulk_t + io.out_bulk_every_t - 1e-12)
            do_output_bulk     = true
            last_output_bulk_t = tt
        end

        if do_output_boundary
            boundary_fields[1].data = geta4(boundary)
            boundary_fields[2].data = getfx2(boundary)
            boundary_fields[3].data = getfy2(boundary)

            B1_1  = getB1(bulkevols[1])
            B2_1  = getB2(bulkevols[1])
            G_1   = getG(bulkevols[1])
            phi_1 = getphi(bulkevols[1])
            @views copyto!(boundary_fields[4].data, B1_1[1,:,:])
            @views copyto!(boundary_fields[5].data, B2_1[1,:,:])
            @views copyto!(boundary_fields[6].data, G_1[1,:,:])

            @views copyto!(phi2, phi_1[1,:,:])
            # phi2 = phi0^3 phi[1,:,:] - phi0 xi^2
            phi2 .= phi0*phi0*phi0 .* phi2 .- phi0 .* xi .* xi

            boundary_fields[7].data = phi2

            # write data
            out_bdry(boundary_fields, params=params)
        end

        if do_output_gauge
            gauge_fields.data = getxi(gauge)
            # write data
            out_gauge(gauge_fields, params=params)
        end

        if do_output_bulk
            @inbounds for i in 1:Nsys
                bulkevols_fields[i][1].data = getB1(bulkevols[i])
                bulkevols_fields[i][2].data = getB2(bulkevols[i])
                bulkevols_fields[i][3].data = getG(bulkevols[i])
                bulkevols_fields[i][4].data = getphi(bulkevols[i])
            end
            # write data
            out_bulk.(bulkevols_fields, params=params)
        end

        nothing
    end
end


function output_writer(bulkconstrains::BulkPartition{Nsys,BulkConstrained{T}}, charts,
                       tinfo::Jecco.TimeInfo, io::InOut, potential::Potential,
                       phi0) where {Nsys,T}
    @assert Nsys == length(charts)

    # NamedTuple with potential parameters
    params = parameters(potential)
    params = merge((phi0=phi0,), params)

    # output structure
    out  = Jecco.Output(io.out_dir, "constrained_", tinfo)

    # output fields
    fields = ntuple(i -> (
        Jecco.Field("S c=$i",    getS(bulkconstrains[i]),    charts[i]),
        Jecco.Field("Fx c=$i",   getFx(bulkconstrains[i]),   charts[i]),
        Jecco.Field("Fy c=$i",   getFy(bulkconstrains[i]),   charts[i]),
        Jecco.Field("B1d c=$i",  getB1d(bulkconstrains[i]),  charts[i]),
        Jecco.Field("B2d c=$i",  getB2d(bulkconstrains[i]),  charts[i]),
        Jecco.Field("Gd c=$i",   getGd(bulkconstrains[i]),   charts[i]),
        Jecco.Field("phid c=$i", getphid(bulkconstrains[i]), charts[i]),
        Jecco.Field("Sd c=$i",   getSd(bulkconstrains[i]),   charts[i]),
        Jecco.Field("A c=$i",    getA(bulkconstrains[i]),    charts[i])
    ), Nsys)

    last_output_t = -io.out_bulkconstrained_every_t

    function (bulkconstrains)
        it = tinfo.it
        tt = tinfo.t

        do_output = false

        if (io.out_bulkconstrained_every > 0 && it % io.out_bulkconstrained_every == 0)
            do_output = true
        end

        if (io.out_bulkconstrained_every_t > 0 &&
            tt >= last_output_t + io.out_bulkconstrained_every_t - 1e-12)
            do_output = true
            last_output_t = tt
        end

        if do_output
            @inbounds for i in 1:Nsys
                fields[i][1].data = getS(bulkconstrains[i])
                fields[i][2].data = getFx(bulkconstrains[i])
                fields[i][3].data = getFy(bulkconstrains[i])
                fields[i][4].data = getB1d(bulkconstrains[i])
                fields[i][5].data = getB2d(bulkconstrains[i])
                fields[i][6].data = getGd(bulkconstrains[i])
                fields[i][7].data = getphid(bulkconstrains[i])
                fields[i][8].data = getSd(bulkconstrains[i])
                fields[i][9].data = getA(bulkconstrains[i])
            end
            # write data
            out.(fields, params=params)
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
    gauge     = getgauge(u)
    bulkevols = getbulkevolvedpartition(u)

    a4  = geta4(boundary)
    fx2 = getfx2(boundary)
    fy2 = getfy2(boundary)
    xi  = getxi(gauge)

    # output fields
    boundary_fields = (
        Jecco.Field("a4",   a4,  chart2D),
        Jecco.Field("fx2",  fx2, chart2D),
        Jecco.Field("fy2",  fy2, chart2D),
    )
    gauge_fields = Jecco.Field("xi", xi, chart2D)
    bulkevols_fields = ntuple(i -> (
        Jecco.Field("B1 c=$i",  getB1(bulkevols[i]),  charts[i]),
        Jecco.Field("B2 c=$i",  getB2(bulkevols[i]),  charts[i]),
        Jecco.Field("G c=$i",   getG(bulkevols[i]),   charts[i]),
        Jecco.Field("phi c=$i", getphi(bulkevols[i]), charts[i])
    ), Nsys)

    function (u::EvolVars)
        boundary  = getboundary(u)
        gauge     = getgauge(u)
        bulkevols = getbulkevolvedpartition(u)

        boundary_fields[1].data = geta4(boundary)
        boundary_fields[2].data = getfx2(boundary)
        boundary_fields[3].data = getfy2(boundary)

        gauge_fields.data = getxi(gauge)

        @inbounds for i in 1:Nsys
            bulkevols_fields[i][1].data = getB1(bulkevols[i])
            bulkevols_fields[i][2].data = getB2(bulkevols[i])
            bulkevols_fields[i][3].data = getG(bulkevols[i])
            bulkevols_fields[i][4].data = getphi(bulkevols[i])
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
        B1GF  = getB1(bulkevols[i])
        B2GF  = getB2(bulkevols[i])
        GGF   = getG(bulkevols[i])
        phiGF = getphi(bulkevols[i])

        B1, chart = get_field(ts, it=it, field="B1 c=$i")
        @assert size(B1) == size(B1GF)
        copyto!(B1GF, B1)

        B2, chart = get_field(ts, it=it, field="B2 c=$i")
        @assert size(B2) == size(B2GF)
        copyto!(B2GF, B2)

        G, chart = get_field(ts, it=it, field="G c=$i")
        @assert size(G) == size(GGF)
        copyto!(GGF, G)

        phi, chart = get_field(ts, it=it, field="phi c=$i")
        @assert size(phi) == size(phiGF)
        copyto!(phiGF, phi)
    end

    nothing
end

# restore all boundary fields
function restore!(boundary::Boundary, ts::OpenPMDTimeSeries, it::Int)
    a4GF  = geta4(boundary)
    fx2GF = getfx2(boundary)
    fy2GF = getfy2(boundary)

    a4, chart = get_field(ts, it=it, field="a4")
    @assert size(a4) == size(a4GF)
    copyto!(a4GF, a4)

    fx2, chart = get_field(ts, it=it, field="fx2")
    @assert size(fx2) == size(fx2GF)
    copyto!(fx2GF, fx2)

    fy2, chart = get_field(ts, it=it, field="fy2")
    @assert size(fy2) == size(fy2GF)
    copyto!(fy2GF, fy2)

    nothing
end

# restore gauge field
function restore!(gauge::Gauge, ts::OpenPMDTimeSeries, it::Int)
    xiGF  = getxi(gauge)
    xi, chart = get_field(ts, it=it, field="xi")
    @assert size(xi) == size(xiGF)
    copyto!(xiGF, xi)
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
    println("INFO: Recovering from checkpoint file $(ts.files[end])")

    restore!(bulkevols, ts, it)
    restore!(boundary, ts, it)
    restore!(gauge, ts, it)

    ts.current_iteration, ts.current_t
end

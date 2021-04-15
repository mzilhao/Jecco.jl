
@time @testset "Homogeneous Black brane tests:" begin

    grid = SpecCartGrid3D(
        x_min            = -5.0,
        x_max            =  5.0,
        x_nodes          =  6,
        y_min            = -5.0,
        y_max            =  5.0,
        y_nodes          =  6,
        u_outer_min      =  0.1,
        u_outer_max      =  1.003,
        u_outer_domains  =  4,
        u_outer_nodes    =  32,
        u_inner_nodes    =  12,
    )

    id = BlackBrane()

    evoleq = AffineNull(
        phi0       = 0.0,
        potential  = ZeroPotential(),
        gaugecondition = AdS5_3_1.ConstantGauge(),
    )

    # atlas of grid configuration and respective SystemPartition
    atlas     = Atlas(grid)
    systems   = SystemPartition(grid)

    # allocate variables
    boundary  = Boundary(grid)
    gauge     = Gauge(grid)
    gauge_t   = similar(gauge)

    bulkevols      = BulkEvolvedPartition(grid)
    bulkconstrains = BulkConstrainedPartition(grid)
    bulkderivs     = BulkDerivPartition(grid)
    horizoncache   = AdS5_3_1.HorizonCache(systems[end], evoleq.gaugecondition.fd_order)

    # initial conditions
    id(bulkconstrains, bulkevols, bulkderivs, boundary, gauge,
       horizoncache, systems, evoleq)

    # function to solve the nested system, given the initial data
    nested = AdS5_3_1.Nested(systems, bulkconstrains, bulkderivs)

    # solve nested system for the constrained variables
    nested(bulkevols, boundary, gauge, evoleq)

    # analyze data

    Du_in  = systems[1].Du
    Du_out = systems[2].Du

    # inner grid
    bulk_in  = bulkconstrains[1]
    chart_in = atlas[1]
    uu_in    = chart_in.coords[1][:]

    # first outer grid
    bulk_out_1  = bulkconstrains[2]
    chart_out_1 = atlas[2]
    uu_out_1    = chart_out_1.coords[1][:]

    # take u-derivatives
    Du_A_in    = Du_in * bulk_in.A
    Du_A_out_1 = Du_out * bulk_out_1.A


    # Compare value of A and its u-derivative on the inner-outer domain interface

    xi0         = gauge.xi[1,1,1]
    A0_in       = bulk_in.A[:,1,1]
    A0_out      = bulk_out_1.A[:,1,1]
    Du_A0_in    = Du_A_in[:,1,1]
    Du_A0_out_1 = Du_A_out_1[:,1,1]

    Aout_in_end = 1/uu_in[end]^2 + 2*xi0/uu_in[end] + xi0^2 + uu_in[end]^2 * A0_in[end]

    Du_Aout_in_end = -2/uu_in[end]^3 - 2*xi0/uu_in[end]^2 + 2*uu_in[end] * A0_in[end] +
        uu_in[end]^2 * Du_A0_in[end]

    # in-out interface
    @test Aout_in_end    ≈ A0_out[1]
    @test Du_Aout_in_end ≈ Du_A0_out_1[1]


    # remaining outer grids
    for i in 2:grid.u_outer_domains
        Du_1    = systems[i].Du
        bulk_1  = bulkconstrains[i]
        chart_1 = atlas[i]
        uu_1    = chart_1.coords[1][:]

        Du_2    = systems[i+1].Du
        bulk_2  = bulkconstrains[i+1]
        chart_2 = atlas[i+1]
        uu_2    = chart_2.coords[1][:]

        # take u-derivatives
        Du_A_1    = Du_1 * bulk_1.A
        Du_A_2    = Du_2 * bulk_2.A

        # 1-2 interface
        @test all(bulk_1.A[end,:,:] .≈ bulk_2.A[1,:,:])
        @test all(Du_A_1[end,:,:]   .≈ Du_A_2[1,:,:])
    end


    # now check if the result matches the expected value. since the
    # configuration is homogeneous, let's compare at
    i = 1; j = 1

    # first let's check inner grid
    A_in    = bulk_in.A[:,i,j]
    A_u_in  = 2 .* uu_in .* bulk_in.A[:,i,j] .+ uu_in .* uu_in .* Du_A_in[:,i,j]

    a4 = -id.energy_dens/0.75
    A_in_exact = a4 ./ (1 .+ 2 .* xi0 .* uu_in .+ xi0^2 .* uu_in .* uu_in)

    A_u_in_exact = 2 .* uu_in .* a4 ./ (1 .+ xi0 .* uu_in).^3

    @test A_in    ≈ A_in_exact
    @test A_u_in  ≈ A_u_in_exact


    # now for the outer grids
    for i in 2:grid.u_outer_domains+1
        Du    = systems[i].Du
        bulk  = bulkconstrains[i]
        chart = atlas[i]
        uu    = chart.coords[1][:]

        # take u-derivative
        Du_A  = Du * bulk.A

        A     = bulk.A[:,i,j]
        A_u   = Du_A[:,i,j]

        A_exact  = 1 ./ uu.^2 .+ 2 .* xi0 ./ uu .+ xi0^2 .+
            uu.^2 .* a4 ./ (1 .+ 2 .* xi0 .* uu .+ xi0^2 .* uu.^2)

        A_u_exact = 2 .* (.- 1 .- uu .* xi0 + a4 .* uu.^4 ./ (1 .+ uu .* xi0).^3 ) ./ uu.^3

        @test A   ≈ A_exact
        @test A_u ≈ A_u_exact

    end

end

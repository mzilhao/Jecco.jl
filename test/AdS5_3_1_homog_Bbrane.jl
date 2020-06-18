using Test

using Jecco
using Jecco.AdS5_3_1

@time @testset "Homogeneous Black brane tests:" begin

    out_doms = 4

    grid = SpecCartGrid3D(
        x_min            = -5.0,
        x_max            =  5.0,
        x_nodes          =  6,
        y_min            = -5.0,
        y_max            =  5.0,
        y_nodes          =  6,
        u_outer_min      =  0.1,
        u_outer_max      =  1.0,
        u_outer_domains  =  out_doms,
        u_outer_nodes    =  30,
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

    # initial conditions
    init_data!(bulkevols, boundary, gauge, systems, id)

    # function to solve the nested system, given the initial data
    nested = Nested(systems, bulkconstrains, bulkderivs)

    # solve nested system for the constrained variables
    nested(bulkevols, boundary, gauge, evoleq)

    #################
    # analyze data  #
    #################

    # gauge
    xi     = gauge.xi[1,1,1] #due to homogeneity; change for non-homog.
    Du_in  = systems[1].Du
    Du_out = systems[2].Du

    # inner grid
    i = 1

    bulk_in = Bulk(bulkevols[i], bulkconstrains[i])

    chart_in = atlas.charts[i]
    uu_in    = chart_in.coords[1][:]

    A_in  = bulk_in.A[:,1,1];
    Au_in = Du_in*A_in;

    # outer grids
    num_out_domains = grid.u_outer_domains

    # FIXME type stability
    bulk_out  = []
    chart_out = []
    uu_out    = []

    A_out          = []
    Au_out         = []

    uu_out_full    = []
    A_out_full_num = []

    for i in 2:num_out_domains+1

        bulk_out_i = Bulk(bulkevols[i], bulkconstrains[i])
        append!(bulk_out, [bulk_out_i])

        chart_out_i = atlas.charts[i]
        append!(chart_out, [chart_out_i])

        uu_out_i    = chart_out_i.coords[1][:]
        append!(uu_out, [uu_out_i])

        A_out_i = bulk_out_i.A[:,1,1]
        append!(A_out,  [A_out_i])
        append!(Au_out, [Du_out*A_out_i])

    end

    #################################################
    # Compare value of A on the interface of domains#
    #################################################

    A_in_end = (1/(uu_in[end]^2) + 2.0*xi/uu_in[end] +   xi^2 + uu_in[end]^2 * A_in[end])
    Au_in_end = (-2/(uu_in[end]^3) - 2.0*xi/(uu_in[end]^2)  + 2.0* uu_in[end]* A_in[end] + (uu_in[end]^2)* Au_in[end])

    # in-out interface
    @test A_in_end  ≈ A_out[1][1]
    @test Au_in_end ≈ Au_out[1][1]

    # mutliple out domains interface
    if num_out_domains>=2
        for i in 2:num_out_domains
            @test A_out[i-1][end]  ≈ A_out[i][1]
            @test Au_out[i-1][end] ≈ Au_out[i][1]
        end
    end

    #######################################################
    # compare the function in,out domains and full domain #
    # be careful of redefs in different domains           #
    #######################################################

    # comparison inner grid only
    A_in_reg_num = zeros(length(A_in))
    A_in_reg_exact = zeros(length(A_in))
    Au_in_reg_num = zeros(length(A_in))
    Au_in_reg_exact = zeros(length(A_in))
    a4 = -id.energy_dens/0.75

    # A regularized inner grid
    for i in 1:length(A_in)
        u = uu_in[i]
        A = A_in[i]
        Au = Au_in[i]
        A_in_reg_num[i]    = xi^2 + u^2 * A
        Au_in_reg_num[i]   = 2.0*u * A + u^2 * Au
        A_in_reg_exact[i]  = xi^2 + u^2 * a4 /(1.0 + 2.0* xi *u + xi^2 *u^2 )
        Au_in_reg_exact[i] =  (2.0 * u * a4)/((1.0 + xi *u)^3)
    end

    @test A_in_reg_num  ≈ A_in_reg_exact
    @test Au_in_reg_num ≈ Au_in_reg_exact

    # comparison outer grids only
    uu_out_full    = []
    A_out_full_num = []
    Au_out_full_num = []
    for i in 1:num_out_domains
        append!(uu_out_full, uu_out[i])
        append!(A_out_full_num, A_out[i])
        append!(Au_out_full_num, Au_out[i])
    end

    #exact
    A_out_full_exact = zeros(length(A_out_full_num))
    Au_out_full_exact = zeros(length(A_out_full_num))
    for i in 1:length(A_out_full_num)
        u = uu_out_full[i]
        A_out_full_exact[i] = 1/(u^2) + 2.0*xi/u + xi^2 + u^2 * a4/(1.0 + 2.0* xi *u + xi^2 *u^2)
        Au_out_full_exact[i] = -2/(u^3) - 2.0*xi/(u^2)  - u^2 * a4 * (2*xi + 2*u*xi)/((1.0 + 2.0* xi *u + xi^2 *u^2)^2) + 2*a4*u/(1.0 + 2.0* xi *u + xi^2 *u^2)
    end

    @test A_out_full_num  ≈ A_out_full_exact
    #@test Au_out_full_num ≈ Au_out_full_exact atol = 1.e-2

end

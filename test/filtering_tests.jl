
@testset "KO Filtering tests:" begin

    # 1D case

    # N=3
    dim = 20
    delta_    = zeros(dim)
    delta_[3] = 1.0

    order = 3
    sigma = 16.0
    ko_filter = Jecco.KO_Filter(order, sigma, dim)

    ko_filter(delta_)

    @test delta_[1:5] ≈ ko_filter.kernel

    # N=5
    dim = 20
    delta_    = zeros(dim)
    delta_[4] = 1.0

    order = 5
    sigma = 64.0
    ko_filter = Jecco.KO_Filter(order, sigma, dim)

    ko_filter(delta_)

    @test delta_[1:7] ≈ ko_filter.kernel


    # 3D case

    dimx = 12
    dimy = 16
    dimz = 10
    delta_         = zeros(dimx, dimy, dimz)
    delta_[:,4,:] .= 1

    order = 5
    sigma = 64.0
    ko_filter = Jecco.KO_Filter{2}(order, sigma, dimx, dimy, dimz)

    ko_filter(delta_)

    myarr = []
    for k in 1:dimz
        for i in 1:dimx
            push!(myarr, delta_[i,1:7,k] ≈ ko_filter.kernel)
        end
    end

    @test all(myarr)

end


@testset "Exponential filtering tests:" begin

    # kernel test

    dim = 49
    γ   = 8.0

    exp_filter = Jecco.Exp_Filter{1}(γ, dim)

    kernel_bench = [ 1.0, 9.99999999999889e-01, 9.99999999377987e-01, 9.99997121269861e-01,
                     9.88950030743732e-01, 2.22034255471165e-16]

    @test exp_filter.kernel[end-5:end] ≈ kernel_bench

end


@testset "KO Dissipation tests:" begin

    # N=3
    DKO = Jecco.KO_Centered{1}(3, 16.0, 0.1, 20)
    @test DKO.stencil_coefs ≈ [-10.0, 40.0, -60.0, 40.0, -10.0]

    delta_    = zeros(20)
    delta_[3] = 1.0
    bla = DKO * delta_
    @test bla[1:5] ≈ DKO.stencil_coefs


    # N=5
    DKO = Jecco.KO_Centered{1}(5, 64.0, 0.1, 20)
    @test DKO.stencil_coefs ≈ [10.0, -60.0, 150.0, -200.0, 150.0, -60.0, 10.0]

    delta_    = zeros(20)
    delta_[4] = 1.0
    bla = DKO * delta_
    @test bla[1:7] ≈ DKO.stencil_coefs


    # N=7
    DKO = Jecco.KO_Centered{1}(7, 256.0, 0.1, 20)
    @test DKO.stencil_coefs ≈ [-10.0, 80.0, -280.0, 560.0, -700.0, 560.0, -280.0, 80.0, -10.0]

    delta_    = zeros(20)
    delta_[5] = 1.0
    bla = DKO * delta_
    @test bla[1:9] ≈ DKO.stencil_coefs


    # 3D case

    dimx = 12
    dimy = 16
    dimz = 10
    delta_         = zeros(dimx, dimy, dimz)
    delta_[:,4,:] .= 1

    order = 5
    sigma = 64.0
    DKO = Jecco.KO_Centered{2}(order, sigma, 1.0, dimy)

    bla = DKO * delta_

    myarr = []
    for k in 1:dimz
        for i in 1:dimx
            push!(myarr, bla[i,1:7,k] ≈ DKO.stencil_coefs)
        end
    end
    @test all(myarr)

end


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
    exp_filter   = Jecco.Exp_Filter(γ, dim)
    kernel_bench = [ 1.0, 9.99999999999889e-01, 9.99999999377987e-01, 9.99997121269861e-01,
                     9.88950030743732e-01, 2.22034255471165e-16]
    @test exp_filter.kernel[end-5:end] ≈ kernel_bench


    # 1D case
    dim = 16
    delta_    = zeros(dim)
    delta_[1] = 1.0

    exp_filter1D = Jecco.Exp_Filter{1}(γ, dim)

    f = copy(delta_)

    exp_filter1D(f)

    @test exp_filter1D.fft_plan * f ≈ exp_filter1D.kernel


    # 3D case
    dimx = 12
    dimy = 16
    dimz = 8
    delta3D_ = zeros(dimx, dimy, dimz)
    delta3D_[:,1,:] .= 1

    exp_filter3D = Jecco.Exp_Filter{2}(γ, dimx, dimy, dimz)

    f = copy(delta3D_)

    exp_filter3D(f)

    tmp = exp_filter3D.fft_plan * f
    myarr = []
    for k in 1:dimz
        for i in 1:dimx
            push!(myarr, tmp[i,:,k] ≈ exp_filter3D.kernel)
        end
    end

    @test all(myarr)
end

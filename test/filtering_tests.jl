
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

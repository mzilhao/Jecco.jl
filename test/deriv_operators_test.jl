
@testset "FD Derivative tests:" begin

    # 1D case
    xmin   = -2.0*pi
    xmax   =  2.0*pi
    xnodes =  600
    ord    =  4

    hx    = (xmax - xmin) / xnodes

    xx = collect(xmin:hx:xmax-hx)
    f  = sin.(xx)

    D1 = CenteredDiff(1, ord, hx, length(xx))
    D2 = CenteredDiff(2, ord, hx, length(xx))

    df1 = D1 * f
    df2 = D2 * f

    @test df1 ≈ cos.(xx) atol=hx^ord
    @test d2f ≈ -f atol=hx^ord
end

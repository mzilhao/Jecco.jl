
@testset "Grids tests:" begin

    umin   = 0.0
    umax   = 2.0
    unodes = 16

    ucoord = SpectralCoord("u", umin, umax, unodes)
    u, = Jecco.cheb(umin, umax, unodes)
    @test ucoord[:] ≈ u

    xmin   = -10.0
    xmax   =  10.0
    xnodes =  20

    xcoord = CartCoord{2}("x", xmin, xmax, xnodes, endpoint=false)
    @test xcoord[:] == collect(-10.0:1.0:9.0)

    ymin   = -20.0
    ymax   =  20.0
    ynodes =  20

    ycoord = CartCoord{3}("y", ymin, ymax, ynodes, endpoint=false)

    grid = Grid(ucoord, xcoord, ycoord)

    u, x, y = grid[:]

    @test x == collect(-10.0:1.0:9.0)
    @test y == collect(-20.0:2.0:19.0)

end

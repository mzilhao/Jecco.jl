
@testset "FD Derivative tests:" begin

    Dx = CenteredDiff{1}(1, 4, 1, 10)
    @test 12*Dx[1,:]  == [0.0, 8.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -8.0]
    @test 12*Dx[2,:]  == [-8.0, 0.0, 8.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
    @test 12*Dx[3,:]  == [1.0, -8.0, 0.0, 8.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    @test 12*Dx[4,:]  == [0.0, 1.0, -8.0, 0.0, 8.0, -1.0, 0.0, 0.0, 0.0, 0.0]
    @test 12*Dx[5,:]  == [0.0, 0.0, 1.0, -8.0, 0.0, 8.0, -1.0, 0.0, 0.0, 0.0]
    @test 12*Dx[6,:]  == [0.0, 0.0, 0.0, 1.0, -8.0, 0.0, 8.0, -1.0, 0.0, 0.0]
    @test 12*Dx[7,:]  == [0.0, 0.0, 0.0, 0.0, 1.0, -8.0, 0.0, 8.0, -1.0, 0.0]
    @test 12*Dx[8,:]  == [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -8.0, 0.0, 8.0, -1.0]
    @test 12*Dx[9,:]  == [-1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -8.0, 0.0, 8.0]
    @test 12*Dx[10,:] == [8.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -8.0, 0.0]

    # 1D case

    xmin   = -2.0*pi
    xmax   =  2.0*pi
    xnodes =  600
    ord    =  4

    hx     = (xmax - xmin) / xnodes

    x  = collect(xmin:hx:xmax-hx)
    f  = sin.(x)

    D1 = CenteredDiff(1, ord, hx, length(x))
    D2 = CenteredDiff(2, ord, hx, length(x))

    df  = D1 * f
    d2f = D2 * f

    @test df  ≈ cos.(x) atol=hx^ord
    @test d2f ≈ -f atol=hx^ord

    # now for the callable, point-wise, methods
    for i in eachindex(f)
        @test D1(f,i) == df[i]
    end

    # 3D case

    ymin   = -1.0*pi
    ymax   =  1.0*pi
    ynodes =  20
    zmin   = -1.0*pi
    zmax   =  1.0*pi
    znodes =  300

    hy     = (ymax - ymin) / ynodes
    hz     = (zmax - zmin) / znodes

    y      = collect(ymin:hy:ymax-hy)
    z      = collect(zmin:hz:zmax-hz)

    f      = [sin.(x1) .* sin.(x2) .* sin.(x3) for x1 in x, x2 in y, x3 in z]
    dxf0   = [cos.(x1) .* sin.(x2) .* sin.(x3) for x1 in x, x2 in y, x3 in z]
    dyf0   = [sin.(x1) .* cos.(x2) .* sin.(x3) for x1 in x, x2 in y, x3 in z]
    dzf0   = [sin.(x1) .* sin.(x2) .* cos.(x3) for x1 in x, x2 in y, x3 in z]
    dxzf0  = [cos.(x1) .* sin.(x2) .* cos.(x3) for x1 in x, x2 in y, x3 in z]

    Dx     = CenteredDiff{1}(1, ord, hx, length(x))
    Dz     = CenteredDiff{3}(1, ord, hz, length(z))

    Dxx    = CenteredDiff{1}(2, ord, hx, length(x))
    Dzz    = CenteredDiff{3}(2, ord, hz, length(z))

    dxf    = Dx * f
    dzf    = Dz * f
    dxzf   = Dx * (Dz * f)

    @test dxf  ≈ dxf0
    @test dzf  ≈ dzf0
    @test dxzf ≈ dxzf0

    d2xf   = Dxx * f
    d2zf   = Dzz * f

    @test d2xf ≈ -f
    @test d2zf ≈ -f


    # now for callable, point-wise, methods
    @test Dx(f,2,10,120)  == dxf[2,10,120]
    @test Dx(f,42,20,300) == dxf[42,20,300]
    @test Dz(f,2,10,120)  == dzf[2,10,120]
    @test Dz(f,42,20,300) == dzf[42,20,300]
end

@testset "Spectral Derivative tests:" begin

    # 1D case
    xmin   = -2.0
    xmax   =  2.0
    xnodes =  32

    x, = Jecco.cheb(xmin, xmax, xnodes)
    f = 0.5 * x.^2

    Dx  = ChebDeriv(1, xmin, xmax, xnodes)
    Dxx = ChebDeriv(2, xmin, xmax, xnodes)

    dxf  = Dx * f
    dxxf = Dxx * f
    @test dxf  ≈ x
    @test dxxf ≈ fill(1.0, size(dxxf))

    # now for the callable, point-wise, methods
    for i in eachindex(f)
        @test Dx(f,i) ≈ dxf[i]
    end


    # 3D case

    ymin   = -1.0
    ymax   =  1.0
    ynodes =  8

    zmin   = -1.0
    zmax   =  1.0
    znodes =  16

    y, = Jecco.cheb(ymin, ymax, ynodes)
    z, = Jecco.cheb(zmin, zmax, znodes)

    Dy  = ChebDeriv{2}(1, ymin, ymax, ynodes)
    Dyy = ChebDeriv{2}(2, ymin, ymax, ynodes)
    Dz  = ChebDeriv{3}(1, zmin, zmax, znodes)
    Dzz = ChebDeriv{3}(2, zmin, zmax, znodes)

    f     = [0.5 * x1.^2 .* cos.(x2) .* sin.(x3) for x1 in x, x2 in y, x3 in z]
    dxf0  = [x1 .* cos.(x2) .* sin.(x3)          for x1 in x, x2 in y, x3 in z]
    dzf0  = [0.5 * x1.^2 .* cos.(x2) .* cos.(x3) for x1 in x, x2 in y, x3 in z]
    dxxf0 = [cos.(x2) .* sin.(x3)                for x1 in x, x2 in y, x3 in z]
    dzzf0 = -copy(f)

    dxf = Dx * f
    dzf = Dz * f
    d2xf = Dxx * f
    d2zf = Dzz * f

    @test dxf ≈ dxf0
    @test dzf ≈ dzf0

    @test d2xf ≈ dxxf0
    @test d2zf ≈ dzzf0

    # now for callable, point-wise, methods
    @test Dx(f,2,4,16)  ≈ dxf[2,4,16]
    @test Dx(f,1,6,12)  ≈ dxf[1,6,12]
    @test Dz(f,2,3,8)   ≈ dzf[2,3,8]
    @test Dz(f,16,8,1)  ≈ dzf[16,8,1]
end

@testset "Cross FD derivative tests:" begin
    # 2D FD case
    xmin   = -2.0*pi
    xmax   =  2.0*pi
    xnodes =  600
    ymin   = -1.0*pi
    ymax   =  1.0*pi
    ynodes =  300
    ord    =  4

    hx     = (xmax - xmin) / xnodes
    hy     = (ymax - ymin) / ynodes

    x      = collect(xmin:hx:xmax-hx)
    y      = collect(ymin:hy:ymax-hy)

    f      = [sin.(x1) .* sin.(x2) for x1 in x, x2 in y]

    Dx     = CenteredDiff{1}(1, ord, hx, length(x))
    Dy     = CenteredDiff{2}(1, ord, hy, length(y))

    dxyf   = Dx * (Dy * f)

    @test Dx(Dy, f,  2,120) ≈ dxyf[2,120]
    @test Dx(Dy, f, 42,300) ≈ dxyf[42,300]
end

@testset "FD cross derivative for general arrays tests:" begin

    xmin   = -2.0*pi
    xmax   =  2.0*pi
    xnodes =  600
    ord    =  4

    ymin   = -1.0
    ymax   =  1.0
    ynodes =  16

    zmin   = -1.0*pi
    zmax   =  1.0*pi
    znodes =  300


    hx     = (xmax - xmin) / xnodes
    hz     = (zmax - zmin) / znodes

    x      = collect(xmin:hx:xmax-hx)
    y,     = Jecco.cheb(ymin, ymax, ynodes)
    z      = collect(zmin:hz:zmax-hz)

    f   = [0.5 * sin.(x1) .* x2.^2 .* sin.(x3)  for x1 in x, x2 in y, x3 in z]

    Dx  = CenteredDiff{1}(1, ord, hx, length(x))
    Dz  = CenteredDiff{3}(1, ord, hz, length(z))

    dxzf  = Dx * (Dz * f)

    @test Dx(Dz, f, 2,16,1)     ≈ dxzf[2,16,1]
    @test Dx(Dz, f, 1,12,2)     ≈ dxzf[1,12,2]
    @test Dx(Dz, f, 100,12,300) ≈ dxzf[100,12,300]
    @test Dx(Dz, f, 600,8,100)  ≈ dxzf[600,8,100]


    g   = [sin.(x3) .* 0.5 * sin.(x2) .* x1.^2 for x3 in z, x1 in y, x2 in x]

    Dz  = CenteredDiff{1}(1, ord, hz, length(z))
    Dx  = CenteredDiff{3}(1, ord, hx, length(x))

    dxzg  = Dx * (Dz * g)

    @test Dz(Dx, g, 100,2,1)     ≈ dxzg[100,2,1]
    @test Dz(Dx, g, 1,1,1)       ≈ dxzg[1,1,1]
    @test Dz(Dx, g, 300,16,600)  ≈ dxzg[300,16,600]
    @test Dz(Dx, g, 10,8,150)    ≈ dxzg[10,8,150]
end

# TG: if I do not add the line below it seems that the Fornberg
# func. is not loaded with Jecco for me locally in the FD_new branch

include("../src/fornberg.jl")

@testset "Fornberg FD weights tests:" begin

    # all stencils provide derivatives at point labelled with 0

    # centered 1st order derivatives i.e. ...,-1,0,1,...
    @test calculate_weights(1,2,1) == [-1//2, 0//1, 1//2]
    @test calculate_weights(1,4,2) == [1//12, -2//3, 0//1, 2//3, -1//12]
    @test calculate_weights(1,6,3) == [-1//60, 3//20, -3//4, 0//1, 3//4, -3//20, 1//60]
    @test calculate_weights(1,8,4) == [1//280, -4//105, 1//5, -4//5, 0//1,
                                       4//5, -1//5, 4//105, -1//280 ]

    # forward 1st order derivatives i.e. 0,1,...
    @test calculate_weights(1,1,0) == [-1//1, 1//1]
    @test calculate_weights(1,2,0) == [-3//2, 2//1, -1//2]
    @test calculate_weights(1,3,0) == [-11//6, 3//1, -3//2, 1//3]
    @test calculate_weights(1,4,0) == [-25//12, 4//1, -3//1, 4//3, -1//4]
    @test calculate_weights(1,5,0) == [-137//60, 5//1, -5//1, 10//3, -5//4, 1//5]
    @test calculate_weights(1,6,0) == [-49//20, 6//1, -15//2, 20//3, -15//4, 6//5, -1//6]

    # forward 1st order derivatives i.e. -1,0,1,2,...
    @test calculate_weights(1,3,1) == [-1//3, -1//2, 1//1, -1//6]
    @test calculate_weights(1,4,1) == [-1//4, -5//6, 3//2, -1//2, 1//12]
    @test calculate_weights(1,5,1) == [-1//5, -13//12, 2//1, -1//1, 1//3, -1//20]
    @test calculate_weights(1,6,1) == [-1//6, -77//60, 5//2, -5//3, 5//6, -1//4, 1//30]
    @test calculate_weights(1,7,1) == [-1//7, -29//20, 3//1, -5//2, 5//3, -3//4, 1//5, -1//42]

    # forward 1st order derivatives i.e. -2,-1,0,1,2,3,...
    @test calculate_weights(1,5,2) == [1//20, -1//2, -1//3, 1//1, -1//4, 1//30]
    @test calculate_weights(1,6,2) == [1//30, -2//5, -7//12, 4//3, -1//2, 2//15, -1//60]

    # backward 1st order derivatives i.e. ...,-1,0
    @test calculate_weights(1,1,1) == [-1//1, 1//1]
    @test calculate_weights(1,2,2) == [1//2, -2//1, 3//2]
    @test calculate_weights(1,3,3) == [-1//3, 3//2, -3//1, 11//6]
    @test calculate_weights(1,4,4) == [1//4, -4//3, 3//1, -4//1, 25//12]

    # backward 1st order derivatives i.e. ...,-2,-1,0,1
    @test calculate_weights(1,3,2) == [1//6, -1//1, 1//2, 1//3]
    @test calculate_weights(1,4,3) == [-1//12, 1//2, -3//2, 5//6, 1//4]
    @test calculate_weights(1,5,4) == [1//20, -1//3, 1//1, -2//1, 13//12, 1//5]

    # backward 1st order derivatives i.e. ...,-3,-2,-1,0,1,2
    @test calculate_weights(1,5,3) == [-1//30, 1//4, -1//1, 1//3, 1//2, -1//20]
    @test calculate_weights(1,6,4) == [1//60, -2//15, 1//2, -4//3, 7//12, 2//5, -1//30]

    # centered 2nd order derivatives i.e. ...,-1,0,1,...
    @test calculate_weights(2,2,1) == [1//1, -2//1, 1//1]
    @test calculate_weights(2,4,2) == [-1//12, 4//3, -5//2, 4//3, -1//12]
    @test calculate_weights(2,6,3) == [1//90, -3//20, 3//2, -49//18,
                                       3//2, -3//20, 1//90]
    @test calculate_weights(2,8,4) == [-1//560, 8//315, -1//5, 8//5, -205//72,
                                       8//5, -1//5, 8//315, -1//560]

    # forward 2nd order derivatives i.e. 0,1,...
    @test calculate_weights(2,2,0) == [1//1, -2//1, 1//1]
    @test calculate_weights(2,3,0) == [2//1, -5//1, 4//1, -1//1]
    @test calculate_weights(2,4,0) == [35//12, -26//3, 19//2, -14//3, 11//12]
    @test calculate_weights(2,5,0) == [15//4, -77//6, 107//6, -13//1, 61//12, -5//6]
    @test calculate_weights(2,6,0) == [203//45, -87//5, 117//4, -254//9, 33//2, -27//5, 137//180]
    @test calculate_weights(2,7,0) == [469//90, -223//10, 879//20, -949//18, 41//1, -201//10,
                                       1019//180, -7//10]

    # forward 2nd order derivatives i.e. -1,0,1,2,...
    @test calculate_weights(2,3,1) == [1//1, -2//1, 1//1, 0//1]
    @test calculate_weights(2,4,1) == [11//12, -5//3, 1//2, 1//3, -1//12]
    @test calculate_weights(2,5,1) == [5//6, -5//4, -1//3, 7//6, -1//2, 1//12]
    @test calculate_weights(2,6,1) == [137//180, -49//60, -17//12, 47//18, -19//12,
                                       31//60, -13//180]

    # forward 2nd order derivatives i.e. -2,-1,0,1,2,3,...
    @test calculate_weights(2,5,2) == [-1//12, 4//3, -5//2, 4//3, -1//12, 0//1]
    @test calculate_weights(2,6,2) == [-13//180, 19//15, -7//3, 10//9,
                                       1//12, -1//15, 1//90]

    # backward 2nd order derivatives i.e. ...,-1,0
    @test calculate_weights(2,2,2) == [1//1, -2//1, 1//1]
    @test calculate_weights(2,3,3) == [-1//1, 4//1, -5//1, 2//1]
    @test calculate_weights(2,4,4) == [11//12, -14//3, 19//2, -26//3, 35//12]

    # backward 2nd order derivatives i.e. ...,-2,-1,0,1
    @test calculate_weights(2,3,2) == [0//1, 1//1, -2//1, 1//1]
    @test calculate_weights(2,4,3) == [-1//12, 1//3, 1//2, -5//3, 11//12]
    @test calculate_weights(2,5,4) == [1//12, -1//2, 7//6, -1//3, -5//4, 5//6]

    # backward 2nd order derivatives i.e. ...,-3,-2,-1,0,1,2
    @test calculate_weights(2,5,3) == [0//1, -1//12, 4//3, -5//2, 4//3, -1//12]
    @test calculate_weights(2,6,4) == [1//90, -1//15, 1//12, 10//9, -7//3, 19//15, -13//180]

    # centered 3rd order derivatives i.e ...,-3,-2,-1,0,1,2,3
    @test calculate_weights(3,6,3) == [1//8, -1//1, 13//8, 0//1, -13//8, 1//1, -1//8]

    # forward 3rd order derivatives i.e 0,1,2,3,...
    @test calculate_weights(3,3,0) == [-1//1, 3//1, -3//1, 1//1]
    @test calculate_weights(3,4,0) == [-5//2, 9//1, -12//1, 7//1, -3//2]

    # forward 3rd order derivatives i.e -1,0,1,2,...
    @test calculate_weights(3,3,1) == [-1//1, 3//1, -3//1, 1//1]
    @test calculate_weights(3,4,1) == [-3//2, 5//1, -6//1, 3//1, -1//2]

    # backward 3rd order derivatives i.e -3,-2,-1,0,...
    @test calculate_weights(3,3,3) == [-1//1, 3//1, -3//1, 1//1]
    @test calculate_weights(3,4,4) == [3//2, -7//1, 12//1, -9//1, 5//2]

    # backward 3rd order derivatives i.e -2,-1,0,1...
    @test calculate_weights(3,3,2) == [-1//1, 3//1, -3//1, 1//1]
    @test calculate_weights(3,4,3) == [1//2, -3//1, 6//1, -5//1, 3//2]

end

# TODO
@testset "EqualSizeStencilFD tests:" begin

    Dx = EqualSizeStencilFD{1}(1, 4, 1.0, 10)
    @test 12*Dx[1,:]  ≈ [-25.0, 48.0, -36.0, 16.0, -3, 0.0, 0.0, 0.0, 0.0, 0.0] atol=1e-15
    @test 12*Dx[2,:]  ≈ [-3.0, -10.0, 18.0, -6.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0] atol=1e-15
    @test 12*Dx[3,:]  ≈ [1.0, -8.0, 0.0, 8.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0] atol=1e-15
    @test 12*Dx[4,:]  ≈ [0.0, 1.0, -8.0, 0.0, 8.0, -1.0, 0.0, 0.0, 0.0, 0.0] atol=1e-15
    @test 12*Dx[5,:]  ≈ [0.0, 0.0, 1.0, -8.0, 0.0, 8.0, -1.0, 0.0, 0.0, 0.0] atol=1e-15
    @test 12*Dx[6,:]  ≈ [0.0, 0.0, 0.0, 1.0, -8.0, 0.0, 8.0, -1.0, 0.0, 0.0] atol=1e-15
    @test 12*Dx[7,:]  ≈ [0.0, 0.0, 0.0, 0.0, 1.0, -8.0, 0.0, 8.0, -1.0, 0.0] atol=1e-15
    @test 12*Dx[8,:]  ≈ [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -8.0, 0.0, 8.0, -1.0] atol=1e-15
    @test 12*Dx[9,:]  ≈ [0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 6.0, -18.0, 10.0, 3.0] atol=1e-15
    @test 12*Dx[10,:] ≈ [0.0, 0.0, 0.0, 0.0, 0.0, 3.0, -16.0, 36.0, -48.0, 25.0] atol=1e-15

end



@testset "Periodic Derivative tests:" begin

    Dx = CenteredDiff{1}(1, 4, 1.0, 10)
    @test 12*Dx[1,:]  ≈ [0.0, 8.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -8.0] atol=1e-15
    @test 12*Dx[2,:]  ≈ [-8.0, 0.0, 8.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0] atol=1e-15
    @test 12*Dx[3,:]  ≈ [1.0, -8.0, 0.0, 8.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0] atol=1e-15
    @test 12*Dx[4,:]  ≈ [0.0, 1.0, -8.0, 0.0, 8.0, -1.0, 0.0, 0.0, 0.0, 0.0] atol=1e-15
    @test 12*Dx[5,:]  ≈ [0.0, 0.0, 1.0, -8.0, 0.0, 8.0, -1.0, 0.0, 0.0, 0.0] atol=1e-15
    @test 12*Dx[6,:]  ≈ [0.0, 0.0, 0.0, 1.0, -8.0, 0.0, 8.0, -1.0, 0.0, 0.0] atol=1e-15
    @test 12*Dx[7,:]  ≈ [0.0, 0.0, 0.0, 0.0, 1.0, -8.0, 0.0, 8.0, -1.0, 0.0] atol=1e-15
    @test 12*Dx[8,:]  ≈ [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -8.0, 0.0, 8.0, -1.0] atol=1e-15
    @test 12*Dx[9,:]  ≈ [-1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -8.0, 0.0, 8.0] atol=1e-15
    @test 12*Dx[10,:] ≈ [8.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -8.0, 0.0] atol=1e-15

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

@testset "EqualSizeStencilFD Derivative tests:" begin

    # 1D case

    xmin   = -2.0*pi
    xmax   =  2.0*pi
    xnodes =  600
    ord    =  6

    hx     = (xmax - xmin) / xnodes

    x  = collect(xmin:hx:xmax-hx)
    f  = sin.(x)

    D1 = EqualSizeStencilFD(1, ord, hx, length(x))
    D2 = EqualSizeStencilFD(2, ord, hx, length(x))

    df  = D1 * f
    d2f = D2 * f

    @test df  ≈ cos.(x)
    @test d2f ≈ -f

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

    Dx     = EqualSizeStencilFD{1}(1, ord, hx, length(x))
    Dz     = EqualSizeStencilFD{3}(1, ord, hz, length(z))

    Dxx    = EqualSizeStencilFD{1}(2, ord, hx, length(x))
    Dzz    = EqualSizeStencilFD{3}(2, ord, hz, length(z))

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

@testset "PeriodicFD cross derivative for general arrays tests:" begin
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


@testset "EqualSizeStencilFD cross derivative for general arrays tests:" begin
    xmin   = -2.0*pi
    xmax   =  2.0*pi
    xnodes =  600
    ord    =  6

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

    Dx  = EqualSizeStencilFD{1}(1, ord, hx, length(x))
    Dz  = EqualSizeStencilFD{3}(1, ord, hz, length(z))

    dxzf  = Dx * (Dz * f)

    @test Dx(Dz, f, 2,16,1)     ≈ dxzf[2,16,1]
    @test Dx(Dz, f, 1,12,2)     ≈ dxzf[1,12,2]
    @test Dx(Dz, f, 100,12,300) ≈ dxzf[100,12,300]
    @test Dx(Dz, f, 600,8,100)  ≈ dxzf[600,8,100]


    g   = [sin.(x3) .* 0.5 * sin.(x2) .* x1.^2 for x3 in z, x1 in y, x2 in x]

    Dz  = EqualSizeStencilFD{1}(1, ord, hz, length(z))
    Dx  = EqualSizeStencilFD{3}(1, ord, hx, length(x))

    dxzg  = Dx * (Dz * g)

    @test Dz(Dx, g, 100,2,1)     ≈ dxzg[100,2,1]
    @test Dz(Dx, g, 1,1,1)       ≈ dxzg[1,1,1]
    @test Dz(Dx, g, 300,16,600)  ≈ dxzg[300,16,600]
    @test Dz(Dx, g, 10,8,150)    ≈ dxzg[10,8,150]
end

# the following test fails; it seems that pointwise it cannot find the
# deriv if one operator is periodic and the other is not.
# TG: Below is the error I get locally:

"""
ERROR: MethodError: no method matching _D_intA_lowB(::Jecco.PeriodicFD{Float64,1,Array{Float64,1}}, ::Jecco.FiniteDiffDeriv{Float64,3,Array{Float64,1},Array{Array{Float64,1},1},Array{Array{Float64,1},1}}, ::Array{Float64,3}, ::Tuple{Int64,Int64,Int64})
Closest candidates are:
  _D_intA_lowB(::Jecco.PeriodicFD{T,N1,S} where S, ::Jecco.PeriodicFD{T,N2,S} where S, ::AbstractArray, ::Any) where {T<:Real, N1, N2} at /home/thanasis/repos/Jecco.jl/src/deriv_operators.jl:469
  _D_intA_lowB(::Jecco.FiniteDiffDeriv{T,N1,S1,S2,S3} where S3 where S2 where S1, ::Jecco.FiniteDiffDeriv{T,N2,S1,S2,S3} where S3 where S2 where S1, ::AbstractArray, ::Any) where {T<:Real, N1, N2} at /home/thanasis/repos/Jecco.jl/src/deriv_operators.jl:709
Stacktrace:
 [1] (::Jecco.PeriodicFD{Float64,1,Array{Float64,1}})(::Jecco.FiniteDiffDeriv{Float64,3,Array{Float64,1},Array{Array{Float64,1},1},Array{Array{Float64,1},1}}, ::Array{Float64,3}, ::Int64, ::Int64, ::Int64) at /home/thanasis/repos/Jecco.jl/src/deriv_operators.jl:409
 [2] top-level scope at REPL[51]:1
"""

@testset "PeriodicFD and EqualSizeStencilFD cross derivative for general arrays tests:" begin
    xmin   = -2.0*pi
    xmax   =  2.0*pi
    xnodes =  600
    ord    =  6

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

    Dx  = EqualSizeStencilFD{1}(1, ord, hx, length(x))
    Dz  = CenteredDiff{3}(1, ord, hz, length(z))

    dxzf  = Dx * (Dz * f)

    @test Dx(Dz, f, 2,16,1)     ≈ dxzf[2,16,1]
    @test Dx(Dz, f, 1,12,2)     ≈ dxzf[1,12,2]
    @test Dx(Dz, f, 100,12,300) ≈ dxzf[100,12,300]
    @test Dx(Dz, f, 600,8,100)  ≈ dxzf[600,8,100]


    g   = [sin.(x3) .* 0.5 * sin.(x2) .* x1.^2 for x3 in z, x1 in y, x2 in x]

    Dz  = CenteredDiff{1}(1, ord, hz, length(z))
    Dx  = EqualSizeStencilFD{3}(1, ord, hx, length(x))

    dxzg  = Dx * (Dz * g)

    @test Dz(Dx, g, 100,2,1)     ≈ dxzg[100,2,1]
    @test Dz(Dx, g, 1,1,1)       ≈ dxzg[1,1,1]
    @test Dz(Dx, g, 300,16,600)  ≈ dxzg[300,16,600]
    @test Dz(Dx, g, 10,8,150)    ≈ dxzg[10,8,150]
end


ABCS  = zeros(4)
AA    = zeros(2,2)
BB    = zeros(2,2)
CC    = zeros(2,2)
SS    = zeros(2)

svars    = tuple(ones(11)...)
fvars    = tuple(ones(37)...)
sdvars   = tuple([ZeroPotential(); ones(70)]...)
bdgvars  = tuple(ones(71)...)
phidvars = tuple([ZeroPotential(); ones(72)]...)
avars    = tuple([ZeroPotential(); ones(76)]...)
xivars   = tuple(ones(84)...)


@testset "S equation outer grid coefficients:" begin
    AdS5_3_1.S_eq_coeff!(ABCS, svars, AdS5_3_1.Outer())

    @test all(ABCS .≈ [6.0, 12.0, 10.381097845541817, 0.0])
end

@testset "Fxy equations outer grid coefficients:" begin
    AdS5_3_1.Fxy_eq_coeff!(AA, BB, CC, SS, fvars, AdS5_3_1.Outer())

    @test all( AA .≈ [5.43656365691809 0.0; 0.0 2.0] )
    @test all( BB .≈ [-12.9449900106386 5.626860407847019; -4.4222650840900215 4.762195691083631] )
    @test all( CC .≈ [149.68137141822632 -24.40497260570832; 42.56474178480922 2.459692111055068] )
    @test all( SS .≈ [0.0, 0.0] )
end

@testset "Sd equation outer grid coefficients:" begin
    AdS5_3_1.Sd_eq_coeff!(ABCS, sdvars, AdS5_3_1.Outer())

    @test all( ABCS .≈ [0.0, -32.61938194150854, 65.23876388301709,1037.7259866620095])
end

@testset "B2d equation outer grid coefficients:" begin
    AdS5_3_1.B2d_eq_coeff!(ABCS, bdgvars, AdS5_3_1.Outer())

    @test all( ABCS .≈ [0.0, -32.61938194150854, 48.92907291226281, 31.353708676211966] )
end

@testset "B1dGd equations outer grid coefficients:" begin
    AdS5_3_1.B1dGd_eq_coeff!(AA, BB, CC, SS, bdgvars, AdS5_3_1.Outer())

    @test all( AA .≈ [0.0 0.0; 0.0 0.0] )
    @test all( BB .≈ [-32.61938194150854 0.0; 0.0 -32.61938194150854] )
    @test all( CC .≈ [73.77180356980473 24.842730657541917; -59.15297244604868 48.92907291226281] )
    @test all( SS .≈ [-322.48341532407966, -82.91438383027555] )
end

@testset "phid equations outer grid coefficients:" begin
    AdS5_3_1.phid_eq_coeff!(ABCS, phidvars, AdS5_3_1.Outer())

    @test all( ABCS .≈ [0.0, -21.74625462767236, 32.61938194150854, 473.5501053413936] )
end

@testset "A equations outer grid coefficients:" begin
    AdS5_3_1.A_eq_coeff!(ABCS, avars, AdS5_3_1.Outer())

    @test all( ABCS .≈ [16.30969097075427, 32.61938194150854, 0.0,-2722.8310332986166 ])
end

@testset "xi_t equation coefficients:" begin
    abccS = AdS5_3_1.xi_t_eq_coeff(xivars, AdS5_3_1.Outer())

    @test all( abccS .≈ (-11.401909375823356, -1.5430806348152437, 6.3890560989306495, -55.88040794752505, 15.338669886560382, -132.77982902946636, -366.5278066874655) )
end

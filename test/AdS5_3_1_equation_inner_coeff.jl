
using Test
using BenchmarkTools

using Jecco
using Jecco.AdS5_3_1

ABCS  = zeros(4)
AA    = zeros(2,2)
BB    = zeros(2,2)
CC    = zeros(2,2)
SS    = zeros(2)

svars  = Jecco.AdS5_3_1.SVars(ones(11)...)
fvars  = Jecco.AdS5_3_1.FVars(ones(37)...)
sdvars = Jecco.AdS5_3_1.SdVars(ZeroPotential(), ones(70)...)


@testset "S equation inner grid coefficients:" begin
    Jecco.AdS5_3_1.S_eq_coeff!(ABCS, svars, Inner())

    @test all(ABCS .≈ [6.0, 48.0, 133.42988060987634, 130.85976121975267])
end
@btime Jecco.AdS5_3_1.S_eq_coeff!($ABCS, $svars, Inner())

@testset "Fxy equations inner grid coefficients:" begin
    Jecco.AdS5_3_1.Fxy_eq_coeff!(AA, BB, CC, SS, fvars, Inner())

    @test all( AA .≈ [48.92907291226281 0.0; 0.0 18.0] )
    @test all( BB .≈ [811.6226411252798 -151.92523101186953; 119.40115727043062 41.42071634074196] )
    @test all( CC .≈ [9249.882276773546 -1318.4601133442727; 2930.339069731014 256.1897346300231] )
    @test all( SS .≈ [9577.64621761626, 3981.7746343640238])
end
@btime Jecco.AdS5_3_1.Fxy_eq_coeff!($AA, $BB, $CC, $SS, $fvars, Inner())

@testset "Sd equation inner grid coefficients:" begin
    Jecco.AdS5_3_1.Sd_eq_coeff!(ABCS, sdvars, Inner())

    @test all( ABCS .≈ [0.0, -880.7233124207306, -2544.311791437666, 279203.2839670398] )
end
@btime Jecco.AdS5_3_1.Sd_eq_coeff!($ABCS, $sdvars, Inner())

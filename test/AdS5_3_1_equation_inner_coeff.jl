
using Test
using BenchmarkTools

using Jecco
using Jecco.AdS5_3_1

ABCS  = zeros(4)

svars = Jecco.AdS5_3_1.SVars(ones(11)...)


@testset "S equation inner grid coefficients:" begin
    Jecco.AdS5_3_1.S_eq_coeff!(ABCS, svars, Inner())
    @test all(ABCS .≈ [6.0, 48.0, 133.42988060987634, 130.85976121975267])
end

@btime Jecco.AdS5_3_1.S_eq_coeff!($ABCS, $svars, Inner())

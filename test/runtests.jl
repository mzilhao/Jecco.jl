
using Jecco
using Test

@time @testset "derivative operators tests:" begin include("deriv_operators_tests.jl") end

@time @testset "grid tests:" begin include("grid_tests.jl") end

@time @testset "output tests:" begin include("output_tests.jl") end

# AdS5_3_1 inner grid coefficients tests:
include("AdS5_3_1_equation_inner_coeff.jl")


nothing

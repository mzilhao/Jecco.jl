
using Jecco
using Test

@time @testset "derivative operators tests:" begin include("deriv_operators_tests.jl") end

@time @testset "grid tests:" begin include("grid_tests.jl") end

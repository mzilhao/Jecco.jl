module Jecco

using LinearAlgebra
using Vivi
using Parameters

include("deriv_operators.jl")
include("output.jl")
include("spectral.jl")
include("KG_3_1/KG_3_1.jl")

export CenteredDiff, ChebDeriv

end # module

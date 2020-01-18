module Jecco

using LinearAlgebra
# using Vivi
using Parameters

include("deriv_operators.jl")
include("grid.jl")
include("output.jl")
include("spectral.jl")
include("KG_3_1/KG_3_1.jl")

export CartCoord, SpectralCoord, Grid, xx
export CenteredDiff, ChebDeriv

end # module

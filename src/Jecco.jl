module Jecco

using HDF5
using LinearAlgebra
using FFTW

include("types.jl")
include("fornberg.jl")
include("deriv_operators.jl")
include("grid.jl")
include("input.jl")
include("output.jl")
include("startup.jl")
include("spectral.jl")
include("filtering.jl")
include("utils.jl")

export AbstractPartition, FlattenedVector

export AbstractCoord, CartesianCoord, GaussLobattoCoord
export Chart, Atlas
export Cartesian, GaussLobatto

export CenteredDiff, EqualSizeStencilFD, ChebDeriv, FourierDeriv
export ChebInterpolator
export KO_Filter, Exp_Filter
export KO_Centered

export OpenPMDTimeSeries, get_field
export FieldTimeSeries, get_coords

include("KG_3_1/KG_3_1.jl")
include("AdS5_3_1/AdS5_3_1.jl")
include("AdS4_3_1/AdS4_3_1.jl")

end # module

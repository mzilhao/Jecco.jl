module AdS5_3_1

using Jecco
using Parameters
using RecursiveArrayTools

abstract type GridType end
struct Inner <: GridType end
struct Outer <: GridType end

export ParamBase, Grid3D, ParamID #, ParamEvol, ParamIO
export Potential, ZeroPotential
export BlackBrane, init!, init
export Inner, Outer, System
export EvolVars
export BulkVars, BoundaryVars, GaugeVars, BaseVars


include("types.jl")

# include("param.jl")
include("system.jl")
include("initial_data.jl")
include("potential.jl")
# include("dphidt.jl")
include("equation_inner_coeff.jl")
include("equation_outer_coeff.jl")
include("inner_to_outer.jl")
include("solve_nested.jl")
# include("rhs.jl")
# include("run.jl")
# include("ibvp.jl")

end

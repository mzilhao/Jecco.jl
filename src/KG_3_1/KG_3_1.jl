module KG_3_1

using Jecco
using LinearAlgebra

import Base.Threads.@threads
import Base.Threads.@spawn

# abstract types and structs used throughout
include("types.jl")

include("system.jl")

# include("potential.jl")

# include("initial_data.jl")


# include("param.jl")

# include("dphidt.jl")
# include("equation_coeff.jl")
# include("solve_nested.jl")
# include("rhs.jl")
# include("run.jl")
# include("ibvp.jl")

# export ParamBase, ParamGrid, ParamID, ParamEvol, ParamIO
# export Potential
# export VV # this will contain the potential
# export System
# export BulkVars, BoundaryVars, AllVars

export SpecCartGrid3D
export System, SystemPartition

end

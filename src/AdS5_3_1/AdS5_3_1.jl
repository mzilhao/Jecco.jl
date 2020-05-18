module AdS5_3_1

using Jecco
using Parameters
using RecursiveArrayTools

# abstract types and structs used throughout
include("types.jl")

include("system.jl")

include("potential.jl")

include("initial_data.jl")

# nested system
include("equation_inner_coeff.jl")
include("equation_outer_coeff.jl")
include("inner_to_outer.jl")
include("solve_nested.jl")


export Grid3D
export Potential, ZeroPotential
export BlackBrane, init_data!, init_data
export Inner, Outer, System
export BulkEvol, Boundary, Gauge, Bulk
export SystemPartition
export BulkVars, BoundaryVars, GaugeVars, BaseVars

end

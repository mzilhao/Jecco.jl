module AdS5_3_1

using Jecco
using OrdinaryDiffEq

# abstract types and structs used throughout
include("types.jl")

include("system.jl")

include("potential.jl")

include("initial_data.jl")

# time marching orders
include("setup_rhs.jl")

# nested system
include("equation_inner_coeff.jl")
include("equation_outer_coeff.jl")
include("inner_to_outer.jl")
include("solve_nested.jl")

# evolution equations
include("compute_boundary_t.jl")
include("compute_gauge_t.jl")
include("compute_bulkevolved_t.jl")

# run the model
include("run.jl")

export SpecCartGrid3D
export Potential, ZeroPotential
export BlackBrane, BlackBraneB1Pert, init_data!, init_data
export ConstantAH
export AffineNull, EvolTest0
export BulkEvolved, BulkConstrained, Boundary, Gauge, Bulk
export System, SystemPartition
export BulkEvolvedPartition, BulkConstrainedPartition, BulkDerivPartition
export Nested
export Integration, InOut, run_model

end

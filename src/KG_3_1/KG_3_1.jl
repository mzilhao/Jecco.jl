module KG_3_1

using Jecco
using LinearAlgebra
using OrdinaryDiffEq

import Base.Threads.@threads
import Base.Threads.@spawn

# abstract types and structs used throughout
include("types.jl")

include("system.jl")

include("potential.jl")

include("initial_data.jl")

# nested system
include("equation_coeff.jl")
include("solve_nested.jl")

# evolution equations
include("compute_bulkevolved_t.jl")

# time marching orders
include("setup_rhs.jl")

# input/output
include("IO.jl")

# run the model
include("run.jl")

# always set the number of BLAS threads to 1 upon loading the module. by default
# it uses a bunch of them and we don't want that since they trample over each
# other when solving the nested systems equations. it's much better to thread
# over the loop. see also the discussion here:
# https://github.com/JuliaLang/julia/issues/33409
#
# this saves us the tedious task of always setting OMP_NUM_THREADS=1 before
# launching julia.
function __init__()
    LinearAlgebra.BLAS.set_num_threads(1)
    nothing
end

export SpecCartGrid3D
export Potential, ConstPotential, SquarePotential
export System, SystemPartition
export BulkEvolvedPartition, BulkConstrainedPartition, Boundary
export Uniform2D, Sine2D
export AffineNull
export Integration, InOut, run_model

end

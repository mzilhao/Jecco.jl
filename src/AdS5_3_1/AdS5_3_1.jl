module AdS5_3_1

using Jecco
using Jecco: mul_col!
using SparseArrays: SparseMatrixCSC
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

# finding the Apparent Horizon
include("equations_AH.jl")
include("find_AH.jl")

# input/output
include("IO.jl")

# run the model
include("run.jl")

# expressions for the Vacuum Expectation Values
include("VEVs.jl")

# random utilities
include("utils.jl")

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
export Potential, ZeroPotential, Phi8Potential
export BlackBrane, BlackBranePert
export PhiGaussian_u
export init_data!, init_data, restore!
export ConstantAH, AHF
export AffineNull, EvolTest0
export BulkEvolved, BulkConstrained, Boundary, Gauge, Bulk
export System, SystemPartition
export BulkEvolvedPartition, BulkConstrainedPartition, BulkDerivPartition
export Nested
export Integration, InOut, run_model
export get_energy, get_Jx, get_Jy, get_px, get_py, get_pz, get_pxy, get_Ophi
export BoundaryTimeSeries, XiTimeSeries, BulkTimeSeries, VEVTimeSeries
export convert_to_mathematica

end

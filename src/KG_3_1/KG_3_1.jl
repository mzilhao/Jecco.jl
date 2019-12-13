module KG_3_1

using Jecco
using Vivi
using Parameters

export ParamBase, ParamGrid, ParamID, ParamEvol, ParamIO
export Potential
export VV # this will contain the potential
export System
export BulkVars, BoundaryVars, AllVars

struct BulkVars{A}
    phi    :: A
    S      :: A
    Sd     :: A
    phid   :: A
    A      :: A
    dphidt :: A
end
BulkVars(phi, S, Sd, phid, A, dphidt) =  BulkVars{typeof(phi)}(phi, S, Sd, phid, A, dphidt)
function BulkVars(phi::Array{<:Number,N}) where {N}
    S      = similar(phi)
    Sd     = similar(phi)
    phid   = similar(phi)
    A      = similar(phi)
    dphidt = similar(phi)
    BulkVars{typeof(phi)}(phi, S, Sd, phid, A, dphidt)
end

BulkVars(phis::Vector) = [BulkVars(phi) for phi in phis]

function Base.getindex(bulk::BulkVars, i::Int)
    phi    = bulk.phi[i]
    S      = bulk.S[i]
    Sd     = bulk.Sd[i]
    phid   = bulk.phid[i]
    A      = bulk.A[i]
    dphidt = bulk.dphidt[i]
    BulkVars{typeof(phi)}(phi, S, Sd, phid, A, dphidt)
end

function Base.getindex(bulk::BulkVars, kr::AbstractRange)
    phi    = bulk.phi[kr]
    S      = bulk.S[kr]
    Sd     = bulk.Sd[kr]
    phid   = bulk.phid[kr]
    A      = bulk.A[kr]
    dphidt = bulk.dphidt[kr]
    BulkVars{typeof(phi)}(phi, S, Sd, phid, A, dphidt)
end

function Base.getindex(bulk::BulkVars, I::Vararg)
    phi    = bulk.phi[I...]
    S      = bulk.S[I...]
    Sd     = bulk.Sd[I...]
    phid   = bulk.phid[I...]
    A      = bulk.A[I...]
    dphidt = bulk.dphidt[I...]
    BulkVars{typeof(phi)}(phi, S, Sd, phid, A, dphidt)
end

function Base.getindex(bulk::BulkVars, ::Colon)
    phi    = bulk.phi[:]
    S      = bulk.S[:]
    Sd     = bulk.Sd[:]
    phid   = bulk.phid[:]
    A      = bulk.A[:]
    dphidt = bulk.dphidt[:]
    BulkVars{typeof(phi)}(phi, S, Sd, phid, A, dphidt)
end

Base.lastindex(bulk::BulkVars) = lastindex(bulk.phi)
Base.lastindex(bulk::BulkVars, i::Int) = lastindex(bulk.phi, i)



struct BoundaryVars{A}
    a4   :: A
end

mutable struct AllVars{T}
    u        :: T

    phi_d0   :: T
    phi_du   :: T
    phi_dxx  :: T
    phi_dyy  :: T

    Sd_d0    :: T

    phid_d0  :: T
    phid_du  :: T

    A_d0     :: T
end
function AllVars{T}() where {T<:AbstractFloat}
    N = 1 + 4 + 1 + 2 + 1
    array = zeros(N)
    AllVars{T}(array...)
end

include("param.jl")
include("system.jl")
include("initial_data.jl")
include("potential.jl")
include("dphidt.jl")
include("equation_coeff.jl")
include("solve_nested.jl")
include("rhs.jl")
include("ibvp.jl")

par_base = ParamBase(
    which_potential = "square",
)

# define potential. TODO: move this somewhere else
VV = Potential(par_base)

end

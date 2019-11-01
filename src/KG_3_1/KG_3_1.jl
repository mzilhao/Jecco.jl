module KG_3_1

using Jecco
using Vivi
using Parameters

export Param
export System
export BulkVars, BoundaryVars, AllVars

@with_kw struct Param
    A0x         :: Float64
    A0y         :: Float64

    tmax        :: Float64
    out_every   :: Int

    xmin        :: Float64
    xmax        :: Float64
    xnodes      :: Int
    ymin        :: Float64
    ymax        :: Float64
    ynodes      :: Int
    umin        :: Float64
    umax        :: Float64
    unodes      :: Int

    # dtfac       :: Float64    = 0.5
    dt          :: Float64

    folder      :: String  = "./data"
    prefix      :: String  = "phi"
end

struct System{C,D,E} <: Vivi.System
    coords :: C
    uderiv :: D
    xderiv :: E
    yderiv :: E
end

function System(coords::CoordSystem)
    # FIXME
    ord    = 4
    BC     = :periodic

    # dx     = xcoord.delta
    # dy     = ycoord.delta
    # dt0    = p.dtfac * min(dx, dy)

    derivs = Vivi.Deriv(coords, (nothing, ord, ord), (nothing, BC, BC))
    uderiv = derivs[1]
    xderiv = derivs[2]
    yderiv = derivs[3]

    System{typeof(coords), typeof(uderiv), typeof(xderiv)}(coords, uderiv, xderiv, yderiv)
end

function System(ucoord::SpectralCoord, xcoord::CartCoord, ycoord::CartCoord)
    coords = Vivi.CoordSystem{Float64}("uxy", [ucoord, xcoord, ycoord])
    System(coords)
end


# TODO: determine it using the metric
function timestep(sys::System, f)
end


struct BulkVars{A}
    phi  :: A
    S    :: A
    Sd   :: A
    phid :: A
    A    :: A
end
BulkVars(phi, S, Sd, phid, A) =  BulkVars{typeof(phi)}(phi, S, Sd, phid, A)
function BulkVars(phi::Array{<:Number,N}) where {N}
    S    = similar(phi)
    Sd   = similar(phi)
    phid = similar(phi)
    A    = similar(phi)
    BulkVars{typeof(phi)}(phi, S, Sd, phid, A)
end

BulkVars(phis::Vector) = [BulkVars(phi) for phi in phis]

function Base.getindex(bulk::BulkVars, i::Int)
    phi   = bulk.phi[i]
    S     = bulk.S[i]
    Sd    = bulk.Sd[i]
    phid  = bulk.phid[i]
    A     = bulk.A[i]
    BulkVars{typeof(phi)}(phi, S, Sd, phid, A)
end

function Base.getindex(bulk::BulkVars, kr::AbstractRange)
    phi   = bulk.phi[kr]
    S     = bulk.S[kr]
    Sd    = bulk.Sd[kr]
    phid  = bulk.phid[kr]
    A     = bulk.A[kr]
    BulkVars{typeof(phi)}(phi, S, Sd, phid, A)
end

function Base.getindex(bulk::BulkVars, I::Vararg)
    phi   = bulk.phi[I...]
    S     = bulk.S[I...]
    Sd    = bulk.Sd[I...]
    phid  = bulk.phid[I...]
    A     = bulk.A[I...]
    BulkVars{typeof(phi)}(phi, S, Sd, phid, A)
end

function Base.getindex(bulk::BulkVars, ::Colon)
    phi   = bulk.phi[:]
    S     = bulk.S[:]
    Sd    = bulk.Sd[:]
    phid  = bulk.phid[:]
    A     = bulk.A[:]
    BulkVars{typeof(phi)}(phi, S, Sd, phid, A)
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


include("initial_data.jl")
include("dphidt.jl")
include("equation_coeff.jl")
include("solve_nested.jl")
include("rhs.jl")
include("ibvp.jl")

end


"""
Extend this type for different `Potential` choices
"""
abstract type Potential end

"""
Extend this type for different `InitialData` choices
"""
abstract type InitialData end


struct BulkEvolved{T,A,S} <: AbstractPartition{3,A}
    x :: S
end
"""
    BulkEvolved{T}(undef, Nu, Nx, Ny)

Construct a container of uninitialized Arrays to hold the bulk variable
evolved in time: phi
"""
function BulkEvolved{T}(::UndefInitializer, Nu::Int, Nx::Int, Ny::Int) where {T<:Real}
    phi = Array{T}(undef, Nu, Nx, Ny)
    x   = (phi=phi,)
    BulkEvolved{T,eltype(x),typeof(x)}(x)
end

getphi(ff::BulkEvolved) = ff.x.phi


struct BulkConstrained{T,A,S} <: AbstractPartition{3,A}
    x :: S
end
"""
    BulkConstrained{T}(undef, Nu, Nx, Ny)

Construct a container of uninitialized Arrays to hold all the bulk variables
that are constrained (not evolved in time): S, phid, Sd, A
"""
function BulkConstrained{T}(::UndefInitializer, Nu::Int, Nx::Int, Ny::Int) where {T<:Real}
    S    = Array{T}(undef, Nu, Nx, Ny)
    phid = Array{T}(undef, Nu, Nx, Ny)
    Sd   = Array{T}(undef, Nu, Nx, Ny)
    A    = Array{T}(undef, Nu, Nx, Ny)
    x    = (S=S, phid=phid, Sd=Sd, A=A)
    BulkConstrained{T,eltype(x),typeof(x)}(x)
end

getS(ff::BulkConstrained)    = ff.x.S
getphid(ff::BulkConstrained) = ff.x.phid
getSd(ff::BulkConstrained)   = ff.x.Sd
getA(ff::BulkConstrained)    = ff.x.A


struct Boundary{T,A,S} <: AbstractPartition{3,A}
    x :: S
end
"""
    Boundary{T}(undef, Nx, Ny)

Construct a container of uninitialized Arrays to hold the boundary
variables a4.

These variables are automatically defined on a `(1,Nx,Ny)` grid, rather than
a `(Nx,Ny)` one.
"""
function Boundary{T}(::UndefInitializer, Nx::Int, Ny::Int) where {T<:Real}
    a4  = Array{T}(undef, 1, Nx, Ny)
    x   = (a4=a4,)
    Boundary{T,eltype(x),typeof(x)}(x)
end

geta4(ff::Boundary) = ff.x.a4


struct BulkDeriv{T}
    Du_phi  :: Array{T,3}
    Du_S    :: Array{T,3}
    Du_Sd   :: Array{T,3}
    Du_A    :: Array{T,3}
    Duu_phi :: Array{T,3}
    Duu_S   :: Array{T,3}
    Duu_A   :: Array{T,3}
end

function BulkDeriv{T}(::UndefInitializer, Nu::Int, Nx::Int, Ny::Int) where {T<:Real}
    Du_phi   = Array{T}(undef, Nu, Nx, Ny)
    Du_S     = Array{T}(undef, Nu, Nx, Ny)
    Du_Sd    = Array{T}(undef, Nu, Nx, Ny)
    Du_A     = Array{T}(undef, Nu, Nx, Ny)
    Duu_phi  = Array{T}(undef, Nu, Nx, Ny)
    Duu_S    = Array{T}(undef, Nu, Nx, Ny)
    Duu_A    = Array{T}(undef, Nu, Nx, Ny)
    BulkDeriv{T}(Du_phi, Du_S, Du_Sd, Du_A, Duu_phi, Duu_S, Duu_A)
end




"""
    BulkPartition{N,A} <: AbstractPartition{N,A}

Container to store (bulk) quantities that may be spread across different grid
partitions. This is to be thought as a `Tuple` (or an `Array`) of `Bulk` objects.
"""
struct BulkPartition{N,A} <: AbstractPartition{N,A}
    x :: NTuple{N,A}
end
BulkPartition(x...) = BulkPartition(tuple(x...))


function BulkEvolved(bulks::BulkPartition{N}) where{N}
    f = ntuple(i -> BulkEvolved(bulks[i]), N)
    BulkPartition(f)
end

function BulkConstrained(bulks::BulkPartition{N}) where{N}
    f = ntuple(i -> BulkConstrained(bulks[i]), N)
    BulkPartition(f)
end



# ### TODO: remove below


# struct BulkVars{A}
#     phi    :: A
#     S      :: A
#     Sd     :: A
#     phid   :: A
#     A      :: A
#     dphidt :: A
# end
# BulkVars(phi, S, Sd, phid, A, dphidt) =  BulkVars{typeof(phi)}(phi, S, Sd, phid, A, dphidt)
# function BulkVars(phi::Array{<:Number,N}) where {N}
#     S      = similar(phi)
#     Sd     = similar(phi)
#     phid   = similar(phi)
#     A      = similar(phi)
#     dphidt = similar(phi)
#     BulkVars{typeof(phi)}(phi, S, Sd, phid, A, dphidt)
# end

# BulkVars(phis::Vector) = [BulkVars(phi) for phi in phis]

# function Base.getindex(bulk::BulkVars, i::Int)
#     phi    = bulk.phi[i]
#     S      = bulk.S[i]
#     Sd     = bulk.Sd[i]
#     phid   = bulk.phid[i]
#     A      = bulk.A[i]
#     dphidt = bulk.dphidt[i]
#     BulkVars{typeof(phi)}(phi, S, Sd, phid, A, dphidt)
# end

# function Base.getindex(bulk::BulkVars, kr::AbstractRange)
#     phi    = bulk.phi[kr]
#     S      = bulk.S[kr]
#     Sd     = bulk.Sd[kr]
#     phid   = bulk.phid[kr]
#     A      = bulk.A[kr]
#     dphidt = bulk.dphidt[kr]
#     BulkVars{typeof(phi)}(phi, S, Sd, phid, A, dphidt)
# end

# function Base.getindex(bulk::BulkVars, I::Vararg)
#     phi    = bulk.phi[I...]
#     S      = bulk.S[I...]
#     Sd     = bulk.Sd[I...]
#     phid   = bulk.phid[I...]
#     A      = bulk.A[I...]
#     dphidt = bulk.dphidt[I...]
#     BulkVars{typeof(phi)}(phi, S, Sd, phid, A, dphidt)
# end

# function Base.getindex(bulk::BulkVars, ::Colon)
#     phi    = bulk.phi[:]
#     S      = bulk.S[:]
#     Sd     = bulk.Sd[:]
#     phid   = bulk.phid[:]
#     A      = bulk.A[:]
#     dphidt = bulk.dphidt[:]
#     BulkVars{typeof(phi)}(phi, S, Sd, phid, A, dphidt)
# end

# Base.lastindex(bulk::BulkVars) = lastindex(bulk.phi)
# Base.lastindex(bulk::BulkVars, i::Int) = lastindex(bulk.phi, i)

# function setup(par_base)
#     global VV = Potential(par_base)
# end


# struct BoundaryVars{A}
#     a4   :: A
# end

# mutable struct AllVars{T}
#     u        :: T

#     phi_d0   :: T
#     phi_du   :: T
#     phi_dxx  :: T
#     phi_dyy  :: T

#     Sd_d0    :: T

#     phid_d0  :: T
#     phid_du  :: T

#     A_d0     :: T
# end
# function AllVars{T}() where {T<:AbstractFloat}
#     N = 1 + 4 + 1 + 2 + 1
#     array = zeros(N)
#     AllVars{T}(array...)
# end

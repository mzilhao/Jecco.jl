
struct BulkEvolved{T}
    phi :: Array{T,3}
end

struct BulkConstrained{T}
    S    :: Array{T,3}
    phid :: Array{T,3}
    Sd   :: Array{T,3}
    A    :: Array{T,3}
end

struct BulkDeriv{T}
    Du_phi  :: Array{T,3}
    Du_S    :: Array{T,3}
    Du_Sd   :: Array{T,3}
    Du_A    :: Array{T,3}
    Duu_phi :: Array{T,3}
    Duu_S   :: Array{T,3}
    Duu_A   :: Array{T,3}
end

struct Boundary{T}
    a4  :: Array{T,3}
end


"""
    BulkEvolved{T}(undef, Nu, Nx, Ny)

Construct a container of uninitialized Arrays to hold the bulk variable
evolved in time phi
"""
function BulkEvolved{T}(::UndefInitializer, Nu::Int, Nx::Int, Ny::Int) where {T<:Real}
    phi = Array{T}(undef, Nu, Nx, Ny)
    BulkEvolved{T}(phi)
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
    BulkConstrained{T}(S, phid, Sd, A)
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
    Boundary{T}(undef, Nx, Ny)

Construct a container of uninitialized Arrays to hold the boundary
variables a4.

These variables are automatically defined on a `(1,Nx,Ny)` grid, rather than
a `(Nx,Ny)`.
"""
function Boundary{T}(::UndefInitializer, Nx::Int, Ny::Int) where {T<:Real}
    a4  = Array{T}(undef, 1, Nx, Ny)
    Boundary{T}(a4)
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

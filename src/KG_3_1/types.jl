
"""
Extend this type for different `Potential` choices
"""
abstract type Potential end

"""
Extend this type for different `InitialData` choices
"""
abstract type InitialData end


struct BulkEvolved{T,A,S} <: FlattenedVector{T,3,A}
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


struct BulkConstrained{T,A,S} <: FlattenedVector{T,3,A}
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


struct Boundary{T,A,S} <: FlattenedVector{T,3,A}
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


struct EvolVars{T,N,A} <: FlattenedVector{T,N,A}
    x :: NTuple{N,A}
end
EvolVars(x::NTuple{N,A}) where {N,A} = EvolVars{eltype(A),N,A}(x)

Base.similar(ff::EvolVars{T,N,S}) where {T,N,S} = EvolVars{T,N,S}(similar.(ff.x))

"""
    EvolVars(boundary::Boundary, gauge::Gauge, bulkevols::NTuple)

Build a container to store all the evolved quantities as elements of an
`NTuple`. The idea is to treat them as a single column vector for the point of
view of the time evolution routine.
"""
function EvolVars(bulkevols::BulkPartition{Nsys}) where {Nsys}
    EvolVars(bulkevols.x)
end


"""
Extend this type for different `Potential` choices
"""
abstract type Potential end

"""
Extend this type for different `InitialData` choices
"""
abstract type InitialData end

"""
Extend this type for different `EvolutionEquations`.
"""
abstract type EvolutionEquations end

Base.@kwdef struct AffineNull{TP<:Potential} <: EvolutionEquations
    potential      :: TP  = ConstPotential()
end

"""
Parameters for the time evolution
"""
Base.@kwdef struct Integration{T,Tdt,S}
    dt              :: Tdt  = :auto
    tmax            :: T
    ODE_method      :: S    = AB3()
    adaptive        :: Bool = false
    # relative tolerance for adaptive integrators
    reltol          :: Float64 = 1e-6
    # filter_poststep :: Bool = true
end


"""
Parameters for Input/Output
"""
Base.@kwdef struct InOut
    # negative values suppress output
    out_bulk_every              :: Int  = -1

    out_bulk_every_t            :: Float64  = -1.0

    checkpoint_every_walltime_hours :: Float64 = -1.0

    # stop and checkpoint upon reaching this walltime
    max_walltime       :: Float64 = 1.e20

    # trigger termination
    termination_from_file :: Bool    = true
    check_file_every      :: Int     = 10
    termination_file      :: String  = "TERMINATE"

    # name of script
    _parfile           :: String  = splitext(basename(Base.source_path()))[1]

    # use name of script by default
    out_dir            :: String  = _parfile
    checkpoint_dir     :: String  = _parfile

    recover            :: Symbol  = :auto
    recover_dir        :: String  = _parfile

    # be very careful with this option! it will remove the whole folder contents
    # if set to true! use only for fast debugging runs
    remove_existing    :: Bool  = false
end


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

function Base.similar(ff::BulkEvolved{T}) where {T}
    phi = similar(getphi(ff))
    x   = (phi=phi,)
    BulkEvolved{T,eltype(x),typeof(x)}(x)
end
function Base.copy(ff::BulkEvolved{T}) where {T}
    phi = copy(getphi(ff))
    x   = (phi=phi,)
    BulkEvolved{T,eltype(x),typeof(x)}(x)
end
function Base.zero(ff::BulkEvolved{T}) where {T}
    phi = zero(getphi(ff))
    x   = (phi=phi,)
    BulkEvolved{T,eltype(x),typeof(x)}(x)
end


struct BulkConstrained{T,A,S} <: FlattenedVector{T,3,A}
    x :: S
end
"""
    BulkConstrained{T}(undef, Nu, Nx, Ny)

Construct a container of uninitialized Arrays to hold all the bulk variables
that are constrained (not evolved in time): S, phid, Sd, A
"""
function BulkConstrained{T}(::UndefInitializer, Nu::Int, Nx::Int, Ny::Int) where {T<:Real}
    phid = Array{T}(undef, Nu, Nx, Ny)
    Sd   = Array{T}(undef, Nu, Nx, Ny)
    A    = Array{T}(undef, Nu, Nx, Ny)
    x    = (phid=phid, Sd=Sd, A=A)
    BulkConstrained{T,eltype(x),typeof(x)}(x)
end

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
    Du_phid :: Array{T,3}
end

function BulkDeriv{T}(::UndefInitializer, Nu::Int, Nx::Int, Ny::Int) where {T<:Real}
    Du_phi   = Array{T}(undef, Nu, Nx, Ny)
    Du_phid  = Array{T}(undef, Nu, Nx, Ny)
    BulkDeriv{T}(Du_phi, Du_phid)
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

Base.similar(bulks::BulkPartition{N,A}) where {N,A} = BulkPartition{N,A}(similar.(bulks.x))


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
Base.copy(ff::EvolVars{T,N,S}) where {T,N,S} = EvolVars{T,N,S}(copy.(ff.x))
Base.zero(ff::EvolVars{T,N,S}) where {T,N,S} = EvolVars{T,N,S}(zero.(ff.x))


"""
    EvolVars(boundary::Boundary, gauge::Gauge, bulkevols::NTuple)

Build a container to store all the evolved quantities as elements of an
`NTuple`. The idea is to treat them as a single column vector for the point of
view of the time evolution routine.
"""
function EvolVars(bulkevols::BulkPartition{Nsys}) where {Nsys}
    EvolVars(bulkevols.x)
end

@inline getbulkevolvedpartition(ff::EvolVars) = BulkPartition(ff.x)

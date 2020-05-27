
"""
Extend this type for different `Potential` choices
"""
abstract type Potential end

"""
Extend this type for different `InitialData` choices
"""
abstract type InitialData end

"""
Extend this type for different `EvolutionEquations`
"""
abstract type EvolutionEquations end

struct EvolTest0 <: EvolutionEquations end

Base.@kwdef struct AffineNull{T,TP<:Potential} <: EvolutionEquations
    phi0          :: T   = 0.0
    potential     :: TP  = ZeroPotential()
end


abstract type AbstractVars{T} <: AbstractVector{T} end

struct BulkEvolved{T} <: AbstractVars{T}
    B1  :: Array{T,3}
    B2  :: Array{T,3}
    G   :: Array{T,3}
    phi :: Array{T,3}
end

struct BulkConstrained{T} <: AbstractVars{T}
    S    :: Array{T,3}
    Fx   :: Array{T,3}
    Fy   :: Array{T,3}
    B1d  :: Array{T,3}
    B2d  :: Array{T,3}
    Gd   :: Array{T,3}
    phid :: Array{T,3}
    Sd   :: Array{T,3}
    A    :: Array{T,3}
end

struct Bulk{T} <: AbstractVars{T}
    B1   :: Array{T,3}
    B2   :: Array{T,3}
    G    :: Array{T,3}
    phi  :: Array{T,3}
    S    :: Array{T,3}
    Fx   :: Array{T,3}
    Fy   :: Array{T,3}
    B1d  :: Array{T,3}
    B2d  :: Array{T,3}
    Gd   :: Array{T,3}
    phid :: Array{T,3}
    Sd   :: Array{T,3}
    A    :: Array{T,3}
end

struct Boundary{T} <: AbstractVars{T}
    a4  :: Array{T,3}
    fx2 :: Array{T,3}
    fy2 :: Array{T,3}
end

struct Gauge{T} <: AbstractVars{T}
    xi  :: Array{T,3}
end

@inline varlist(::BulkEvolved)     = (:B1, :B2, :G, :phi)
@inline varlist(::BulkConstrained) = (:S, :Fx, :Fy, :B1d, :B2d, :Gd, :phid, :Sd, :A)
@inline varlist(::Bulk)            = (:B1, :B2, :G, :phi, :S, :Fx, :Fy, :B1d, :B2d,
                                      :Gd, :phid, :Sd, :A)
@inline varlist(::Boundary)        = (:a4, :fx2, :fy2)
@inline varlist(::Gauge)           = (:xi,)


"""
    BulkEvolved{T}(undef, Nu, Nx, Ny)

Construct a container of uninitialized Arrays to hold all the bulk variables
that are evolved in time: B1, B2, G, phi
"""
function BulkEvolved{T}(::UndefInitializer, Nu::Int, Nx::Int, Ny::Int) where {T<:Real}
    B1  = Array{T}(undef, Nu, Nx, Ny)
    B2  = Array{T}(undef, Nu, Nx, Ny)
    G   = Array{T}(undef, Nu, Nx, Ny)
    phi = Array{T}(undef, Nu, Nx, Ny)
    BulkEvolved{T}(B1, B2, G, phi)
end

"""
    BulkConstrained{T}(undef, Nu, Nx, Ny)

Construct a container of uninitialized Arrays to hold all the bulk variables
that are constrained (not evolved in time): S, Fx, Fy, B1d, B2d, Gd, phid, Sd, A
"""
function BulkConstrained{T}(::UndefInitializer, Nu::Int, Nx::Int, Ny::Int) where {T<:Real}
    S    = Array{T}(undef, Nu, Nx, Ny)
    Fx   = Array{T}(undef, Nu, Nx, Ny)
    Fy   = Array{T}(undef, Nu, Nx, Ny)
    B1d  = Array{T}(undef, Nu, Nx, Ny)
    B2d  = Array{T}(undef, Nu, Nx, Ny)
    Gd   = Array{T}(undef, Nu, Nx, Ny)
    phid = Array{T}(undef, Nu, Nx, Ny)
    Sd   = Array{T}(undef, Nu, Nx, Ny)
    A    = Array{T}(undef, Nu, Nx, Ny)
    BulkConstrained{T}(S, Fx, Fy, B1d, B2d, Gd, phid, Sd, A)
end

"""
    Bulk(bulkevol::BulkEvolved, bulkconstrain::BulkConstrained)

Construct a container to hold all the bulk variables where the evolved variables
point to the given bulkevol struct and the constrained ones point to the bulkconstrain struct
"""
function Bulk(bulkevol::BulkEvolved{T}, bulkconstrain::BulkConstrained{T}) where {T}
    B1    = bulkevol.B1
    B2    = bulkevol.B2
    G     = bulkevol.G
    phi   = bulkevol.phi
    S     = bulkconstrain.S
    Fx    = bulkconstrain.Fx
    Fy    = bulkconstrain.Fy
    B1d   = bulkconstrain.B1d
    B2d   = bulkconstrain.B2d
    Gd    = bulkconstrain.Gd
    phid  = bulkconstrain.phid
    Sd    = bulkconstrain.Sd
    A     = bulkconstrain.A
    Bulk{T}(B1, B2, G, phi, S, Fx, Fy, B1d, B2d, Gd, phid, Sd, A)
end

"""
    Boundary{T}(undef, Nx, Ny)

Construct a container of uninitialized Arrays to hold all the boundary
variables: a4, fx2, fy2

These variables are automatically defined on a `(1,Nx,Ny)` grid, rather than
a `(Nx,Ny)` one, so that the same Dx and Dy differential operators defined for
the bulk quantities can also straightforwardly apply on them. Remember that the
axis along which the operator applies is defined on the operator itself. So, by
defining things this way, the Dx operator (which acts along the 2nd index) will
also do the correct thing when acting on a4, fx2, fy2.
"""
function Boundary{T}(::UndefInitializer, Nx::Int, Ny::Int) where {T<:Real}
    a4  = Array{T}(undef, 1, Nx, Ny)
    fx2 = Array{T}(undef, 1, Nx, Ny)
    fy2 = Array{T}(undef, 1, Nx, Ny)
    Boundary{T}(a4, fx2, fy2)
end

"""
    Gauge{T}(undef, Nx, Ny)

Construct a container with an uninitialized Array to hold the gauge variable xi

As in the Boundary struct, xi is automatically defined on a `(1,Nx,Ny)` grid,
rather than a `(Nx,Ny)` one, so that the same Dx and Dy differential operators
defined for the bulk quantities can also straightforwardly apply on it.
"""
function Gauge{T}(::UndefInitializer, Nx::Int, Ny::Int) where {T<:Real}
    xi  = Array{T}(undef, 1, Nx, Ny)
    Gauge{T}(xi)
end


getB1(ff::BulkEvolved)       = ff.B1
getB2(ff::BulkEvolved)       = ff.B2
getG(ff::BulkEvolved)        = ff.G
getphi(ff::BulkEvolved)      = ff.phi

getS(ff::BulkConstrained)    = ff.S
getFx(ff::BulkConstrained)   = ff.Fx
getFy(ff::BulkConstrained)   = ff.Fy
getB1d(ff::BulkConstrained)  = ff.B1d
getB2d(ff::BulkConstrained)  = ff.B2d
getGd(ff::BulkConstrained)   = ff.Gd
getphid(ff::BulkConstrained) = ff.phid
getA(ff::BulkConstrained)    = ff.A

getB1(ff::Bulk)              = ff.B1
getB2(ff::Bulk)              = ff.B2
getG(ff::Bulk)               = ff.G
getphi(ff::Bulk)             = ff.phi
getS(ff::Bulk)               = ff.S
getFx(ff::Bulk)              = ff.Fx
getFy(ff::Bulk)              = ff.Fy
getB1d(ff::Bulk)             = ff.B1d
getB2d(ff::Bulk)             = ff.B2d
getGd(ff::Bulk)              = ff.Gd
getphid(ff::Bulk)            = ff.phid
getA(ff::Bulk)               = ff.A

geta4(ff::Boundary)          = ff.a4
getfx2(ff::Boundary)         = ff.fx2
getfy2(ff::Boundary)         = ff.fy2

getxi(ff::Gauge)             = ff.xi


@inline function Base.length(ff::AbstractVars)
    vars = varlist(ff)
    sum_l = 0
    for x in vars
        f = getproperty(ff,x)   # f will point to each variable in the ff struct
        sum_l += length(f)
    end
    sum_l
end
@inline Base.size(ff::AbstractVars) = (length(ff),)

# indexing. this is just a linear indexing through all the arrays. adapted from
# RecursiveArrayTools
@inline Base.firstindex(ff::AbstractVars) = 1
@inline Base.lastindex(ff::AbstractVars) = length(ff)

@inline function Base.getindex(ff::AbstractVars, i::Int)
    vars = varlist(ff)
    @inbounds for x in vars
        f  = getproperty(ff,x)
        i -= length(f)
        if i <= 0
            return f[length(f)+i]
        end
    end
end
@inline function Base.setindex!(ff::AbstractVars, v, i::Int)
    vars = varlist(ff)
    @inbounds for x in vars
        f  = getproperty(ff,x)
        i -= length(f)
        if i <= 0
            f[length(f)+i] = v
            break
        end
    end
end

"""
    unpack(ff::AbstractVars)

Return `Array` with all the different fields in the structure
"""
function unpack(ff::AbstractVars)
    vars = varlist(ff)
    [getproperty(ff,x) for x in vars]
end


struct BulkPartition{N,A}
    _x :: NTuple{N,A}
end
"""
    BulkPartition(x...)

Container to store (bulk) quantities that may be spread across different grid
partitions. This is to be thought as a `Tuple` (or an `Array`) of `Bulk` objects.
"""
BulkPartition(x...) = BulkPartition(tuple(x...))

@inline Base.iterate(ff::BulkPartition)         = iterate(ff._x)
@inline Base.iterate(ff::BulkPartition, i::Int) = iterate(ff._x, i)

@inline Base.length(ff::BulkPartition{N}) where{N} = N
@inline Base.firstindex(ff::BulkPartition) = 1
@inline Base.lastindex(ff::BulkPartition)  = length(ff)

@inline Base.getindex(ff::BulkPartition, i::Int) = ff._x[i]

function BulkEvolved(bulks::BulkPartition{N}) where{N}
    f = ntuple(i -> BulkEvolved(bulks[i]), N)
    BulkPartition(f)
end

function BulkConstrained(bulks::BulkPartition{N}) where{N}
    f = ntuple(i -> BulkConstrained(bulks[i]), N)
    BulkPartition(f)
end



struct EvolVars{T,N,A} <: AbstractVector{T}
    x :: NTuple{N,A}
end
EvolVars(x::NTuple{N,A}) where {N,A} = EvolVars{eltype(A),N,A}(x)
EvolVars(x...) = EvolVars(tuple(x...))


function EvolVars(ff::AbstractVector{A}) where{A<:AbstractArray}
    x = Tuple(ff)
    EvolVars{eltype(A),length(x),A}(x)
end

"""
    EvolVars(boundary::Boundary, gauge::Gauge, bulkevols::NTuple)

Build a container to store all the evolved quantities as elements of an
`NTuple`. The idea is to treat them as a single column vector for the point of
view of the time evolution routine. Inspired in `ArrayPartition` from
`RecursiveArrayTools`
"""
function EvolVars(boundary::Boundary, gauge::Gauge,
                  bulkevols::BulkPartition{Nsys}) where {Nsys}
    f1 = unpack(boundary)
    f2 = unpack(gauge)

    # this technique allocates a bit of memory and is not fully type-stable, but
    # i think it's not a problem since we only call this function once
    f3 = [unpack(bulkevol) for bulkevol in bulkevols]
    EvolVars([f1; f2; f3...])
end

Base.similar(ff::EvolVars{T,N,S}) where {T,N,S} = EvolVars{T,N,S}(similar.(ff.x))

@inline Base.length(ff::EvolVars) = sum((length(x) for x in ff.x))
@inline Base.size(ff::EvolVars)   = (length(ff),)

# indexing. this is just a linear indexing through all the arrays. adapted from
# RecursiveArrayTools
@inline Base.firstindex(ff::EvolVars) = 1
@inline Base.lastindex(ff::EvolVars)  = length(ff)

@inline function Base.getindex(ff::EvolVars, i::Int)
    @inbounds for j in 1:length(ff.x)
        f  = ff.x[j]
        i -= length(f)
        if i <= 0
            return f[length(f)+i]
        end
    end
end
@inline function Base.setindex!(ff::EvolVars, v, i::Int)
    @inbounds for j in 1:length(ff.x)
        f  = ff.x[j]
        i -= length(f)
        if i <= 0
            f[length(f)+i] = v
            break
        end
    end
end

@inline npartitions(::EvolVars{T,N}) where {T,N} = N
@inline getudomains(::EvolVars{T,N}) where {T,N} = div(N-4, 4)

geta4(ff::EvolVars)   = ff.x[1]
getfx2(ff::EvolVars)  = ff.x[2]
getfy2(ff::EvolVars)  = ff.x[3]
getxi(ff::EvolVars)   = ff.x[4]

function getB1(ff::EvolVars, i::Int)
    Nsys = getudomains(ff)
    @assert i > 0
    @assert i <= Nsys
    ff.x[5 + (i-1)*4]
end

function getB2(ff::EvolVars, i::Int)
    Nsys = getudomains(ff)
    @assert i > 0
    @assert i <= Nsys
    ff.x[6 + (i-1)*4]
end

function getG(ff::EvolVars, i::Int)
    Nsys = getudomains(ff)
    @assert i > 0
    @assert i <= Nsys
    ff.x[7 + (i-1)*4]
end

function getphi(ff::EvolVars, i::Int)
    Nsys = getudomains(ff)
    @assert i > 0
    @assert i <= Nsys
    ff.x[8 + (i-1)*4]
end

@inline getboundary(ff::EvolVars) = Boundary(geta4(ff), getfx2(ff), getfy2(ff))
@inline getgauge(ff::EvolVars)    = Gauge(getxi(ff))

@inline getbulkevolved(ff::EvolVars, i::Int) =
    BulkEvolved(getB1(ff,i), getB2(ff,i), getG(ff,i), getphi(ff,i))

function getbulkevolvedpartition(ff::EvolVars)
    Nsys = getudomains(ff)
    f = ntuple(i -> getbulkevolved(ff,i), Nsys)
    BulkPartition(f)
end

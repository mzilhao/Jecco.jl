
"""
Extend this type for different potential choices
"""
abstract type Potential end

"""
Extend this type for all IBVP parameters, including different initial
conditions. It needs the fields `phi0` and `potential`.
"""
abstract type IBVP end

"""
Extend this type for different evolution equations
"""
abstract type AbstractEvolEq end


abstract type AbstractVars{T} <: AbstractVector{T} end

abstract type EvolVars{T} <: AbstractVars{T} end


struct BulkEvol{T} <: AbstractVars{T}
    B1  :: Array{T,3}
    B2  :: Array{T,3}
    G   :: Array{T,3}
    phi :: Array{T,3}
end

struct Boundary{T} <: AbstractVars{T}
    a4  :: Array{T,3}
    fx2 :: Array{T,3}
    fy2 :: Array{T,3}
end

struct Gauge{T} <: AbstractVars{T}
    xi  :: Array{T,3}
end

struct Bulk{T,N} <: AbstractVars{T}
    B1   :: Array{T,N}
    B2   :: Array{T,N}
    G    :: Array{T,N}
    phi  :: Array{T,N}
    S    :: Array{T,N}
    Fx   :: Array{T,N}
    Fy   :: Array{T,N}
    B1d  :: Array{T,N}
    B2d  :: Array{T,N}
    Gd   :: Array{T,N}
    phid :: Array{T,N}
    Sd   :: Array{T,N}
    A    :: Array{T,N}
end

@inline varlist(::BulkEvol) = [:B1, :B2, :G, :phi]
@inline varlist(::Boundary) = [:a4, :fx2, :fy2]
@inline varlist(::Gauge)    = [:xi]
@inline varlist(::Bulk)     = [:B1, :B2, :G, :phi, :S, :Fx, :Fy, :B1d, :B2d,
                               :Gd, :phid, :Sd, :A]

"""
    BulkEvol{T}(undef, Nu, Nx, Ny)

Construct a container of uninitialized Arrays to hold all the bulk variables
that are evolved in time: B1, B2, G, phi
"""
function BulkEvol{T}(::UndefInitializer, Nu::Int, Nx::Int, Ny::Int) where {T<:Real}
    B1  = Array{T}(undef, Nu, Nx, Ny)
    B2  = Array{T}(undef, Nu, Nx, Ny)
    G   = Array{T}(undef, Nu, Nx, Ny)
    phi = Array{T}(undef, Nu, Nx, Ny)
    BulkEvol{T}(B1, B2, G, phi)
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

"""
    Bulk{T}(undef, Nxx...)

Construct a container of uninitialized Arrays to hold all the bulk variables:
B1, B2, G, phi, S, Fx, Fy, B1d, B2d, Gd, phid, Sd, A

"""
function Bulk{T}(::UndefInitializer, Nxx::Vararg{Int,N}) where {T<:Real,N}
    B1   = Array{T}(undef, Nxx...)
    B2   = Array{T}(undef, Nxx...)
    G    = Array{T}(undef, Nxx...)
    phi  = Array{T}(undef, Nxx...)
    S    = Array{T}(undef, Nxx...)
    Fx   = Array{T}(undef, Nxx...)
    Fy   = Array{T}(undef, Nxx...)
    B1d  = Array{T}(undef, Nxx...)
    B2d  = Array{T}(undef, Nxx...)
    Gd   = Array{T}(undef, Nxx...)
    phid = Array{T}(undef, Nxx...)
    Sd   = Array{T}(undef, Nxx...)
    A    = Array{T}(undef, Nxx...)
    Bulk{T,N}(B1, B2, G, phi, S, Fx, Fy, B1d, B2d, Gd, phid, Sd, A)
end
"""
    Bulk(bulkevol::BulkEvol)

Construct a container to hold all the bulk variables, but where the evolved ones
point to the given bulkevol struct
"""
function Bulk(ff::BulkEvol{T}) where {T}
    B1    = ff.B1
    B2    = ff.B2
    G     = ff.G
    phi   = ff.phi
    S     = similar(B1)
    Fx    = similar(B1)
    Fy    = similar(B1)
    B1d   = similar(B1)
    B2d   = similar(B1)
    Gd    = similar(B1)
    phid  = similar(B1)
    Sd    = similar(B1)
    A     = similar(B1)
    Bulk(B1, B2, G, phi, S, Fx, Fy, B1d, B2d, Gd, phid, Sd, A)
end

getB1(ff::BulkEvol)  = ff.B1
getB2(ff::BulkEvol)  = ff.B2
getG(ff::BulkEvol)   = ff.G
getphi(ff::BulkEvol) = ff.phi

geta4(ff::Boundary)  = ff.a4
getfx2(ff::Boundary) = ff.fx2
getfy2(ff::Boundary) = ff.fy2

getxi(ff::Gauge)     = ff.xi

getB1(ff::Bulk)      = ff.B1
getB2(ff::Bulk)      = ff.B2
getG(ff::Bulk)       = ff.G
getphi(ff::Bulk)     = ff.phi
getS(ff::Bulk)       = ff.S
getFx(ff::Bulk)      = ff.Fx
getFy(ff::Bulk)      = ff.Fy
getB1d(ff::Bulk)     = ff.B1d
getB2d(ff::Bulk)     = ff.B2d
getGd(ff::Bulk)      = ff.Gd
getphid(ff::Bulk)    = ff.phid
getA(ff::Bulk)       = ff.A

Base.similar(ff::BulkEvol) = BulkEvol(similar(ff.B1), similar(ff.B2), similar(ff.G), similar(ff.phi))
Base.similar(ff::Boundary) = Boundary(similar(ff.a4), similar(ff.fx2), similar(ff.fy2))
Base.similar(ff::Gauge)    = Gauge(similar(ff.xi))
Base.similar(ff::Bulk)     =
    BulkAll(similar(ff.B1), similar(ff.B2), similar(ff.G), similar(ff.phi),
            similar(ff.S), similar(ff.Fx), similar(ff.Fy), similar(ff.B1d),
            similar(ff.B2d), similar(ff.Gd), similar(phid), similar(Sd),
            similar(ff.A))

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


struct EvolPartition{T,N,A} <: EvolVars{T}
    x :: NTuple{N,A}
end
function EvolPartition(ff::AbstractVector{A}) where{A<:AbstractArray}
    x = Tuple(ff)
    EvolPartition{eltype(A),length(x),A}(x)
end

"""
    EvolPartition(boundary::Boundary, gauge::Gauge, bulkevols::NTuple)

Build a container to store all the evolved quantities as elements of an
`NTuple`. The idea is to treat them as a single column vector for the point of
view of the time evolution routine. Inspired in `ArrayPartition` from
`RecursiveArrayTools`
"""
function EvolPartition(boundary::Boundary{T}, gauge::Gauge{T},
                       bulkevols::NTuple{Nsys,BulkEvol{T}}) where {T,Nsys}
    f1 = unpack(boundary)
    f2 = unpack(gauge)
    f3 = [unpack(bulkevol) for bulkevol in bulkevols]
    EvolPartition([f1; f2; f3...])
end

Base.similar(ff::EvolPartition{T,N,S}) where {T,N,S} = EvolPartition{T,N,S}(similar.(ff.x))

@inline Base.length(ff::EvolPartition) = sum((length(x) for x in ff.x))
@inline Base.size(ff::EvolPartition)   = (length(ff),)

# indexing. this is just a linear indexing through all the arrays. adapted from
# RecursiveArrayTools
@inline Base.firstindex(ff::EvolPartition) = 1
@inline Base.lastindex(ff::EvolPartition)  = length(ff)

@inline function Base.getindex(ff::EvolPartition, i::Int)
    @inbounds for j in 1:length(ff.x)
        f  = ff.x[j]
        i -= length(f)
        if i <= 0
            return f[length(f)+i]
        end
    end
end
@inline function Base.setindex!(ff::EvolPartition, v, i::Int)
    @inbounds for j in 1:length(ff.x)
        f  = ff.x[j]
        i -= length(f)
        if i <= 0
            f[length(f)+i] = v
            break
        end
    end
end

@inline npartitions(::EvolPartition{T,N}) where {T,N} = N
@inline get_udomains(::EvolPartition{T,N}) where {T,N} = div(N-4, 4)

geta4(ff::EvolPartition)   = ff.x[1]
getfx2(ff::EvolPartition)  = ff.x[2]
getfy2(ff::EvolPartition)  = ff.x[3]
getxi(ff::EvolPartition)   = ff.x[4]

function getB1(ff::EvolPartition, i::Int)
    Nsys = get_udomains(ff)
    @assert i > 0
    @assert i <= Nsys
    ff.x[5 + (i-1)*4]
end

function getB2(ff::EvolPartition, i::Int)
    Nsys = get_udomains(ff)
    @assert i > 0
    @assert i <= Nsys
    ff.x[6 + (i-1)*4]
end

function getG(ff::EvolPartition, i::Int)
    Nsys = get_udomains(ff)
    @assert i > 0
    @assert i <= Nsys
    ff.x[7 + (i-1)*4]
end

function getphi(ff::EvolPartition, i::Int)
    Nsys = get_udomains(ff)
    @assert i > 0
    @assert i <= Nsys
    ff.x[8 + (i-1)*4]
end

@inline getboundary(ff::EvolPartition) = Boundary(geta4(ff), getfx2(ff), getfy2(ff))
@inline getgauge(ff::EvolPartition) = Gauge(getxi(ff))

@inline getbulkevol(ff::EvolPartition, i::Int) =
    BulkEvol(getB1(ff,i), getB2(ff,i), getG(ff,i), getphi(ff,i))

function getbulkevols(ff::EvolPartition)
    Nsys = get_udomains(ff)
    [getbulkevol(ff,i) for i in 1:Nsys]
end

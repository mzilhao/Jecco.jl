
abstract type GridType end
struct Inner <: GridType end
struct Outer <: GridType end

"""
Extend this type for different potential choices
"""
abstract type Potential end

"""
Extend this type for different initial conditions
"""
abstract type IBVP{T} end


abstract type AbstractVars{T} <: AbstractVector{T} end

abstract type EvolVars{T} <: AbstractVars{T} end

struct BulkEvol{T} <: EvolVars{T}
    B1  :: Array{T,3}
    B2  :: Array{T,3}
    G   :: Array{T,3}
    phi :: Array{T,3}
end

struct Boundary{T} <: EvolVars{T}
    a4  :: Array{T,3}
    fx2 :: Array{T,3}
    fy2 :: Array{T,3}
end

struct Gauge{T} <: EvolVars{T}
    xi  :: Array{T,3}
end

@inline varlist(::BulkEvol) = [:B1, :B2, :G, :phi]
@inline varlist(::Boundary) = [:a4, :fx2, :fy2]
@inline varlist(::Gauge)    = [:xi]


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

Base.similar(ff::BulkEvol) = BulkEvol(similar(ff.B1), similar(ff.B2), similar(ff.G), similar(ff.phi))
Base.similar(ff::Boundary) = Boundary(similar(ff.a4), similar(ff.fx2), similar(ff.fy2))
Base.similar(ff::Gauge)    = Gauge(similar(ff.xi))


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

@inline function Base.getindex(evol::AbstractVars, i::Int)
    vars = varlist(evol)
    @inbounds for x in vars
        f  = getproperty(evol,x)
        i -= length(f)
        if i <= 0
            return f[length(f)+i]
        end
    end
end
@inline function Base.setindex!(evol::AbstractVars, v, i::Int)
    vars = varlist(evol)
    @inbounds for x in vars
        f  = getproperty(evol,x)
        i -= length(f)
        if i <= 0
            f[length(f)+i] = v
            break
        end
    end
end

getB1(ff::BulkEvol)  = ff.B1
getB2(ff::BulkEvol)  = ff.B2
getG(ff::BulkEvol)   = ff.G
getphi(ff::BulkEvol) = ff.phi

geta4(ff::Boundary)  = ff.a4
getfx2(ff::Boundary) = ff.fx2
getfy2(ff::Boundary) = ff.fy2

getxi(ff::Gauge)     = ff.xi


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

@inline varlist(::Bulk) = [:B1, :B2, :G, :phi, :S, :Fx, :Fy, :B1d, :B2d,
                           :Gd, :phid, :Sd, :A]

"""
    Bulk{T}(undef, Nu, Nx, Ny)

Construct a container of uninitialized Arrays to hold all the bulk variables:
B1, B2, G, phi, S, Fx, Fy, B1d, B2d, Gd, phid, Sd, A

"""
function Bulk{T}(::UndefInitializer, Nu::Int, Nx::Int, Ny::Int) where {T<:Real}
    B1   = Array{T}(undef, Nu, Nx, Ny)
    B2   = Array{T}(undef, Nu, Nx, Ny)
    G    = Array{T}(undef, Nu, Nx, Ny)
    phi  = Array{T}(undef, Nu, Nx, Ny)
    S    = Array{T}(undef, Nu, Nx, Ny)
    Fx   = Array{T}(undef, Nu, Nx, Ny)
    Fy   = Array{T}(undef, Nu, Nx, Ny)
    B1d  = Array{T}(undef, Nu, Nx, Ny)
    B2d  = Array{T}(undef, Nu, Nx, Ny)
    Gd   = Array{T}(undef, Nu, Nx, Ny)
    phid = Array{T}(undef, Nu, Nx, Ny)
    Sd   = Array{T}(undef, Nu, Nx, Ny)
    A    = Array{T}(undef, Nu, Nx, Ny)
    Bulk{T}(B1, B2, G, phi, S, Fx, Fy, B1d, B2d, Gd, phid, Sd, A)
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
    Bulk{T}(B1, B2, G, phi, S, Fx, Fy, B1d, B2d, Gd, phid, Sd, A)
end

getB1(ff::Bulk)   = ff.B1
getB2(ff::Bulk)   = ff.B2
getG(ff::Bulk)    = ff.G
getphi(ff::Bulk)  = ff.phi
getS(ff::Bulk)    = ff.S
getFx(ff::Bulk)   = ff.Fx
getFy(ff::Bulk)   = ff.Fy
getB1d(ff::Bulk)  = ff.B1d
getB2d(ff::Bulk)  = ff.B2d
getGd(ff::Bulk)   = ff.Gd
getphid(ff::Bulk) = ff.phid
getA(ff::Bulk)    = ff.A

Base.similar(ff::Bulk) =
    BulkAll(similar(ff.B1), similar(ff.B2), similar(ff.G), similar(ff.phi),
            similar(ff.S), similar(ff.Fx), similar(ff.Fy), similar(ff.B1d),
            similar(ff.B2d), similar(ff.Gd), similar(phid), similar(Sd),
            similar(ff.A))



struct EvolPartition{T,N,S<:NTuple} <: AbstractVector{T}
    x :: S
end

function EvolPartition(ff::AbstractVector{A}) where{A}
    x = Tuple(ff)
    EvolPartition{eltype(A),length(x),typeof(x)}(x)
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




#####


struct BulkVars{T}
    B1       :: T
    B2       :: T
    G        :: T
    phi      :: T
    S        :: T
    Fx       :: T
    Fy       :: T
    B1d      :: T
    B2d      :: T
    Gd       :: T
    phid     :: T
    Sd       :: T
    A        :: T
end

function BulkVars(::Type{T}, Nxx::Vararg) where{T}
    B1     = zeros(T, Nxx...)
    B2     = copy(B1)
    G      = copy(B1)
    phi    = copy(B1)
    S      = copy(B1)
    Fx     = copy(B1)
    Fy     = copy(B1)
    B1d    = copy(B1)
    B2d    = copy(B1)
    Gd     = copy(B1)
    phid   = copy(B1)
    Sd     = copy(B1)
    A      = copy(B1)
    BulkVars{typeof(B1)}(B1, B2, G, phi, S, Fx, Fy, B1d, B2d, Gd, phid, Sd, A)
end

function BulkVars(ff::EvolVars)
    B1    = ff.B1
    B2    = ff.B2
    G     = ff.G
    phi   = ff.phi
    a4    = ff.a4
    fx2   = ff.fx2
    fy2   = ff.fy2
    xi    = ff.xi
    S     = similar(B1)
    Fx    = similar(B1)
    Fy    = similar(B1)
    B1d   = similar(B1)
    B2d   = similar(B1)
    Gd    = similar(B1)
    phid  = similar(B1)
    Sd    = similar(B1)
    A     = similar(B1)
    BulkVars{typeof(B1)}(gridtype, B1, B2, G, phi, S, Fx, Fy, B1d, B2d, Gd, phid, Sd, A)
end

function BulkVars(B1::T, B2::T, G::T, phi::T) where {T}
    S      = similar(B1)
    Fx     = similar(B1)
    Fy     = similar(B1)
    B1d    = similar(B1)
    B2d    = similar(B1)
    Gd     = similar(B1)
    phid   = similar(B1)
    Sd     = similar(B1)
    A      = similar(B1)

    BulkVars{T}(B1, B2, G, phi, S, Fx, Fy, B1d, B2d, Gd, phid, Sd, A)
end

#function BulkVars(gridtypes::Vector{GT}, ffs::Vector{T}) where {T<:EvolVars}
#    [BulkVars(gridtype, ff) for (gridtype, ff) in (gridtypes, ffs)]
#end

#function BulkVars(systems::Vector{T1}, ffs::Vector{T2}) where {T1<:System,T2<:EvolVars}
#    [BulkVars(sys.gridtype, ff) for (sys, ff) in (systems, ffs)]
#end


struct GaugeVars{A,T}
    xi    :: A
    kappa :: T
end

function GaugeVars(xi::Array{T,N}, kappa::T) where {T<:Number,N}
    GaugeVars{typeof(xi), typeof(kappa)}(xi, kappa)
end

struct BaseVars{PT,T}
    potential :: PT
    phi0      :: T
end


struct BoundaryVars{T}
    a4   :: T
    fx2  :: T
    fy2  :: T
end
BoundaryVars(a4, fx2, fy2) = BoundaryVars{typeof(a4)}(a4, fx2, fy2)


# TODO: use named tuples for these?

struct SVars{T<:Real}
    phi0     :: T

    u        :: T

    xi       :: T

    B1       :: T
    B1p      :: T

    B2       :: T
    B2p      :: T

    G        :: T
    Gp       :: T

    phi      :: T
    phip     :: T
end

struct FVars{T<:Real}
    phi0     :: T

    u        :: T

    xi       :: T
    xi_x     :: T
    xi_y     :: T

    B1       :: T
    B2       :: T
    G        :: T
    phi      :: T
    S        :: T

    B1p      :: T
    B2p      :: T
    Gp       :: T
    phip     :: T
    Sp       :: T

    B1pp     :: T
    B2pp     :: T
    Gpp      :: T
    Spp      :: T

    B1_x     :: T
    B2_x     :: T
    G_x      :: T
    phi_x    :: T
    S_x      :: T

    B1_y     :: T
    B2_y     :: T
    G_y      :: T
    phi_y    :: T
    S_y      :: T

    B1p_x    :: T
    B2p_x    :: T
    Gp_x     :: T
    Sp_x     :: T

    B1p_y    :: T
    B2p_y    :: T
    Gp_y     :: T
    Sp_y     :: T
end


struct SdVars{T<:Real,PT}
    potential:: PT

    phi0     :: T

    u        :: T

    xi       :: T
    xi_x     :: T
    xi_y     :: T
    xi_xx    :: T
    xi_yy    :: T
    xi_xy    :: T

    B1       :: T
    B2       :: T
    G        :: T
    phi      :: T
    S        :: T
    Fx       :: T
    Fy       :: T

    B1p      :: T
    B2p      :: T
    Gp       :: T
    phip     :: T
    Sp       :: T
    Fxp      :: T
    Fyp      :: T

    B1pp     :: T
    B2pp     :: T
    Gpp      :: T
    phipp    :: T
    Spp      :: T
    Fxpp     :: T
    Fypp     :: T

    B1_x     :: T
    B2_x     :: T
    G_x      :: T
    phi_x    :: T
    S_x      :: T
    Fx_x     :: T
    Fy_x     :: T

    B1_y     :: T
    B2_y     :: T
    G_y      :: T
    phi_y    :: T
    S_y      :: T
    Fx_y     :: T
    Fy_y     :: T

    B1p_x    :: T
    B2p_x    :: T
    Gp_x     :: T
    phip_x   :: T
    Sp_x     :: T
    Fxp_x    :: T
    Fyp_x    :: T

    B1p_y    :: T
    B2p_y    :: T
    Gp_y     :: T
    phip_y   :: T
    Sp_y     :: T
    Fxp_y    :: T
    Fyp_y    :: T

    B1_xx    :: T
    B2_xx    :: T
    G_xx     :: T
    phi_xx   :: T
    S_xx     :: T

    B1_yy    :: T
    B2_yy    :: T
    G_yy     :: T
    phi_yy   :: T
    S_yy     :: T

    B2_xy    :: T
    G_xy     :: T
    S_xy     :: T
end


struct BdGVars{T<:Real}
    phi0     :: T

    u        :: T

    xi       :: T
    xi_x     :: T
    xi_y     :: T
    xi_xx    :: T
    xi_yy    :: T
    xi_xy    :: T

    B1       :: T
    B2       :: T
    G        :: T
    phi      :: T
    S        :: T
    Fx       :: T
    Fy       :: T
    Sd       :: T

    B1p      :: T
    B2p      :: T
    Gp       :: T
    phip     :: T
    Sp       :: T
    Fxp      :: T
    Fyp      :: T

    B1pp     :: T
    B2pp     :: T
    Gpp      :: T
    phipp    :: T
    Spp      :: T
    Fxpp     :: T
    Fypp     :: T

    B1_x     :: T
    B2_x     :: T
    G_x      :: T
    phi_x    :: T
    S_x      :: T
    Fx_x     :: T
    Fy_x     :: T

    B1_y     :: T
    B2_y     :: T
    G_y      :: T
    phi_y    :: T
    S_y      :: T
    Fx_y     :: T
    Fy_y     :: T

    B1p_x    :: T
    B2p_x    :: T
    Gp_x     :: T
    phip_x   :: T
    Sp_x     :: T
    Fxp_x    :: T
    Fyp_x    :: T

    B1p_y    :: T
    B2p_y    :: T
    Gp_y     :: T
    phip_y   :: T
    Sp_y     :: T
    Fxp_y    :: T
    Fyp_y    :: T

    B1_xx    :: T
    B2_xx    :: T
    G_xx     :: T
    phi_xx   :: T
    S_xx     :: T

    B1_yy    :: T
    B2_yy    :: T
    G_yy     :: T
    phi_yy   :: T
    S_yy     :: T

    B2_xy    :: T
    G_xy     :: T
    S_xy     :: T
end


struct phidVars{T<:Real,PT}
    potential:: PT

    phi0     :: T

    u        :: T

    xi       :: T
    xi_x     :: T
    xi_y     :: T
    xi_xx    :: T
    xi_yy    :: T
    xi_xy    :: T

    B1       :: T
    B2       :: T
    G        :: T
    phi      :: T
    S        :: T
    Fx       :: T
    Fy       :: T
    Sd       :: T

    B1p      :: T
    B2p      :: T
    Gp       :: T
    phip     :: T
    Sp       :: T
    Fxp      :: T
    Fyp      :: T

    B1pp     :: T
    B2pp     :: T
    Gpp      :: T
    phipp    :: T
    Spp      :: T
    Fxpp     :: T
    Fypp     :: T

    B1_x     :: T
    B2_x     :: T
    G_x      :: T
    phi_x    :: T
    S_x      :: T
    Fx_x     :: T
    Fy_x     :: T

    B1_y     :: T
    B2_y     :: T
    G_y      :: T
    phi_y    :: T
    S_y      :: T
    Fx_y     :: T
    Fy_y     :: T

    B1p_x    :: T
    B2p_x    :: T
    Gp_x     :: T
    phip_x   :: T
    Sp_x     :: T
    Fxp_x    :: T
    Fyp_x    :: T

    B1p_y    :: T
    B2p_y    :: T
    Gp_y     :: T
    phip_y   :: T
    Sp_y     :: T
    Fxp_y    :: T
    Fyp_y    :: T

    B1_xx    :: T
    B2_xx    :: T
    G_xx     :: T
    phi_xx   :: T
    S_xx     :: T

    B1_yy    :: T
    B2_yy    :: T
    G_yy     :: T
    phi_yy   :: T
    S_yy     :: T

    B2_xy    :: T
    G_xy     :: T
    phi_xy   :: T
    S_xy     :: T
end


struct AVars{T<:Real,PT}
    potential:: PT

    phi0     :: T

    u        :: T

    xi       :: T
    xi_x     :: T
    xi_y     :: T
    xi_xx    :: T
    xi_yy    :: T
    xi_xy    :: T

    B1       :: T
    B2       :: T
    G        :: T
    phi      :: T
    S        :: T
    Fx       :: T
    Fy       :: T
    Sd       :: T
    B1d      :: T
    B2d      :: T
    Gd       :: T
    phid     :: T

    B1p      :: T
    B2p      :: T
    Gp       :: T
    phip     :: T
    Sp       :: T
    Fxp      :: T
    Fyp      :: T

    B1pp     :: T
    B2pp     :: T
    Gpp      :: T
    phipp    :: T
    Spp      :: T
    Fxpp     :: T
    Fypp     :: T

    B1_x     :: T
    B2_x     :: T
    G_x      :: T
    phi_x    :: T
    S_x      :: T
    Fx_x     :: T
    Fy_x     :: T

    B1_y     :: T
    B2_y     :: T
    G_y      :: T
    phi_y    :: T
    S_y      :: T
    Fx_y     :: T
    Fy_y     :: T

    B1p_x    :: T
    B2p_x    :: T
    Gp_x     :: T
    phip_x   :: T
    Sp_x     :: T
    Fxp_x    :: T
    Fyp_x    :: T

    B1p_y    :: T
    B2p_y    :: T
    Gp_y     :: T
    phip_y   :: T
    Sp_y     :: T
    Fxp_y    :: T
    Fyp_y    :: T

    B1_xx    :: T
    B2_xx    :: T
    G_xx     :: T
    phi_xx   :: T
    S_xx     :: T

    B1_yy    :: T
    B2_yy    :: T
    G_yy     :: T
    phi_yy   :: T
    S_yy     :: T

    B2_xy    :: T
    G_xy     :: T
    phi_xy   :: T
    S_xy     :: T
end

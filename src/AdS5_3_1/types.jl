
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


struct EvolVars{T} <: AbstractVector{T}
    B1  :: Array{T,3}
    B2  :: Array{T,3}
    G   :: Array{T,3}
    phi :: Array{T,3}
    a4  :: Array{T,3}
    fx2 :: Array{T,3}
    fy2 :: Array{T,3}
    xi  :: Array{T,3}
end

"""
    EvolVars{T}(undef, Nu, Nx, Ny)

Construct a container of uninitialized Arrays to hold all the variables that are
evolved in time: B1, B2, G, phi, a4, fx2, fy2, xi.

xi, a4 and fx2, fy2 are automatically defined on a `(1,Nx,Ny)` grid, rather than
a `(Nx,Ny)` one, so that the same Dx and Dy differential operators defined for
the bulk quantities can also straightforwardly apply on them. Remember that the
axis along which the operator applies is defined on the operator itself. So, by
defining things this way, the Dx operator (which acts along the 2nd index) will
also do the correct thing when acting on a4, fx2, fy2 and xi.
"""
function EvolVars{T}(::UndefInitializer, Nu::Int, Nx::Int, Ny::Int) where {T<:Real}
    B1  = Array{T}(undef, Nu, Nx, Ny)
    B2  = Array{T}(undef, Nu, Nx, Ny)
    G   = Array{T}(undef, Nu, Nx, Ny)
    phi = Array{T}(undef, Nu, Nx, Ny)
    a4  = Array{T}(undef,  1, Nx, Ny)
    fx2 = Array{T}(undef,  1, Nx, Ny)
    fy2 = Array{T}(undef,  1, Nx, Ny)
    xi  = Array{T}(undef,  1, Nx, Ny)
    EvolVars{T}(B1, B2, G, phi, a4, fx2, fy2, xi)
end

Base.similar(ff::EvolVars) = EvolVars(similar(ff.B1), similar(ff.B2), similar(ff.G), similar(ff.phi),
                                      similar(ff.a4), similar(ff.fx2), similar(ff.fy2), similar(ff.xi))

Base.length(ff::EvolVars) = length(ff.B1) + length(ff.B2) + length(ff.G) + length(ff.phi) +
    length(ff.a4) + length(ff.fx2) + length(ff.fy2) + length(ff.xi)
Base.size(ff::EvolVars) = (length(ff),)

# indexing. this is just a linear indexing through all the arrays

@inline Base.firstindex(ff::EvolVars) = 1
@inline Base.lastindex(ff::EvolVars) = length(ff)

@inline function Base.getindex(evol::EvolVars, i::Int)
    vars = [:B1, :B2, :G, :phi, :a4, :fx2, :fy2, :xi]
    @inbounds for x in vars
        f  = getproperty(evol,x)   # f will point to B1, then B2, then G, etc...
        i -= length(f)
        if i <= 0
            return f[length(f)+i]
        end
    end
end

@inline function Base.setindex!(evol::EvolVars, v, i::Int)
    vars = [:B1, :B2, :G, :phi, :a4, :fx2, :fy2, :xi]
    @inbounds for x in vars
        f  = getproperty(evol,x)   # f will point to B1, then B2, then G, etc...
        i -= length(f)
        if i <= 0
            f[length(f)+i] = v
            break
        end
    end
end

getB1(evol::EvolVars)  = evol.B1
getB2(evol::EvolVars)  = evol.B2
getG(evol::EvolVars)   = evol.G
getphi(evol::EvolVars) = evol.phi
geta4(evol::EvolVars)  = evol.a4
getfx2(evol::EvolVars) = evol.fx2
getfy2(evol::EvolVars) = evol.fy2
getxi(evol::EvolVars)  = evol.xi



struct BulkVars{GT<:GridType,T}
    gridtype :: GT
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

function BulkVars(gridtype::GT, ::Type{T}, Nxx::Vararg) where{GT<:GridType,T}
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
    BulkVars{GT,typeof(B1)}(gridtype, B1, B2, G, phi, S, Fx, Fy, B1d, B2d, Gd, phid, Sd, A)
end

function BulkVars(gridtype::GT, ff::EvolVars) where {GT}
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
    BulkVars{GT,typeof(B1)}(gridtype, B1, B2, G, phi, S, Fx, Fy, B1d, B2d, Gd, phid, Sd, A)
end

function BulkVars(gridtype::GT, B1::T, B2::T, G::T, phi::T) where {GT<:GridType,T}
    S      = similar(B1)
    Fx     = similar(B1)
    Fy     = similar(B1)
    B1d    = similar(B1)
    B2d    = similar(B1)
    Gd     = similar(B1)
    phid   = similar(B1)
    Sd     = similar(B1)
    A      = similar(B1)

    BulkVars{GT,T}(gridtype, B1, B2, G, phi, S, Fx, Fy, B1d, B2d, Gd, phid, Sd, A)
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

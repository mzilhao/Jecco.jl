module AdS5_3_1

using Jecco
using Parameters

abstract type GridType end
struct Inner <: GridType end
struct Outer <: GridType end

export ParamBase, ParamGrid, ParamID, ParamEvol, ParamIO
export Potential
export VV # this will contain the potential
export Inner, Outer, System
export BulkVars, BoundaryVars, GaugeVars, BaseVars


# TODO: remove d*dt fields from this struct ?

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
    dB1dt    :: T
    dB2dt    :: T
    dGdt     :: T
    dphidt   :: T
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
    dB1dt  = copy(B1)
    dB2dt  = copy(B1)
    dGdt   = copy(B1)
    dphidt = copy(B1)

    BulkVars{GT,typeof(B1)}(gridtype, B1, B2, G, phi, S, Fx, Fy, B1d, B2d, Gd, phid, Sd, A,
                            dB1dt, dB2dt,dGdt, dphidt)
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
    dB1dt  = similar(B1)
    dB2dt  = similar(B1)
    dGdt   = similar(B1)
    dphidt = similar(B1)

    BulkVars{GT,T}(gridtype, B1, B2, G, phi, S, Fx, Fy, B1d, B2d, Gd, phid, Sd, A,
                   dB1dt, dB2dt,dGdt, dphidt)
end


struct GaugeVars{A,T}
    xi    :: A
    kappa :: T
end

function GaugeVars(xi::Array{T,N}, kappa::T) where {T<:Number,N}
    GaugeVars{typeof(xi), typeof(kappa)}(xi, kappa)
end


function setup(par_base)
    global VV = Potential(par_base)
end

struct BaseVars{T}
    phi0  :: T
end


struct BoundaryVars{T}
    a4   :: T
    fx2  :: T
    fy2  :: T
end
BoundaryVars(a4, fx2, fy2) = BoundaryVars{typeof(a4)}(a4, fx2, fy2)


abstract type AbstractVars{GT<:GridType,T<:Real} end


struct SVars{GT<:GridType,T<:Real} <: AbstractVars{GT,T}
    gridtype :: GT

    u        :: T

    phi0     :: T

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

struct FxyVars{GT<:GridType,T<:Real} <: AbstractVars{GT,T}
    gridtype :: GT

    u        :: T

    phi0     :: T

    xi       :: T
    xi_x     :: T
    xi_y     :: T

    B1       :: T
    B1p      :: T
    B1_x     :: T
    B1_y     :: T
    B1pp     :: T
    B1p_x    :: T
    B1p_y    :: T

    B2       :: T
    B2p      :: T
    B2_x     :: T
    B2_y     :: T
    B2pp     :: T
    B2p_x    :: T
    B2p_y    :: T

    G        :: T
    Gp       :: T
    G_x      :: T
    G_y      :: T
    Gpp      :: T
    Gp_x     :: T
    Gp_y     :: T

    phi      :: T
    phip     :: T
    phi_x    :: T
    phi_y    :: T

    S        :: T
    Sp       :: T
    S_x      :: T
    S_y      :: T
    Spp      :: T
    Sp_x     :: T
    Sp_y     :: T
end


struct SdVars{GT<:GridType,T<:Real} <: AbstractVars{GT,T}
    gridtype :: GT

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




#= Notation

for any function f we're using the following notation (let _x denote partial
derivative with respect to x)

fp  = f_r = -u^2 f_u
fd  = \dot f
ft  = \tilde f = f_x - (Fx + xi_x) f_r
fh  = \hat f   = f_y - (Fy + xi_y) f_r

=#


mutable struct AllVars{GT<:GridType,T<:Real}
    gridtype :: GT

    u        :: T

    phi0     :: T

    xi       :: T
    xi_x     :: T
    xi_y     :: T
    xi_xx    :: T
    xi_xy    :: T
    xi_yy    :: T

    B1       :: T
    B1p      :: T
    B1t      :: T
    B1h      :: T
    B1b      :: T
    B1s      :: T
    B1pt     :: T
    B1ph     :: T

    B2       :: T
    B2p      :: T
    B2t      :: T
    B2h      :: T
    B2b      :: T
    B2s      :: T
    B2pt     :: T
    B2ph     :: T

    G        :: T
    Gp       :: T
    Gt       :: T
    Gh       :: T
    Gb       :: T
    Gs       :: T
    Gpt      :: T
    Gph      :: T

    phi      :: T
    phip     :: T
    phit     :: T
    phih     :: T
    phib     :: T
    phis     :: T
    phipt    :: T
    phiph    :: T

    S        :: T
    Sp       :: T
    St       :: T
    Sh       :: T
    Sb       :: T
    Ss       :: T
    Spt      :: T
    Sph      :: T

    Fx       :: T
    Fxp      :: T
    Fxt      :: T
    Fxh      :: T
    Fxb      :: T
    Fxs      :: T
    Fxpt     :: T
    Fxph     :: T

    Fy       :: T
    Fyp      :: T
    Fyt      :: T
    Fyh      :: T
    Fyb      :: T
    Fys      :: T
    Fypt     :: T
    Fyph     :: T

    Sd       :: T
    B1d      :: T
    B2d      :: T
    Gd       :: T
    phid     :: T

    B2c      :: T
    Gc       :: T
    Sc       :: T
    phic     :: T

    Spp      :: T
end
function AllVars(gridtype::GT, ::Type{T}) where {GT<:GridType,T<:Real}
    N = 2 + 6 + 8*7 + 5 + 4 + 1
    array = zeros(T,N)
    AllVars{GT,T}(gridtype, array...)
end


#= Notation

for any function f we're using the following convention: _x denotes partial
derivative with respect to x and

fp  = f_r  = -u^2 f_u
fpp = f_rr = 2u^3 f_u + u^4 f_uu

=#


include("param.jl")
include("system.jl")
# include("initial_data.jl")
include("potential.jl")
# include("dphidt.jl")
include("equation_inner_coeff.jl")
include("equation_outer_coeff.jl")
include("inner_to_outer.jl")
include("solve_nested.jl")
# include("rhs.jl")
# include("run.jl")
# include("ibvp.jl")

end

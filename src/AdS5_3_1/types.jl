
abstract type GridType end
struct Inner <: GridType end
struct Outer <: GridType end


struct EvolVars{T}
    B1       :: T
    B2       :: T
    G        :: T
    phi      :: T
    a4       :: T
    fx2      :: T
    fy2      :: T
    xi       :: T
end

getB1(evol::EvolVars)  = evol.B1
getB2(evol::EvolVars)  = evol.B2
getG(evol::EvolVars)   = evol.G
getphi(evol::EvolVars) = evol.phi
geta4(evol::EvolVars)  = evol.a4
getfx2(evol::EvolVars) = evol.fx2
getfy2(evol::EvolVars) = evol.fy2
getxi(evol::EvolVars)  = evol.xi

getB1(evols::AbstractVector{EvolVars{T}})  where T = VectorOfArray([evol.B1  for evol in evols])
getB2(evols::AbstractVector{EvolVars{T}})  where T = VectorOfArray([evol.B2  for evol in evols])
getG(evols::AbstractVector{EvolVars{T}})   where T = VectorOfArray([evol.G   for evol in evols])
getphi(evols::AbstractVector{EvolVars{T}}) where T = VectorOfArray([evol.phi for evol in evols])
geta4(evols::AbstractVector{EvolVars{T}})  where T = VectorOfArray([evol.a4  for evol in evols])
getfx2(evols::AbstractVector{EvolVars{T}}) where T = VectorOfArray([evol.fx2 for evol in evols])
getfy2(evols::AbstractVector{EvolVars{T}}) where T = VectorOfArray([evol.fy2 for evol in evols])
getxi(evols::AbstractVector{EvolVars{T}})  where T = VectorOfArray([evol.xi  for evol in evols])

pack(B1s, B2s, Gs, phis, a4s, fx2s, fy2s, xis) =
    ArrayPartition(B1s, B2s, Gs, phis, a4s, fx2s, fy2s, xis)

function pack(evols::AbstractVector{EvolVars{T}}) where T
    B1s  = getB1(evols)
    B2s  = getB2(evols)
    Gs   = getG(evols)
    phis = getphi(evols)
    a4s  = geta4(evols)
    fx2s = getfx2(evols)
    fy2s = getfy2(evols)
    xis  = getxi(evols)
    pack(B1s, B2s, Gs, phis, a4s, fx2s, fy2s, xis)
end

getB1(f::ArrayPartition)  = f.x[1]
getB2(f::ArrayPartition)  = f.x[2]
getG(f::ArrayPartition)   = f.x[3]
getphi(f::ArrayPartition) = f.x[4]
geta4(f::ArrayPartition)  = f.x[5]
getfx2(f::ArrayPartition) = f.x[6]
getfy2(f::ArrayPartition) = f.x[7]
getxi(f::ArrayPartition)  = f.x[8]

unpack(f::ArrayPartition) = f.x


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

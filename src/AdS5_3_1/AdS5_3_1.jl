module AdS5_3_1

using Jecco
using Parameters

export ParamBase, ParamGrid, ParamID, ParamEvol, ParamIO
export Potential
export VV # this will contain the potential
export System
export BulkVars, BoundaryVars, AllVars

# Note: in the future we may promote this to something like BulkVars{Ng,T}, to
# dispatch on Ng (the type of equations to be solved on each grid)

struct BulkVars{T}
    B1     :: T
    B2     :: T
    G      :: T
    phi    :: T
    S      :: T
    Fx     :: T
    Fy     :: T
    B1d    :: T
    B2d    :: T
    Gd     :: T
    phid   :: T
    Sd     :: T
    A      :: T
    dB1dt  :: T
    dB2dt  :: T
    dGdt   :: T
    dphidt :: T
end

BulkVars(B1, B2, G, phi, S, Fx, Fy, B1d, B2d, Gd, phid, Sd, A, dB1dt, dB2dt,
         dGdt, dphidt) = BulkVars{typeof(B1)}(B1, B2, G, phi, S, Fx, Fy, B1d, B2d,
                                              Gd, phid, Sd, A, dB1dt, dB2dt, dGdt, dphidt)

function BulkVars(Nxx::Vararg)
    B1     = zeros(Nxx...)
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

    BulkVars{typeof(B1)}(B1, B2, G, phi, S, Fx, Fy, B1d, B2d, Gd, phid, Sd, A,
                         dB1dt, dB2dt,dGdt, dphidt)
end


function BulkVars(B1::Array{T,N}, B2::Array{T,N}, G::Array{T,N},
                  phi::Array{T,N}) where {T<:Number,N}
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

    BulkVars{typeof(B1)}(B1, B2, G, phi, S, Fx, Fy, B1d, B2d, Gd, phid, Sd, A,
                         dB1dt, dB2dt,dGdt, dphidt)
end


# TODO

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


function setup(par_base)
    global VV = Potential(par_base)
end


struct BoundaryVars{T}
    a4   :: T
    fx2  :: T
    fy2  :: T
end


#= Notation

for any function f we're using the following notation (let _x denote partial
derivative with respect to x)

fp  = f_r = -u^2 f_u
fd  = \dot f
ft  = \tilde f = f_x - Fx f_r
fh  = \hat f   = f_y - Fy f_r

=#

mutable struct AllVars{T}
    u        :: T

    B1       :: T
    B1p      :: T
    B1t      :: T
    B1h      :: T
    B1tp     :: T
    B1hp     :: T

    B2       :: T
    B2p      :: T
    B2t      :: T
    B2h      :: T
    B2tp     :: T
    B2hp     :: T

    G        :: T
    Gp       :: T
    Gt       :: T
    Gh       :: T
    Gtp      :: T
    Ghp      :: T

    phi      :: T
    phip     :: T
    phit     :: T
    phih     :: T
    phitp    :: T
    phihp    :: T

    S        :: T
    Sp       :: T
    St       :: T
    Sh       :: T
    Stp      :: T
    Shp      :: T

    Fx       :: T
    Fxp      :: T
    Fxt      :: T
    Fxh      :: T
    Fxtp     :: T
    Fxhp     :: T

    Fy       :: T
    Fyp      :: T
    Fyt      :: T
    Fyh      :: T
    Fytp     :: T
    Fyhp     :: T
end
function AllVars{T}() where {T<:AbstractFloat}
    N = 1 + 6*7
    array = zeros(N)
    AllVars{T}(array...)
end


include("param.jl")
include("system.jl")
# include("initial_data.jl")
include("potential.jl")
# include("dphidt.jl")
include("equation_coeff.jl")
include("solve_nested.jl")
# include("rhs.jl")
# include("run.jl")
# include("ibvp.jl")

end

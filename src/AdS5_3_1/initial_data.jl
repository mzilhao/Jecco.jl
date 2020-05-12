
"""
Extend this type for different initial conditions
"""
abstract type IBVP{T} end

@with_kw struct BlackBrane{T} <: IBVP{T}
    energy_dens   :: T   = 1.0
    AH_pos        :: T   = 1.0
end

# struct EvolVars{GT<:GridType,T}
#     gridtype :: GT
#     B1       :: T
#     B2       :: T
#     G        :: T
#     phi      :: T
#     a4       :: T
#     fx2      :: T
#     fy2      :: T
#     xi       :: T
# end

"""
Needs a field f
"""
abstract type AbstractEvolVars end


struct EvolVars{A} <: AbstractEvolVars
    f :: A
    function EvolVars(B1::T1, B2::T1, G::T1, phi::T1, a4::T2, fx2::T2, fy2::T2, xi::T2) where {T1,T2}
        f = VectorOfArray([B1, B2, G, phi, a4, fx2, fy2, xi])
        new{typeof(f)}(f)
    end
end

getB1(f::EvolVars)   = f.f[1]
getB2(f::EvolVars)   = f.f[2]
getG(f::EvolVars)    = f.f[3]
getphi(f::EvolVars)  = f.f[4]
geta4(f::EvolVars)   = f.f[5]
getfx2(f::EvolVars)  = f.f[6]
getfy2(f::EvolVars)  = f.f[7]
getxi(f::EvolVars)   = f.f[8]


function init!(f::EvolVars, sys::System, ibvp::BlackBrane)
    # Nu, Nx, Ny = size(sys)
    # ucoord = sys.ucoord
    # xcoord = sys.xcoord
    # ycoord = sys.ycoord

    a40     = -ibvp.energy_dens/0.75
    AH_pos  = ibvp.AH_pos

    xi0 = (-a40)^0.25 - 1/AH_pos

    B1  = getB1(f)
    B2  = getB2(f)
    G   = getG(f)
    phi = getphi(f)
    a4  = geta4(f)
    fx2 = getfx2(f)
    fy2 = getfy2(f)
    xi  = getxi(f)

    fill!(B1,  0)
    fill!(B2,  0)
    fill!(G,   0)
    fill!(phi, 0)
    fill!(a4, a40)
    fill!(fx2, 0)
    fill!(fy2, 0)
    fill!(xi, xi0)

    f
end

function init(sys::System, ibvp::IBVP{T}) where{T}
    Nu, Nx, Ny = size(sys)

    B1  = zeros(T, Nu, Nx, Ny)
    B2  = zeros(T, Nu, Nx, Ny)
    G   = zeros(T, Nu, Nx, Ny)
    phi = zeros(T, Nu, Nx, Ny)

    # Note: it's important that xi, a4 and f2 are defined on a 1*Nx*Ny grid,
    # rather than a Nx*Ny one, so that the same Dx and Dy differential operators
    # defined for the bulk quantities can also straightforwardly apply on them.
    # Remember that the axis along which the operator applies is defined on the
    # operator itself. So, by defining things this way, the Dx operator (which
    # acts along the 2nd index) will also do the correct thing when acting on
    # xi, a4 or f2.
    a4  = zeros(T, 1, Nx, Ny)
    fx2 = zeros(T, 1, Nx, Ny)
    fy2 = zeros(T, 1, Nx, Ny)
    xi  = zeros(T, 1, Nx, Ny)

    f = EvolVars(B1, B2, G, phi, a4, fx2, fy2, xi)
    init!(f, sys, ibvp)
end

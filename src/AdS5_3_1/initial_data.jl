
"""
Extend this type for different initial conditions
"""
abstract type IBVP{T} end

@with_kw struct BlackBrane{T} <: IBVP{T}
    energy_dens   :: T   = 1.0
    AH_pos        :: T   = 1.0
end

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


function init!(f::EvolVars, sys::System, ibvp::BlackBrane)
    # Nu, Nx, Ny = size(sys)
    # ucoord = sys.ucoord
    # xcoord = sys.xcoord
    # ycoord = sys.ycoord

    a40     = -ibvp.energy_dens/0.75
    AH_pos  = ibvp.AH_pos

    xi0 = (-a40)^0.25 - 1/AH_pos

    B1  = f.B1
    B2  = f.B2
    G   = f.G
    phi = f.phi
    a4  = f.a4
    fx2 = f.fx2
    fy2 = f.fy2
    xi  = f.xi

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

init(systems::Vector, ibvp::IBVP) = [init(sys, ibvp) for sys in systems]

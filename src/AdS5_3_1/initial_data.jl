
"""
Extend this type for different initial conditions
"""
abstract type IBVP{T} end

@with_kw struct BlackBrane{T} <: IBVP{T}
    energy_dens   :: T   = 1.0
    AH_pos        :: T   = 1.0
end

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

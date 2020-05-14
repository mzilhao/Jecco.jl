
@with_kw struct BlackBrane{T} <: IBVP{T}
    energy_dens   :: T   = 1.0
    AH_pos        :: T   = 1.0
end

function init_data!(f::EvolVars, sys::System, ibvp::BlackBrane)
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

function init_data(sys::System, ibvp::IBVP{T}) where{T}
    Nu, Nx, Ny = size(sys)
    ff = EvolVars{T}(undef, Nu, Nx, Ny)
    init_data!(ff, sys, ibvp)
end

init_data(systems::Vector, ibvp::IBVP) = [init_data(sys, ibvp) for sys in systems]

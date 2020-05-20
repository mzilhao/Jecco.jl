
Base.@kwdef struct BlackBrane{T,TP<:Potential} <: IBVP
    energy_dens   :: T   = 1.0
    AH_pos        :: T   = 1.0
    phi0          :: T   = 0.0
    potential     :: TP  = ZeroPotential()
end

Base.@kwdef struct IDTest0{T,TP<:Potential} <: IBVP
    b14_0     :: T  = 0.0
    b24_0     :: T  = 0.0
    g4_0      :: T  = 0.0
    phi0      :: T  = 0.0
    phi2_0    :: T  = 0.0
    a4_0      :: T  = 0.0
    fx2_0     :: T  = 0.0
    fy2_0     :: T  = 0.0
    xi_0      :: T  = 0.0
    potential :: TP = ZeroPotential()
end


function init_data!(bulks, systems::SystemPartition{Nsys},
                    ibvp::IBVP) where {Nsys,T<:BulkEvol}
    # the Ref() makes its argument a scalar with respect to broadcast
    init_data!.(bulks, systems, Ref(ibvp))
end


# BlackBrane initial data

function init_data!(ff::BulkEvol, sys::System, ibvp::BlackBrane)
    B1  = getB1(ff)
    B2  = getB2(ff)
    G   = getG(ff)
    phi = getphi(ff)

    fill!(B1,  0)
    fill!(B2,  0)
    fill!(G,   0)
    fill!(phi, 0)

    ff
end

function init_data!(ff::Boundary, sys::System, ibvp::BlackBrane)
    a40 = -ibvp.energy_dens/0.75

    a4  = geta4(ff)
    fx2 = getfx2(ff)
    fy2 = getfy2(ff)

    fill!(a4, a40)
    fill!(fx2, 0)
    fill!(fy2, 0)

    ff
end

function init_data!(ff::Gauge, sys::System, ibvp::BlackBrane)
    a40     = -ibvp.energy_dens/0.75
    AH_pos  = ibvp.AH_pos
    xi0     = (-a40)^0.25 - 1/AH_pos

    xi  = getxi(ff)

    fill!(xi, xi0)

    ff
end


# IDTest0

function init_data!(ff::BulkEvol, sys::System{Outer}, ibvp::IDTest0)
    Nu, Nx, Ny = size(sys)
    ucoord = sys.ucoord
    xcoord = sys.xcoord
    ycoord = sys.ycoord

    B1  = getB1(ff)
    B2  = getB2(ff)
    G   = getG(ff)
    phi = getphi(ff)

    b14_0  = ibvp.b14_0
    b24_0  = ibvp.b24_0

    g4_0   = ibvp.g4_0

    phi0   = ibvp.phi0
    phi2_0 = ibvp.phi2_0

    for j in 1:Ny
        for i in 1:Nx
            for a in 1:Nu
                u = ucoord[a]
                x = xcoord[i]
                y = ycoord[j]
                B1[a,i,j]  = u^4 * b14_0
                B2[a,i,j]  = u^4 * b24_0
                phi[a,i,j] = phi0 * u + phi2_0 * u^3
                G[a,i,j]   = u^4 * g4_0
            end
        end
    end

    ff
end

function init_data!(ff::BulkEvol, sys::System{Inner}, ibvp::IDTest0)
    # Nu, Nx, Ny = size(sys)
    # ucoord = sys.ucoord
    # xcoord = sys.xcoord
    # ycoord = sys.ycoord

    B1  = getB1(ff)
    B2  = getB2(ff)
    G   = getG(ff)
    phi = getphi(ff)

    b14_0  = ibvp.b14_0
    b24_0  = ibvp.b24_0

    g4_0   = ibvp.g4_0

    phi0   = ibvp.phi0
    phi2_0 = ibvp.phi2_0

    fill!(B1,  b14_0)
    fill!(B2,  b24_0)
    fill!(G,   g4_0)
    fill!(phi, phi2_0)

    ff
end

function init_data!(ff::Boundary, sys::System{Inner}, ibvp::IDTest0)
    # _, Nx, Ny = size(sys)
    # xcoord = sys.xcoord
    # ycoord = sys.ycoord

    a4  = geta4(ff)
    fx2 = getfx2(ff)
    fy2 = getfy2(ff)

    a4_0   = ibvp.a4_0
    fx2_0  = ibvp.fx2_0
    fy2_0  = ibvp.fy2_0

    fill!(a4,  a4_0)
    fill!(fx2, fx2_0)
    fill!(fy2, fy2_0)

    ff
end

function init_data!(ff::Gauge, sys::System{Outer}, ibvp::IDTest0)
    # _, Nx, Ny = size(sys)
    # xcoord = sys.xcoord
    # ycoord = sys.ycoord

    xi   = getxi(ff)
    xi_0 = ibvp.xi_0
    fill!(xi,  xi_0)

    ff
end

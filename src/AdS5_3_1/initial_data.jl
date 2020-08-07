
Base.@kwdef struct IDTest0{T,TP<:Potential} <: InitialData
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

Base.@kwdef struct BlackBrane{T,TP<:Potential} <: InitialData
    energy_dens   :: T   = 1.0
    AH_pos        :: T   = 1.0
    phi0          :: T   = 0.0
    potential     :: TP  = ZeroPotential()
end

Base.@kwdef struct BlackBraneB1Pert{T,TP<:Potential} <: InitialData
    energy_dens   :: T   = 1.0
    AH_pos        :: T   = 1.0
    phi0          :: T   = 0.0
    potential     :: TP  = ZeroPotential()
    amp           :: T   = 1.e-1
end

Base.@kwdef struct BlackBraneB1_a4_pert{T,TP<:Potential} <: InitialData
    energy_dens   :: T   = 1.0
    AH_pos        :: T   = 1.0
    phi0          :: T   = 0.0
    potential     :: TP  = ZeroPotential()
    B1_amp        :: T   = 0.0
    a4_amp        :: T   = 5.e-2
    a4_k          :: Int = 1
end

Base.@kwdef struct BlackBranePert_B1B2G{T,TP<:Potential} <: InitialData
    energy_dens   :: T   = 1.0
    AH_pos        :: T   = 1.0
    phi0          :: T   = 0.0
    potential     :: TP  = ZeroPotential()
    amp           :: T   = 1.e-1
end

Base.@kwdef struct PhiGaussian_u{T,TP<:Potential} <: InitialData
    energy_dens   :: T   = 1.0
    AH_pos        :: T   = 1.0
    phi0          :: T   = 0.0
    phi2          :: T   = 0.0
    amp           :: T   = 0.0
    u0            :: T   = 0.0
    sigma         :: T   = 0.1
    potential     :: TP
end

function init_data!(bulkevols, boundary::Boundary, gauge::Gauge,
                    systems::SystemPartition, id::InitialData)

    init_data!(bulkevols, systems, id)
    init_data!(boundary, systems[1],   id)
    init_data!(gauge,    systems[end], id)

    nothing
end


function init_data!(bulkevols, systems::SystemPartition{Nsys},
                    id::InitialData) where {Nsys,T<:BulkEvolved}
    # the Ref() makes its argument a scalar with respect to broadcast
    init_data!.(bulkevols, systems, Ref(id))
end


# BlackBrane initial data

function init_data!(ff::BulkEvolved, sys::System, id::BlackBrane)
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

function init_data!(ff::Boundary, sys::System, id::BlackBrane)
    a40 = -id.energy_dens/0.75

    a4  = geta4(ff)
    fx2 = getfx2(ff)
    fy2 = getfy2(ff)

    fill!(a4, a40)
    fill!(fx2, 0)
    fill!(fy2, 0)

    ff
end

function init_data!(ff::Gauge, sys::System, id::BlackBrane)
    a40     = -id.energy_dens/0.75
    AH_pos  = id.AH_pos
    xi0     = (-a40)^0.25 - 1/AH_pos

    xi  = getxi(ff)

    fill!(xi, xi0)

    ff
end

# BlackBraneB1Pert initial data

function init_data!(ff::BulkEvolved, sys::System{Inner}, id::BlackBraneB1Pert)
    B1  = getB1(ff)
    B2  = getB2(ff)
    G   = getG(ff)
    phi = getphi(ff)

    fill!(B2,  0)
    fill!(G,   0)
    fill!(phi, 0)

    # add the perturbation on B1 id
    pert_amp = id.amp
    # number of maxima in each direction
    nx   = 1
    ny   = 2

    Nu, Nx, Ny = size(sys)
    uu = sys.ucoord
    xx = sys.xcoord
    xmin = xx[1]
    dx   = xx[2] - xx[1] 
    xmax = xx[end] + dx
    yy   = sys.ycoord
    ymin = yy[1]
    dy   = yy[2] - yy[1]  
    ymax = yy[end] + dy
    
    for j in 1:Ny
        for i in 1:Nx
            for a in 1:Nu
                x = xx[i]
                y = yy[j]
                B1[a,i,j]  = pert_amp * sin(2 * π * nx * (xmax-x)/(xmax-xmin) ) * sin(-2 * π * ny * (ymax-y)/(ymax-ymin) )
            end
        end
    end
        
    ff
end

function init_data!(ff::BulkEvolved, sys::System{Outer}, id::BlackBraneB1Pert)
    B1  = getB1(ff)
    B2  = getB2(ff)
    G   = getG(ff)
    phi = getphi(ff)

    fill!(B2,  0)
    fill!(G,   0)
    fill!(phi, 0)

    # add the perturbation on B1 id
    pert_amp = id.amp
    # number of maxima in each direction
    nx   = 1
    ny   = 2

    Nu, Nx, Ny = size(sys)
    uu = sys.ucoord
    xx = sys.xcoord
    xmin = xx[1]
    dx   = xx[2] - xx[1] 
    xmax = xx[end] + dx
    yy   = sys.ycoord
    ymin = yy[1]
    dy   = yy[2] - yy[1]  
    ymax = yy[end] + dy

    for j in 1:Ny
        for i in 1:Nx
            for a in 1:Nu
                x = xx[i]
                y = yy[j]
                u = uu[a]
                B1[a,i,j]  = u^4 * pert_amp * sin(2 * π * nx * (xmax-x)/(xmax-xmin) ) * sin(-2 * π * ny * (ymax-y)/(ymax-ymin) )
            end
        end
    end
        
    ff
end

function init_data!(ff::Boundary, sys::System, id::BlackBraneB1Pert)
    a40 = -id.energy_dens/0.75

    a4  = geta4(ff)
    fx2 = getfx2(ff)
    fy2 = getfy2(ff)

    fill!(a4, a40)
    fill!(fx2, 0)
    fill!(fy2, 0)

    ff
end

function init_data!(ff::Gauge, sys::System, id::BlackBraneB1Pert)
    a40     = -id.energy_dens/0.75
    AH_pos  = id.AH_pos
    xi0     = (-a40)^0.25 - 1/AH_pos

    xi  = getxi(ff)

    fill!(xi, xi0)

    ff
end


# BlackBraneB1_a4_pert initial data

function init_data!(ff::BulkEvolved, sys::System{Inner}, id::BlackBraneB1_a4_pert)
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

function init_data!(ff::BulkEvolved, sys::System{Outer}, id::BlackBraneB1_a4_pert)
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

function init_data!(ff::Boundary, sys::System{Inner}, id::BlackBraneB1_a4_pert)

    # add the perturbation on a4 id
    # amplitude
    amp = id.a4_amp
    # number of maxima
    kx  = id.a4_k

    # get coordinates
    _, Nx, Ny = size(sys)
    xx = sys.xcoord
    yy = sys.ycoord

    xmin = xx[1]
    dx   = xx[2] - xx[1] 
    xmax = xx[end] + dx
    xmid = 0.0#xx[Nx/2 + 1]
    #yy   = sys.ycoord
    #ymin = yy[1]
    #dy   = yy[2] - yy[1]  
    #ymax = yy[end] + dy
    
    a40 = -id.energy_dens/0.75

    a4  = geta4(ff)
    fx2 = getfx2(ff)
    fy2 = getfy2(ff)

    fill!(a4, a40)
    fill!(fx2, 0)
    fill!(fy2, 0)
    
    for j in 1:Ny
        for i in 1:Nx
            x = xx[i]
            a4[1,i,j]  += -a40 * amp * cos(2 * π * kx * (x-xmid)/(xmax-xmin) )
        end
    end

    ff
end

function init_data!(ff::Gauge, sys::System, id::BlackBraneB1_a4_pert)
    a40     = -id.energy_dens/0.75
    AH_pos  = id.AH_pos
    xi0     = (-a40)^0.25 - 1/AH_pos

    xi  = getxi(ff)

    fill!(xi, xi0)

    ff
end

# BlackBranePert_B1B2G initial data

function init_data!(ff::BulkEvolved, sys::System{Inner}, id::BlackBranePert_B1B2G)
    B1  = getB1(ff)
    B2  = getB2(ff)
    G   = getG(ff)
    phi = getphi(ff)

    fill!(phi, 0)

    # add the perturbation on B1 id
    pert_amp = id.amp
    # number of maxima in each direction
    nx   = 1
    ny   = 2

    Nu, Nx, Ny = size(sys)
    uu = sys.ucoord
    xx = sys.xcoord
    xmin = xx[1]
    dx   = xx[2] - xx[1] 
    xmax = xx[end] + dx
    yy   = sys.ycoord
    ymin = yy[1]
    dy   = yy[2] - yy[1]  
    ymax = yy[end] + dy
    
    for j in 1:Ny
        for i in 1:Nx
            for a in 1:Nu
                x = xx[i]
                y = yy[j]
                B1[a,i,j]  = pert_amp * sin(2 * π * nx * (xmax-x)/(xmax-xmin) ) * sin(-2 * π * ny * (ymax-y)/(ymax-ymin) )
            end
        end
    end

    for j in 1:Ny
        for i in 1:Nx
            for a in 1:Nu
                x = xx[i]
                y = yy[j]
                B2[a,i,j]  = pert_amp * sin(2 * π * nx * (xmax-x)/(xmax-xmin) ) * sin(-2 * π * ny * (ymax-y)/(ymax-ymin) )
            end
        end
    end

    for j in 1:Ny
        for i in 1:Nx
            for a in 1:Nu
                x = xx[i]
                y = yy[j]
                G[a,i,j]  = pert_amp * sin(2 * π * nx * (xmax-x)/(xmax-xmin) ) * sin(-2 * π * ny * (ymax-y)/(ymax-ymin) )
            end
        end
    end
    
    ff
end

function init_data!(ff::BulkEvolved, sys::System{Outer}, id::BlackBranePert_B1B2G)
    B1  = getB1(ff)
    B2  = getB2(ff)
    G   = getG(ff)
    phi = getphi(ff)
    
    fill!(phi, 0)
    
    # add the perturbation on B1 id
    pert_amp = id.amp
    # number of maxima in each direction
    nx   = 1
    ny   = 2

    Nu, Nx, Ny = size(sys)
    uu = sys.ucoord
    xx = sys.xcoord
    xmin = xx[1]
    dx   = xx[2] - xx[1] 
    xmax = xx[end] + dx
    yy   = sys.ycoord
    ymin = yy[1]
    dy   = yy[2] - yy[1]  
    ymax = yy[end] + dy

    for j in 1:Ny
        for i in 1:Nx
            for a in 1:Nu
                x = xx[i]
                y = yy[j]
                u = uu[a]
                B1[a,i,j]  = u^4 * pert_amp * sin(2 * π * nx * (xmax-x)/(xmax-xmin) ) * sin(-2 * π * ny * (ymax-y)/(ymax-ymin) )
            end
        end
    end

    
    for j in 1:Ny
        for i in 1:Nx
            for a in 1:Nu
                x = xx[i]
                y = yy[j]
                u = uu[a]
                B2[a,i,j]  = u^4 * pert_amp * sin(2 * π * nx * (xmax-x)/(xmax-xmin) ) * sin(-2 * π * ny * (ymax-y)/(ymax-ymin) )
            end
        end
    end

    
    for j in 1:Ny
        for i in 1:Nx
            for a in 1:Nu
                x = xx[i]
                y = yy[j]
                u = uu[a]
                G[a,i,j]  = u^4 * pert_amp * sin(2 * π * nx * (xmax-x)/(xmax-xmin) ) * sin(-2 * π * ny * (ymax-y)/(ymax-ymin) )
            end
        end
    end

    ff
end

function init_data!(ff::Boundary, sys::System, id::BlackBranePert_B1B2G)
    a40 = -id.energy_dens/0.75

    a4  = geta4(ff)
    fx2 = getfx2(ff)
    fy2 = getfy2(ff)

    fill!(a4, a40)
    fill!(fx2, 0)
    fill!(fy2, 0)

    ff
end

function init_data!(ff::Gauge, sys::System, id::BlackBranePert_B1B2G)
    a40     = -id.energy_dens/0.75
    AH_pos  = id.AH_pos
    xi0     = (-a40)^0.25 - 1/AH_pos

    xi  = getxi(ff)

    fill!(xi, xi0)

    ff
end


# PhiGaussian_u initial data

function analytic_phi(u, x, y, id::PhiGaussian_u)
    phi0   = id.phi0
    phi2   = id.phi2
    amp    = id.amp
    sigma  = id.sigma
    u0     = id.u0

    phi03  = phi0 * phi0 * phi0
    sigma2 = sigma * sigma

    myfunc = amp * exp( -(u-u0)*(u-u0) / (2*sigma2) )

    (phi2 + u * myfunc) / phi03
end

analytic_B1(u, x, y, id::PhiGaussian_u) = 0
analytic_B2(u, x, y, id::PhiGaussian_u) = 0
analytic_G(u, x, y, id::PhiGaussian_u) = 0

function init_data!(bulkevols, boundary::Boundary, gauge::Gauge,
                    systems::SystemPartition, id::PhiGaussian_u)

    init_data!(boundary, systems[1],   id)
    init_data!(gauge,    systems[end], id)
    init_data!(bulkevols, gauge, systems, id)

    nothing
end

function init_data!(bulkevols, gauge::Gauge, systems::SystemPartition,
                    id::PhiGaussian_u) where {Nsys,T<:BulkEvolved}
    # the Ref() makes its argument a scalar with respect to broadcast
    init_data!.(bulkevols, Ref(gauge), systems, Ref(id))
end


function init_data!(ff::Boundary, sys::System, id::PhiGaussian_u)
    a4  = geta4(ff)
    fx2 = getfx2(ff)
    fy2 = getfy2(ff)

    epsilon = id.energy_dens
    phi0    = id.phi0
    phi2    = id.phi2
    oophiM2 = id.potential.oophiM2
    phi04   = phi0 * phi0 * phi0 * phi0

    # TODO: does this change with nonzero phiQ ?
    a40 = (-epsilon - phi0 * phi2 - phi04 * oophiM2 / 4 + 7 * phi04 / 36) / 0.75

    fill!(a4, a40)
    fill!(fx2, 0)
    fill!(fy2, 0)

    ff
end

function init_data!(ff::Gauge, sys::System, id::PhiGaussian_u)
    epsilon = id.energy_dens
    phi0    = id.phi0
    phi2    = id.phi2
    AH_pos  = id.AH_pos
    oophiM2 = id.potential.oophiM2
    phi04   = phi0 * phi0 * phi0 * phi0

    # TODO: does this change with nonzero phiQ ?
    a40 = (-epsilon - phi0 * phi2 - phi04 * oophiM2 / 4 + 7 * phi04 / 36) / 0.75

    # FIXME: this is only valid for the conformal case! this routine needs to be
    # fixed once we can search for the AH
    xi0 = (-a40)^0.25 - 1/AH_pos

    xi  = getxi(ff)

    fill!(xi, xi0)

    ff
end

function init_data!(bulk::BulkEvolved, gauge::Gauge, sys::System{Inner},
                    id::PhiGaussian_u)
    Nu, Nx, Ny = size(sys)
    xx = sys.xcoord
    yy = sys.ycoord
    uu = sys.ucoord

    B1  = getB1(bulk)
    B2  = getB2(bulk)
    G   = getG(bulk)
    phi = getphi(bulk)
    xi  = getxi(gauge)

    phi0 = id.phi0

    for j in 1:Ny
        for i in 1:Nx
            for a in 1:Nu
                u       = uu[a]
                x       = xx[i]
                y       = yy[j]
                xi_ij   = xi[1,i,j]
                aux     = 1 + xi_ij * u
                aux3    = aux * aux * aux
                aux4    = aux * aux3
                u_old   = u / aux
                B1_old  = analytic_B1(u_old, x, y, id)
                B2_old  = analytic_B2(u_old, x, y, id)
                G_old   = analytic_G(u_old, x, y, id)
                phi_old = analytic_phi(u_old, x, y, id)

                B1[a,i,j]  = B1_old / aux4
                B2[a,i,j]  = B2_old / aux4
                G[a,i,j]   = G_old  / aux4
                phi[a,i,j] = xi_ij * xi_ij / (phi0 * phi0 * aux) +
                    phi_old / aux3
            end
        end
    end

    bulk
end

function init_data!(bulk::BulkEvolved, gauge::Gauge, sys::System{Outer},
                    id::PhiGaussian_u)
    Nu, Nx, Ny = size(sys)
    xx = sys.xcoord
    yy = sys.ycoord
    uu = sys.ucoord

    B1  = getB1(bulk)
    B2  = getB2(bulk)
    G   = getG(bulk)
    phi = getphi(bulk)
    xi  = getxi(gauge)

    phi0 = id.phi0

    for j in 1:Ny
        for i in 1:Nx
            for a in 1:Nu
                u         = uu[a]
                x         = xx[i]
                y         = yy[j]
                xi_ij     = xi[1,i,j]
                aux       = 1 + xi_ij * u
                aux3      = aux * aux * aux
                aux4      = aux * aux3
                u_old     = u / aux
                B1_old    = analytic_B1(u_old, x, y, id)
                B2_old    = analytic_B2(u_old, x, y, id)
                G_old     = analytic_G(u_old, x, y, id)
                phi_old   = analytic_phi(u_old, x, y, id)
                B1_inner  = B1_old / aux4
                B2_inner  = B2_old / aux4
                G_inner   = G_old  / aux4
                phi_inner = xi_ij * xi_ij / (phi0 * phi0 * aux) +
                    phi_old / aux3

                B1[a,i,j]  = u^4 * B1_inner
                B2[a,i,j]  = u^4 * B2_inner
                G[a,i,j]   = u^4 * G_inner
                phi[a,i,j] = u * phi0 - u^2 * phi0 * xi_ij + u^3 * phi0^3 * phi_inner
            end
        end
    end

    bulk
end


# IDTest0

function init_data!(ff::BulkEvolved, sys::System{Outer}, id::IDTest0)
    Nu, Nx, Ny = size(sys)
    ucoord = sys.ucoord
    xcoord = sys.xcoord
    ycoord = sys.ycoord

    B1  = getB1(ff)
    B2  = getB2(ff)
    G   = getG(ff)
    phi = getphi(ff)

    b14_0  = id.b14_0
    b24_0  = id.b24_0

    g4_0   = id.g4_0

    phi0   = id.phi0
    phi2_0 = id.phi2_0

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

function init_data!(ff::BulkEvolved, sys::System{Inner}, id::IDTest0)
    # Nu, Nx, Ny = size(sys)
    # ucoord = sys.ucoord
    # xcoord = sys.xcoord
    # ycoord = sys.ycoord

    B1  = getB1(ff)
    B2  = getB2(ff)
    G   = getG(ff)
    phi = getphi(ff)

    b14_0  = id.b14_0
    b24_0  = id.b24_0

    g4_0   = id.g4_0

    phi0   = id.phi0
    phi2_0 = id.phi2_0

    fill!(B1,  b14_0)
    fill!(B2,  b24_0)
    fill!(G,   g4_0)
    fill!(phi, phi2_0)

    ff
end

function init_data!(ff::Boundary, sys::System{Inner}, id::IDTest0)
    # _, Nx, Ny = size(sys)
    # xcoord = sys.xcoord
    # ycoord = sys.ycoord

    a4  = geta4(ff)
    fx2 = getfx2(ff)
    fy2 = getfy2(ff)

    a4_0   = id.a4_0
    fx2_0  = id.fx2_0
    fy2_0  = id.fy2_0

    fill!(a4,  a4_0)
    fill!(fx2, fx2_0)
    fill!(fy2, fy2_0)

    ff
end

function init_data!(ff::Gauge, sys::System{Outer}, id::IDTest0)
    # _, Nx, Ny = size(sys)
    # xcoord = sys.xcoord
    # ycoord = sys.ycoord

    xi   = getxi(ff)
    xi_0 = id.xi_0
    fill!(xi,  xi_0)

    ff
end

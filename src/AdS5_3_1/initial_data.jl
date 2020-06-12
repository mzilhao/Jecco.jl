
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

Base.@kwdef struct BlackBranePert_B1B2G{T,TP<:Potential} <: InitialData
    energy_dens   :: T   = 1.0
    AH_pos        :: T   = 1.0
    phi0          :: T   = 0.0
    potential     :: TP  = ZeroPotential()
    amp           :: T   = 1.e-1
end

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

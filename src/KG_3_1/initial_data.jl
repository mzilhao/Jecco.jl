
function initial_data(sys, p::ParamID)
    if p.ID_type == "sine2D"
        return sine2D(sys, p)
    elseif p.ID_type == "uniform2D"
        return uniform2D(sys, p)
    else
        error("Unknown initial data type.")
    end
end

function uniform2D(sys::System, p::ParamID)
    ucoord = sys.ucoord
    xcoord = sys.xcoord
    ycoord = sys.ycoord

    Nu = length(ucoord)
    Nx = length(xcoord)
    Ny = length(ycoord)
    phif  = zeros(Nu, Nx, Ny)

    # TODO: make parameter
    phi2 = 1.0

    for j in 1:Ny
        for i in 1:Nx
            for a in 1:Nu
                phif[a,i,j] = phi2
            end
        end
    end

    phif
end

uniform2D(systems::Array, p::ParamID) = [uniform2D(sys, p) for sys in systems]

sine2D(x, y, Lx::Real, Ly::Real, kx::Integer, ky::Integer) =
             sin( 2*π * kx / Lx * x ) * sin( 2*π * ky / Ly * y )

function sine2D(sys::System, p::ParamID)
    ucoord = sys.ucoord
    xcoord = sys.xcoord
    ycoord = sys.ycoord

    Nu = length(ucoord)
    Nx = length(xcoord)
    Ny = length(ycoord)
    phif  = zeros(Nu, Nx, Ny)

    Lx    = p.Lx
    Ly    = p.Ly

    kx = 2
    ky = 4

    for j in 1:Ny
        for i in 1:Nx
            for a in 1:Nu
                x = xcoord[i]
                y = ycoord[j]
                phif[a,i,j] = sine2D(x, y, Lx, Ly, kx, ky)
            end
        end
    end

    phif
end

sine2D(systems::Array, p::ParamID) = [sine2D(sys, p) for sys in systems]

function ones2D(sys::System)
    ucoord = sys.ucoord
    xcoord = sys.xcoord
    ycoord = sys.ycoord

    Nu = length(ucoord)
    Nx = length(xcoord)
    Ny = length(ycoord)

    ones(Nx, Ny)
end

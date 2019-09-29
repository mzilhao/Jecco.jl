
sine2D(x, y, Lx::Real, Ly::Real, kx::Integer, ky::Integer) =
             sin( 2*π * kx / Lx * x ) * sin( 2*π * ky / Ly * y )

function sine2D(sys::System, p::Param)
    coords = sys.coords

    uu, xx, yy = Vivi.xx(coords)
    Nu = length(uu)
    Nx = length(xx)
    Ny = length(yy)
    phif  = zeros(Nu, Nx, Ny)

    Lx    = p.xmax - p.xmin
    Ly    = p.ymax - p.ymin

    kx = 2
    ky = 4

    for j in eachindex(yy)
        for i in eachindex(xx)
            for a in eachindex(uu)
                x = xx[i]
                y = yy[j]
                phif[a,i,j] = sine2D(x, y, Lx, Ly, kx, ky)
            end
        end
    end

    phif
end

function ones2D(sys::System)
    coords = sys.coords

    uu, xx, yy = Vivi.xx(coords)
    Nu = length(uu)
    Nx = length(xx)
    Ny = length(yy)

    ones(Nx, Ny)
end


function cosine2D(sys::System, p::Param)
    coords = sys.coords
    derivs = sys.derivs

    uu, xx, yy = Vivi.xx(coords)
    Nu = length(uu)
    Nx = length(xx)
    Ny = length(yy)
    phif  = zeros(Nu, Nx, Ny)

    Lx    = xx[end] - xx[1]
    Ly    = yy[end] - yy[1]
    xmid  = 0.5 * (xx[1] + xx[end])
    ymid  = 0.5 * (yy[1] + yy[end])

    kx = 2.0
    ky = 4.0

    for j in eachindex(yy)
        for i in eachindex(xx)
            for a in eachindex(uu)
                phif[a,i,j] = cos( 2*π * kx / Lx * (xx[i] - xmid) ) *
                    cos( 2*π * ky / Ly * (yy[j] - ymid) )
            end
        end
    end

    phif
end

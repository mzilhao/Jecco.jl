
@with_kw struct ParamGrid
    xmin        :: Float64
    xmax        :: Float64
    xnodes      :: Int
    ymin        :: Float64
    ymax        :: Float64
    ynodes      :: Int
    umin        :: Float64
    umax        :: Float64
    udomains    :: Int     = 1
    unodes      :: Int # number of points per domain
end

struct System{C,D,E} <: Vivi.System
    coords :: C
    uderiv :: D
    xderiv :: E
    yderiv :: E
end

function System(coords::CoordSystem)
    # FIXME
    ord    = 4
    BC     = :periodic

    # dx     = xcoord.delta
    # dy     = ycoord.delta
    # dt0    = p.dtfac * min(dx, dy)

    derivs = Vivi.Deriv(coords, (nothing, ord, ord), (nothing, BC, BC))
    uderiv = derivs[1]
    xderiv = derivs[2]
    yderiv = derivs[3]

    System{typeof(coords), typeof(uderiv), typeof(xderiv)}(coords, uderiv, xderiv, yderiv)
end

function System(ucoord::SpectralCoord, xcoord::CartCoord, ycoord::CartCoord)
    coords = Vivi.CoordSystem{Float64}("uxy", [ucoord, xcoord, ycoord])
    System(coords)
end

function create_sys(p::ParamGrid)
    Nsys = p.udomains
    delta_udom = (p.umax - p.umin) / Nsys

    ucoord = [Vivi.SpectralCoord("u", p.umin + (i-1)*delta_udom, p.umin + i*delta_udom,
                                 p.unodes) for i in 1:Nsys]
    xcoord  = Vivi.CartCoord("x", p.xmin, p.xmax, p.xnodes, endpoint=false)
    ycoord  = Vivi.CartCoord("y", p.ymin, p.ymax, p.ynodes, endpoint=false)

    systems = [System(ucoord[i], xcoord, ycoord) for i in 1:Nsys]

    systems
end

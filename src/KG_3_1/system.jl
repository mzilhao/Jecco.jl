
# TODO:
#  - identify instances where "coords" is used. see what the need for it is, ..

struct System{C,Du,Dx,Dy} <: Vivi.System
    coords :: C
    Du     :: Du
    Duu    :: Du
    Dx     :: Dx
    Dxx    :: Dx
    Dy     :: Dy
    Dyy    :: Dy
end

function System(ucoord::SpectralCoord, xcoord::CartCoord, ycoord::CartCoord)
    coords = Vivi.CoordSystem{Float64}("uxy", [ucoord, xcoord, ycoord])

    ord = 4

    Du  = ChebDeriv{1}(1, ucoord.min, ucoord.max, Int64(ucoord.nodes))
    Duu = ChebDeriv{1}(2, ucoord.min, ucoord.max, Int64(ucoord.nodes))

    Dx  = CenteredDiff{2}(1, ord, xcoord.delta, Int64(xcoord.nodes))
    Dxx = CenteredDiff{2}(2, ord, xcoord.delta, Int64(xcoord.nodes))

    Dy  = CenteredDiff{3}(1, ord, ycoord.delta, Int64(ycoord.nodes))
    Dyy = CenteredDiff{3}(2, ord, ycoord.delta, Int64(ycoord.nodes))

    System{typeof(coords), typeof(Du), typeof(Dx),
           typeof(Dy)}(coords, Du, Duu, Dx, Dxx, Dy, Dyy)
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

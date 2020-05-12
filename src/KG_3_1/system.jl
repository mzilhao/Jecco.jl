
struct System{Cu,Cx,Cy,Du,Dx,Dy}
    ucoord :: Cu
    xcoord :: Cx
    ycoord :: Cy
    Du     :: Du
    Duu    :: Du
    Dx     :: Dx
    Dxx    :: Dx
    Dy     :: Dy
    Dyy    :: Dy
end

function System(ucoord::GaussLobattoCoord,
                xcoord::CartesianCoord,
                ycoord::CartesianCoord)
    ord = 4

    Du  = ChebDeriv{1}(1, ucoord.min, ucoord.max, ucoord.nodes)
    Duu = ChebDeriv{1}(2, ucoord.min, ucoord.max, ucoord.nodes)

    Dx  = CenteredDiff{2}(1, ord, Jecco.delta(xcoord), xcoord.nodes)
    Dxx = CenteredDiff{2}(2, ord, Jecco.delta(xcoord), xcoord.nodes)

    Dy  = CenteredDiff{3}(1, ord, Jecco.delta(ycoord), ycoord.nodes)
    Dyy = CenteredDiff{3}(2, ord, Jecco.delta(ycoord), ycoord.nodes)

    System{typeof(ucoord), typeof(xcoord), typeof(ycoord), typeof(Du), typeof(Dx),
           typeof(Dy)}(ucoord, xcoord, ycoord, Du, Duu, Dx, Dxx, Dy, Dyy)
end

size(sys::System) = (sys.ucoord.nodes, sys.xcoord.nodes, sys.ycoord.nodes)

function create_systems(p::ParamGrid)
    Nsys = p.udomains
    delta_udom = (p.umax - p.umin) / Nsys

    ucoords = [GaussLobatto{1}("u", p.umin + (i-1)*delta_udom, p.umin + i*delta_udom,
                               p.unodes) for i in 1:Nsys]
    xcoord  = Cartesian{2}("x", p.xmin, p.xmax, p.xnodes, endpoint=false)
    ycoord  = Cartesian{3}("y", p.ymin, p.ymax, p.ynodes, endpoint=false)

    systems = [System(ucoords[i], xcoord, ycoord) for i in 1:Nsys]

    systems
end

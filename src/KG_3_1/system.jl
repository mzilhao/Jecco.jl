
struct System{G,Du,Dx,Dy}
    grid   :: G
    Du     :: Du
    Duu    :: Du
    Dx     :: Dx
    Dxx    :: Dx
    Dy     :: Dy
    Dyy    :: Dy
end

function System(ucoord::AbstractCoord{T,1,GaussLobatto},
                xcoord::AbstractCoord{T,2,Cartesian},
                ycoord::AbstractCoord{T,3,Cartesian}) where {T<:Real}
    grid = Grid(ucoord, xcoord, ycoord)

    ord = 4

    Du  = ChebDeriv{1}(1, ucoord.min, ucoord.max, ucoord.nodes)
    Duu = ChebDeriv{1}(2, ucoord.min, ucoord.max, ucoord.nodes)

    Dx  = CenteredDiff{2}(1, ord, delta(xcoord), xcoord.nodes)
    Dxx = CenteredDiff{2}(2, ord, delta(xcoord), xcoord.nodes)

    Dy  = CenteredDiff{3}(1, ord, delta(ycoord), ycoord.nodes)
    Dyy = CenteredDiff{3}(2, ord, delta(ycoord), ycoord.nodes)

    System{typeof(grid), typeof(Du), typeof(Dx),
           typeof(Dy)}(grid, Du, Duu, Dx, Dxx, Dy, Dyy)
end

function create_systems(p::ParamGrid)
    Nsys = p.udomains
    delta_udom = (p.umax - p.umin) / Nsys

    ucoords = [SpectralCoord{1}("u", p.umin + (i-1)*delta_udom, p.umin + i*delta_udom,
                               p.unodes) for i in 1:Nsys]
    xcoord  = CartCoord{2}("x", p.xmin, p.xmax, p.xnodes, endpoint=false)
    ycoord  = CartCoord{3}("y", p.ymin, p.ymax, p.ynodes, endpoint=false)

    systems = [System(ucoords[i], xcoord, ycoord) for i in 1:Nsys]

    systems
end

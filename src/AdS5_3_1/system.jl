
import Base: size

abstract type GridType end
abstract type Inner <: GridType end
abstract type Outer <: GridType end

abstract type AbstractSystem{GT} end

grid_type(sys::AbstractSystem{GT}) where{GT<:GridType} = GT

struct System{GT,Cu,Cx,Cy,Du,Dx,Dy} <: AbstractSystem{GT}
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

function System{GT}(ucoord::GaussLobattoCoord,
                    xcoord::CartesianCoord,
                    ycoord::CartesianCoord) where {T<:Real,GT<:GridType}
    ord = 4

    Du  = ChebDeriv{1}(1, ucoord.min, ucoord.max, ucoord.nodes)
    Duu = ChebDeriv{1}(2, ucoord.min, ucoord.max, ucoord.nodes)

    Dx  = CenteredDiff{2}(1, ord, Jecco.delta(xcoord), xcoord.nodes)
    Dxx = CenteredDiff{2}(2, ord, Jecco.delta(xcoord), xcoord.nodes)

    Dy  = CenteredDiff{3}(1, ord, Jecco.delta(ycoord), ycoord.nodes)
    Dyy = CenteredDiff{3}(2, ord, Jecco.delta(ycoord), ycoord.nodes)

    System{GT,typeof(ucoord), typeof(xcoord), typeof(ycoord), typeof(Du), typeof(Dx),
           typeof(Dy)}(ucoord, xcoord, ycoord, Du, Duu, Dx, Dxx, Dy, Dyy)
end

size(sys::System) = (sys.ucoord.nodes, sys.xcoord.nodes, sys.ycoord.nodes)

function create_systems(p::ParamGrid)
    u_inner_coord = GaussLobatto{1}("u", 0.0, p.u_outer_min, p.u_inner_nodes)

    N_outer_sys = p.u_outer_domains
    delta_udom  = (p.u_outer_max - p.u_outer_min) / N_outer_sys

    u_outer_coords =
        [GaussLobatto{1}("u", p.u_outer_min + (i-1)*delta_udom,
                          p.u_outer_min + i*delta_udom, p.u_outer_nodes)
         for i in 1:N_outer_sys]

    xcoord  = Cartesian{2}("x", p.x_min, p.x_max, p.x_nodes, endpoint=false)
    ycoord  = Cartesian{3}("y", p.y_min, p.y_max, p.y_nodes, endpoint=false)

    inner_system = System{Inner}(u_inner_coord, xcoord, ycoord)

    outer_systems = [System{Outer}(u_outer_coords[i], xcoord, ycoord)
                     for i in 1:N_outer_sys]

    [inner_system; outer_systems]
end

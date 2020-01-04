
abstract type CoordType end

abstract type Cartesian    <: CoordType end
abstract type GaussLobatto <: CoordType end

abstract type AbstractCoord{T,N,C} end

coord_axis(A::AbstractCoord{T,N,C}) where {T,N,C} = N
coord_type(A::AbstractCoord{T,N,C}) where {T,N,C} = C


struct Coord{T<:Real,N,C,T2} <: AbstractCoord{T,N,C}
    name  :: String
    min   :: T
    max   :: T
    nodes :: Int
    delta :: T2
end

struct CartCoord{T,N} end

function CartCoord{N}(name::String, xmin::T, xmax::T,
                      nodes::Integer) where {T<:Real,N}
    delta = (xmax - xmin) / (nodes - 1)
    Coord{T,N,Cartesian,typeof(delta)}(name, xmin, xmax, nodes, delta)
end
function CartCoord{N}(name::String, xmin::T, xmax::T, nodes::Integer;
                      endpoint::Bool=true) where {T<:Real,N}
    min_ = xmin
    if (!endpoint)
        h    = (xmax - xmin) / nodes
        max_ = xmax - h
    else
        h    = (xmax - xmin) / (nodes - 1)
        max_ = xmax
    end
    Coord{T,N,Cartesian,typeof(h)}(name, min_, max_, nodes, h)
end
CartCoord(args...) = CartCoord{1}(args...)


struct SpectralCoord{T,N} end

function SpectralCoord{N}(name::String, xmin::T, xmax::T,
                          nodes::Integer) where {T<:Real,N}
    Coord{T,N,GaussLobatto,Missing}(name, xmin, xmax, nodes, missing)
end

SpectralCoord(args...) = SpectralCoord{1}(args...)


@inline function xx(coord::AbstractCoord{T,N,Cartesian},
                    i::Int) where {T<:Real,N}
    coord.min + (i - 1) * coord.delta
end

@inline function xx(coord::AbstractCoord{T,N,GaussLobatto},
                    j::Int) where {T<:Real,N}
    xmin  = coord.min
    xmax  = coord.max
    M     = coord.nodes - 1.0
    xj    = -cos( (j - 1.0) * pi / M)
    0.5 * (xmax + xmin + (xmax - xmin) * xj)
end

@inline xx(coord::AbstractCoord) = [xx(coord, i) for i in 1:coord.nodes]


struct Grid{A}
    ndim    :: Int
    coords  :: A

    function Grid{A}(ndim, coords) where {A}
        ndim = length(coords)
        for a in 1:ndim
            @assert(coord_axis(coords[a]) == a, "wrong order in grid array")
        end
        new(ndim, coords)
    end
end

function Grid(coords::Vector)
    ndim = length(coords)
    Grid{typeof(coords)}(ndim, coords)
end

function Grid(coord::AbstractCoord)
    ndim   = 1
    coords = [coord]
    Grid{typeof(coords)}(ndim, coords)
end


@inline function xx(grid::Grid)
    [xx(grid.coords[a]) for a in 1:grid.ndim]
end

@inline function name(grid::Grid)
    [grid.coords[a].name for a in 1:grid.ndim]
end

@inline function min(grid::Grid)
    [grid.coords[a].min for a in 1:grid.ndim]
end

@inline function max(grid::Grid)
    [grid.coords[a].max for a in 1:grid.ndim]
end

@inline function nodes(grid::Grid)
    [grid.coords[a].nodes for a in 1:grid.ndim]
end

@inline function delta(grid::Grid)
    [grid.coords[a].delta for a in 1:grid.ndim]
end

@inline function coord_type(grid::Grid)
    [coord_type(grid.coords[a]) for a in 1:grid.ndim]
end


xmin        = -5.0
xmax        =  5.0
xnodes      =  128
ymin        = -5.0
ymax        =  5.0
ynodes      =  128
umin        =  0.0
umax        =  1.0
unodes      =  64


ucoord  = SpectralCoord{1}("u", umin, umax, unodes)
xcoord  = CartCoord{2}("x", xmin, xmax, xnodes, endpoint=false)
ycoord  = CartCoord{3}("y", ymin, ymax, ynodes, endpoint=false)


grid = Grid([ucoord, xcoord, ycoord])

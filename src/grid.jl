
abstract type CoordType end

abstract type Cartesian    <: CoordType end
abstract type GaussLobatto <: CoordType end

# TODO: consider making type of product of GaussLobatto grids? just need to
# redefine derivative ops, it may be better than having the array of ugrids...
# would have to do something about the inversion of the operator, though... it's
# better to invert one 64 x 64 matrix than 4 matrices 16 x 16. so would have to
# come up with a smart way of doing that as well...

abstract type AbstractCoord{T,N,C} end

coord_elem(A::AbstractCoord{T,N,C}) where {T,N,C} = T
coord_axis(A::AbstractCoord{T,N,C}) where {T,N,C} = N
coord_type(A::AbstractCoord{T,N,C}) where {T,N,C} = C

function coord_type(name::String)
    if name == "Cartesian"
        return Cartesian
    elseif name == "GaussLobatto"
        return GaussLobatto
    else
        error("Unknown coord type")
    end
end


struct Coord{T<:Real,N,C} <: AbstractCoord{T,N,C}
    name  :: String
    min   :: T
    max   :: T
    nodes :: Int
end

struct CartCoord{T,N} end

function CartCoord{N}(name::String, xmin::T, xmax::T, nodes::Integer;
                      endpoint::Bool=true) where {T<:Real,N}
    min_ = xmin
    if (!endpoint)
        h    = (xmax - xmin) / nodes
        max_ = xmax - h
    else
        max_ = xmax
    end
    Coord{T,N,Cartesian}(name, min_, max_, nodes)
end
CartCoord(args...) = CartCoord{1}(args...)


struct SpectralCoord{T,N} end

function SpectralCoord{N}(name::String, xmin::T, xmax::T,
                          nodes::Integer) where {T<:Real,N}
    Coord{T,N,GaussLobatto}(name, xmin, xmax, nodes)
end

SpectralCoord(args...) = SpectralCoord{1}(args...)


@inline function delta(coord::AbstractCoord{T,N,Cartesian}) where {T<:Real,N}
    (coord.max - coord.min) / (coord.nodes - 1)
end

# a Gauss-Lobatto grid has non-uniform spacing. not sure if the best is to
# return "missing", not define the method, or just have it return NaN, as now.
@inline delta(coord::AbstractCoord{T,N,GaussLobatto}) where {T<:Real,N} = NaN


@inline function Base.getindex(coord::AbstractCoord{T,N,Cartesian}, i::Int) where {T<:Real,N}
    h = delta(coord)
    coord.min + (i - 1) * h
end

@inline function Base.getindex(coord::AbstractCoord{T,N,GaussLobatto}, j::Int) where {T<:Real,N}
    xmin  = coord.min
    xmax  = coord.max
    M     = coord.nodes - 1.0
    xj    = -cos( (j - 1.0) * pi / M)
    0.5 * (xmax + xmin + (xmax - xmin) * xj)
end

@inline Base.firstindex(coord::AbstractCoord) = 1
@inline Base.lastindex(coord::AbstractCoord)  = coord.nodes
@inline Base.length(coord::AbstractCoord) = coord.nodes

@inline Base.getindex(coord::AbstractCoord, ::Colon) = [coord[i] for i in 1:coord.nodes]


struct Grid{T<:Real,A}
    ndim    :: Int
    coords  :: A

    function Grid{T,A}(ndim, coords) where {T,A}
        ndim = length(coords)
        for a in 1:ndim
            @assert(coord_axis(coords[a]) == a, "wrong order in grid array")
        end
        new(ndim, coords)
    end
end

function Grid(coords::Tuple)
    ndim = length(coords)
    T    = coord_elem(coords[1])
    Grid{T,typeof(coords)}(ndim, coords)
end

function Grid(coords::Vararg{AbstractCoord,N}) where {N}
    ndim = length(coords)
    T    = coord_elem(coords[1])
    Grid{T,typeof(coords)}(ndim, coords)
end

function Grid(coord::AbstractCoord)
    ndim   = 1
    coords = (coord)
    T      = coord_elem(coord)
    Grid{T,typeof(coords)}(ndim, coords)
end

@inline function Base.getindex(grid::Grid, idx::Vararg{Int,N}) where {N}
    [grid.coords[a][idx[a]] for a in 1:N]
end

@inline function Base.getindex(grid::Grid, ::Colon)
    [grid.coords[a][:] for a in 1:grid.ndim]
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
    [delta(grid.coords[a]) for a in 1:grid.ndim]
end

@inline function coord_type(grid::Grid)
    [coord_type(grid.coords[a]) for a in 1:grid.ndim]
end

@inline Base.size(grid::Grid) = Tuple(nodes(grid))

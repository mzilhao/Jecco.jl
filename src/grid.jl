
abstract type AbstractCoord{N} end

coord_axis(A::AbstractCoord{N}) where{N} = N

struct CartesianCoord{T<:Real,N} <: AbstractCoord{N}
    name  :: String
    min   :: T
    max   :: T
    nodes :: Int
end

struct GaussLobattoCoord{T<:Real,N} <: AbstractCoord{N}
    name  :: String
    min   :: T
    max   :: T
    nodes :: Int
end

coord_type(coord::CartesianCoord)    = "Cartesian"
coord_type(coord::GaussLobattoCoord) = "GaussLobatto"

struct Cartesian{N} end
function Cartesian{N}(name::String, xmin::T, xmax::T, nodes::Integer;
                      endpoint::Bool=true) where {T<:Real,N}
    min_ = xmin
    if (!endpoint)
        h    = (xmax - xmin) / nodes
        max_ = xmax - h
    else
        max_ = xmax
    end
    CartesianCoord{T,N}(name, min_, max_, nodes)
end
Cartesian(args...) = Cartesian{1}(args...)

struct GaussLobatto{N} end
function GaussLobatto{N}(name::String, xmin::T, xmax::T,
                          nodes::Integer) where {T<:Real,N}
    GaussLobattoCoord{T,N}(name, xmin, xmax, nodes)
end
GaussLobatto(args...) = GaussLobatto{1}(args...)


@inline function delta(coord::CartesianCoord) where {T<:Real,N}
    (coord.max - coord.min) / (coord.nodes - 1)
end

# a Gauss-Lobatto grid has non-uniform spacing. not sure if the best is to
# return "missing", not define the method, or just have it return NaN, as now.
@inline delta(coord::GaussLobattoCoord) = NaN

@inline function Base.getindex(coord::CartesianCoord, i::Int)
    h = delta(coord)
    coord.min + (i - 1) * h
end

@inline function Base.getindex(coord::GaussLobattoCoord, j::Int)
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


struct Grid{A}
    ndim    :: Int
    coords  :: A

    function Grid{A}(coords) where {A}
        ndim = length(coords)
        for a in 1:ndim
            @assert(coord_axis(coords[a]) == a, "wrong order in grid array")
        end
        new(ndim, coords)
    end
end
Grid(coords) = Grid{typeof(coords)}(coords)

function Grid(coords::Tuple)
    ndim = length(coords)
    Grid{typeof(coords)}(ndim, coords)
end

function Grid(coords::Vararg{AbstractCoord,N}) where {N}
    ndim = length(coords)
    Grid{typeof(coords)}(ndim, coords)
end

function Grid(coord::AbstractCoord)
    ndim   = 1
    coords = (coord)
    Grid{typeof(coords)}(ndim, coords)
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

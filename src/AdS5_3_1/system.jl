
"""
3D grid with the following configuration. x and y are periodic coordinates,
uniformly spaced. The u coordinate is split into an "inner" and "outer" portion;
we therefore will have an inner and an outer grid. The outer grid is further
decomposed into several u-domains, and in each of these domains the u-direction
uses Gauss-Lobatto points.

Since the boundary conditions for the nested system are always specified at u=0,
the inner grid necessarily starts at u=0 and finishes at u=u_outer_min.
"""
@with_kw struct Grid3D{T<:Real}
    x_min            :: T
    x_max            :: T
    x_nodes          :: Int
    y_min            :: T
    y_max            :: T
    y_nodes          :: Int
    u_outer_min      :: T
    u_outer_max      :: T
    u_outer_domains  :: Int     = 1
    u_outer_nodes    :: Int # number of points per domain
    u_inner_nodes    :: Int
end

# TODO: add a new method "Grid" acting on the above type


struct System{GT,Cu,Cx,Cy,Du,Dx,Dy}
    gridtype :: GT
    ucoord   :: Cu
    xcoord   :: Cx
    ycoord   :: Cy
    Du       :: Du
    Duu      :: Du
    Dx       :: Dx
    Dxx      :: Dx
    Dy       :: Dy
    Dyy      :: Dy
end

function System(gridtype::GT, ucoord::GaussLobattoCoord,
                xcoord::CartesianCoord, ycoord::CartesianCoord) where {GT<:GridType}
    ord = 4

    Du  = ChebDeriv{1}(1, ucoord.min, ucoord.max, ucoord.nodes)
    Duu = ChebDeriv{1}(2, ucoord.min, ucoord.max, ucoord.nodes)

    Dx  = CenteredDiff{2}(1, ord, Jecco.delta(xcoord), xcoord.nodes)
    Dxx = CenteredDiff{2}(2, ord, Jecco.delta(xcoord), xcoord.nodes)

    Dy  = CenteredDiff{3}(1, ord, Jecco.delta(ycoord), ycoord.nodes)
    Dyy = CenteredDiff{3}(2, ord, Jecco.delta(ycoord), ycoord.nodes)

    System{GT,typeof(ucoord), typeof(xcoord), typeof(ycoord), typeof(Du), typeof(Dx),
           typeof(Dy)}(gridtype, ucoord, xcoord, ycoord, Du, Duu, Dx, Dxx, Dy, Dyy)
end

Base.size(sys::System) = (sys.ucoord.nodes, sys.xcoord.nodes, sys.ycoord.nodes)


struct SystemPartition{N,S<:Tuple}
    _x :: S
end

function SystemPartition(x::AbstractVector{S}) where {S<:System}
    _x = Tuple(x)
    N = length(_x)
    SystemPartition{N,typeof(_x)}(_x)
end

@inline Base.iterate(ff::SystemPartition)         = iterate(ff._x)
@inline Base.iterate(ff::SystemPartition, i::Int) = iterate(ff._x, i)

@inline Base.length(ff::SystemPartition{N}) where{N} = N
@inline Base.firstindex(ff::SystemPartition) = 1
@inline Base.lastindex(ff::SystemPartition)  = length(ff)

@inline Base.getindex(ff::SystemPartition, i::Int) = ff._x[i]

"""
    SystemPartition(p::Grid3D)

Create a `SystemPartition`, which is a `Tuple` of `System` where each element
corresponds to a different u-domain
"""
function SystemPartition(p::Grid3D)
    u_inner_coord = GaussLobatto{1}("u", 0.0, p.u_outer_min, p.u_inner_nodes)

    N_outer_sys = p.u_outer_domains
    delta_udom  = (p.u_outer_max - p.u_outer_min) / N_outer_sys

    u_outer_coords =
        [GaussLobatto{1}("u", p.u_outer_min + (i-1)*delta_udom,
                          p.u_outer_min + i*delta_udom, p.u_outer_nodes)
         for i in 1:N_outer_sys]

    xcoord  = Cartesian{2}("x", p.x_min, p.x_max, p.x_nodes, endpoint=false)
    ycoord  = Cartesian{3}("y", p.y_min, p.y_max, p.y_nodes, endpoint=false)

    inner_system = System(Inner(), u_inner_coord, xcoord, ycoord)

    outer_systems = [System(Outer(), u_outer_coords[i], xcoord, ycoord)
                     for i in 1:N_outer_sys]

    SystemPartition([inner_system; outer_systems])
end


"""
    Boundary(p::Grid3D)

Create a `Boundary` struct with arrays of `size = (1,p.x_nodes,p.y_nodes)`
"""
function Boundary(p::Grid3D{T}) where {T}
    Nx = p.x_nodes
    Ny = p.y_nodes
    Boundary{T}(undef, Nx, Ny)
end

"""
    Gauge(p::Grid3D)

Create a `Gauge` struct with arrays of `size = (1,p.x_nodes,p.y_nodes)`
"""
function Gauge(p::Grid3D{T}) where {T}
    Nx = p.x_nodes
    Ny = p.y_nodes
    Gauge{T}(undef, Nx, Ny)
end

"""
    BulkEvols(p::Grid3D)

Create an `NTuple` of `BulkEvol` (with `length = 1 + p.u_outer_domains`) of elements
`BulkEvol`. The first `BulkEvol` has arrays of `size = (p.u_inner_nodes,
p.x_nodes, p.y_nodes)`, and the remaining ones have `size = (p.u_outer_nodes,
p.x_nodes, p.y_nodes)`
"""
function BulkEvols(p::Grid3D{T}) where {T}
    Nx = p.x_nodes
    Ny = p.y_nodes
    Nu_in  = p.u_inner_nodes
    Nu_out = p.u_outer_nodes
    N_outer_sys = p.u_outer_domains

    bulk_in  = [BulkEvol{T}(undef, Nu_in, Nx, Ny)]
    bulk_out = [BulkEvol{T}(undef, Nu_out, Nx, Ny) for i in 1:N_outer_sys]
    [bulk_in; bulk_out]
end

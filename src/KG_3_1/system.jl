
"""
3D grid with the following configuration. x and y are periodic coordinates,
uniformly spaced. The u coordinate is decomposed into several u-domains, and in
each of these domains the u-direction uses Gauss-Lobatto points.
"""
Base.@kwdef struct SpecCartGrid3D{T<:Real}
    x_min            :: T
    x_max            :: T
    x_nodes          :: Int
    y_min            :: T
    y_max            :: T
    y_nodes          :: Int
    u_min            :: T       = 0.0
    u_max            :: T
    u_domains        :: Int     = 1
    u_nodes          :: Int # number of points per domain
    fd_order         :: Int     = 4
end

function Jecco.Atlas(grid::SpecCartGrid3D{T}) where {T}
    Nsys = grid.u_domains
    delta_udom = (grid.u_max - grid.u_min) / Nsys

    ucoords = [GaussLobatto{1}("u", grid.u_min + (i-1)*delta_udom, grid.u_min +
                               i*delta_udom, grid.u_nodes) for i in 1:Nsys]

    xcoord  = Cartesian{2}("x", grid.x_min, grid.x_max, grid.x_nodes, endpoint=false)
    ycoord  = Cartesian{3}("y", grid.y_min, grid.y_max, grid.y_nodes, endpoint=false)

    charts  = [Chart(ucoords[i], xcoord, ycoord) for i in 1:Nsys]

    Atlas(charts)
end

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

function System(ucoord::GaussLobattoCoord, xcoord::CartesianCoord,
                ycoord::CartesianCoord, ord::Int)
    Du  = ChebDeriv{1}(1, ucoord.min, ucoord.max, ucoord.nodes)
    Duu = ChebDeriv{1}(2, ucoord.min, ucoord.max, ucoord.nodes)

    Dx  = CenteredDiff{2}(1, ord, Jecco.delta(xcoord), xcoord.nodes)
    Dxx = CenteredDiff{2}(2, ord, Jecco.delta(xcoord), xcoord.nodes)

    Dy  = CenteredDiff{3}(1, ord, Jecco.delta(ycoord), ycoord.nodes)
    Dyy = CenteredDiff{3}(2, ord, Jecco.delta(ycoord), ycoord.nodes)

    System{typeof(ucoord), typeof(xcoord), typeof(ycoord), typeof(Du), typeof(Dx),
           typeof(Dy)}(ucoord, xcoord, ycoord, Du, Duu, Dx, Dxx, Dy, Dyy)
end

Base.size(sys::System) = (sys.ucoord.nodes, sys.xcoord.nodes, sys.ycoord.nodes)


struct SystemPartition{N,A} <: AbstractPartition{N,A}
    x :: NTuple{N,A}
end

function SystemPartition(x_::AbstractVector{S}) where {S<:System}
    x = Tuple(x_)
    N = length(x)
    SystemPartition{N,eltype(x)}(x)
end

"""
    SystemPartition(grid::SpecCartGrid3D)

Create a `SystemPartition`, which is a `Tuple` of `System` where each element
corresponds to a different u-domain
"""
function SystemPartition(grid::SpecCartGrid3D{T}) where {T}
    Nsys = grid.u_domains
    delta_udom = (grid.u_max - grid.u_min) / Nsys

    ucoords = [GaussLobatto{1}("u", grid.u_min + (i-1)*delta_udom, grid.u_min +
                               i*delta_udom, grid.u_nodes) for i in 1:Nsys]

    xcoord  = Cartesian{2}("x", grid.x_min, grid.x_max, grid.x_nodes, endpoint=false)
    ycoord  = Cartesian{3}("y", grid.y_min, grid.y_max, grid.y_nodes, endpoint=false)

    systems = [System(ucoords[i], xcoord, ycoord, grid.fd_order) for i in 1:Nsys]

    SystemPartition(systems)
end

"""
    Boundary(grid::SpecCartGrid3D)

Create a `Boundary` struct with arrays of `size = (1,grid.x_nodes,grid.y_nodes)`
"""
function Boundary(grid::SpecCartGrid3D{T}) where {T}
    Nx = grid.x_nodes
    Ny = grid.y_nodes
    Boundary{T}(undef, Nx, Ny)
end

"""
    BulkEvolvedPartition(grid::SpecCartGrid3D)

Returns `BulkPartition`, of `length = grid.u_domains`, with
elements of type `BulkEvolved`.
"""
function BulkEvolvedPartition(grid::SpecCartGrid3D{T}) where {T}
    Nx   = grid.x_nodes
    Ny   = grid.y_nodes
    Nu   = grid.u_nodes
    Nsys = grid.u_domains

    bulk  = [BulkEvolved{T}(undef, Nu, Nx, Ny) for i in 1:Nsys]

    BulkPartition(bulk...)
end

"""
    BulkConstrainedPartition(grid::SpecCartGrid3D)

Returns `BulkPartition`, of `length = grid.u_domains`, with
elements of type `BulkConstrained`.
"""
function BulkConstrainedPartition(grid::SpecCartGrid3D{T}) where {T}
    Nx   = grid.x_nodes
    Ny   = grid.y_nodes
    Nu   = grid.u_nodes
    Nsys = grid.u_domains

    bulk  = [BulkConstrained{T}(undef, Nu, Nx, Ny) for i in 1:Nsys]

    BulkPartition(bulk...)
end

function BulkDerivPartition(grid::SpecCartGrid3D{T}) where {T}
    Nx   = grid.x_nodes
    Ny   = grid.y_nodes
    Nu   = grid.u_nodes
    Nsys = grid.u_domains

    bulk = [BulkDeriv{T}(undef, Nu, Nx, Ny) for i in 1:Nsys]

    BulkPartition(bulk...)
end

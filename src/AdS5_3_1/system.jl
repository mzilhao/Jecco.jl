
abstract type GridType end
struct Inner <: GridType end
struct Outer <: GridType end

"""
3D grid with the following configuration. x and y are periodic coordinates,
uniformly spaced. The u coordinate is split into an "inner" and "outer" portion;
we therefore will have an inner and an outer grid. The outer grid is further
decomposed into several u-domains, and in each of these domains the u-direction
uses Gauss-Lobatto points.

Since the boundary conditions for the nested system are always specified at u=0,
the inner grid necessarily starts at u=0 and finishes at u=u_outer_min.
"""
Base.@kwdef struct SpecCartGrid3D{T<:Real}
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
    fd_order         :: Int     = 4
    filter_γ         :: T       = 8.0
    sigma_diss       :: T       = 0.2
end

function Jecco.Atlas(grid::SpecCartGrid3D{T}) where {T}
    u_inner_coord = GaussLobatto{1}("u", zero(T), grid.u_outer_min, grid.u_inner_nodes)

    N_outer_sys = grid.u_outer_domains
    delta_udom  = (grid.u_outer_max - grid.u_outer_min) / N_outer_sys

    u_outer_coords =
        [GaussLobatto{1}("u", grid.u_outer_min + (i-1)*delta_udom,
                          grid.u_outer_min + i*delta_udom, grid.u_outer_nodes)
         for i in 1:N_outer_sys]

    xcoord  = Cartesian{2}("x", grid.x_min, grid.x_max, grid.x_nodes, endpoint=false)
    ycoord  = Cartesian{3}("y", grid.y_min, grid.y_max, grid.y_nodes, endpoint=false)

    inner_chart  = Chart(u_inner_coord, xcoord, ycoord)
    outer_charts = [Chart(u_outer_coords[i], xcoord, ycoord)
                     for i in 1:N_outer_sys]

    Atlas([inner_chart; outer_charts])
end

struct Filters{F1,F2,F3,F4,F5}
    exp_filter    :: F1
    ko_filter_x   :: F2
    ko_filter_y   :: F3
    ko_filter2D_x :: F4
    ko_filter2D_y :: F5
end
function Filters(filter_gamma::T, KO_order::Int, sigma_diss::T,
                 Nu::Int, Nx::Int, Ny::Int) where {T}
    exp_filter    = Exp_Filter{1}(filter_gamma, Nu, Nx, Ny)

    ko_filter_x   = KO_Filter{2}(KO_order, sigma_diss, Nu, Nx, Ny)
    ko_filter_y   = KO_Filter{3}(KO_order, sigma_diss, Nu, Nx, Ny)

    ko_filter2D_x = KO_Filter{2}(KO_order, sigma_diss, 1, Nx, Ny)
    ko_filter2D_y = KO_Filter{3}(KO_order, sigma_diss, 1, Nx, Ny)

    Filters{typeof(exp_filter), typeof(ko_filter_x), typeof(ko_filter_y), typeof(ko_filter2D_x),
            typeof(ko_filter2D_y)}(exp_filter, ko_filter_x, ko_filter_y,
                                   ko_filter2D_x, ko_filter2D_y)
end

struct System{GT,Cu,Cx,Cy,Du,Dx,Dy,TI,TF}
    gridtype    :: GT
    ucoord      :: Cu
    xcoord      :: Cx
    ycoord      :: Cy
    Du          :: Du
    Duu         :: Du
    Dx          :: Dx
    Dxx         :: Dx
    Dy          :: Dy
    Dyy         :: Dy
    uinterp     :: TI
    filters     :: TF
end

function System(gridtype::GT, ucoord::GaussLobattoCoord,
                xcoord::CartesianCoord, ycoord::CartesianCoord, ord::Int,
                filter_gamma::T, sigma_diss::T) where {GT<:GridType,T<:Real}
    Du  = ChebDeriv{1}(1, ucoord.min, ucoord.max, ucoord.nodes)
    Duu = ChebDeriv{1}(2, ucoord.min, ucoord.max, ucoord.nodes)

    Dx  = CenteredDiff{2}(1, ord, Jecco.delta(xcoord), xcoord.nodes)
    Dxx = CenteredDiff{2}(2, ord, Jecco.delta(xcoord), xcoord.nodes)

    Dy  = CenteredDiff{3}(1, ord, Jecco.delta(ycoord), ycoord.nodes)
    Dyy = CenteredDiff{3}(2, ord, Jecco.delta(ycoord), ycoord.nodes)

    uinterp = ChebInterpolator(ucoord.min, ucoord.max, ucoord.nodes)

    KO_order = ord + 1
    filters  = Filters(filter_gamma, KO_order, sigma_diss, ucoord.nodes,
                       xcoord.nodes, ycoord.nodes)

    System{GT, typeof(ucoord), typeof(xcoord), typeof(ycoord), typeof(Du),
           typeof(Dx), typeof(Dy), typeof(uinterp),
           typeof(filters)}(gridtype, ucoord, xcoord, ycoord,
                            Du, Duu, Dx, Dxx, Dy, Dyy, uinterp, filters)
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
    SystemPartition(grid::SpecCartGrid3D)

Create a `SystemPartition`, which is a `Tuple` of `System` where each element
corresponds to a different u-domain
"""
function SystemPartition(grid::SpecCartGrid3D{T}) where {T}
    u_inner_coord = GaussLobatto{1}("u", zero(T), grid.u_outer_min, grid.u_inner_nodes)

    N_outer_sys = grid.u_outer_domains
    delta_udom  = (grid.u_outer_max - grid.u_outer_min) / N_outer_sys

    u_outer_coords =
        [GaussLobatto{1}("u", grid.u_outer_min + (i-1)*delta_udom,
                          grid.u_outer_min + i*delta_udom, grid.u_outer_nodes)
         for i in 1:N_outer_sys]

    xcoord  = Cartesian{2}("x", grid.x_min, grid.x_max, grid.x_nodes, endpoint=false)
    ycoord  = Cartesian{3}("y", grid.y_min, grid.y_max, grid.y_nodes, endpoint=false)

    inner_system = System(Inner(), u_inner_coord, xcoord, ycoord, grid.fd_order,
                          grid.filter_γ, grid.sigma_diss)

    outer_systems = [System(Outer(), u_outer_coords[i], xcoord, ycoord, grid.fd_order,
                            grid.filter_γ, grid.sigma_diss) for i in 1:N_outer_sys]

    SystemPartition([inner_system; outer_systems])
end


"""
    Boundary(grid::SpecCartGrid3D)

Create a `Boundary` struct with arrays of `size = (1,p.x_nodes,p.y_nodes)`
"""
function Boundary(grid::SpecCartGrid3D{T}) where {T}
    Nx = grid.x_nodes
    Ny = grid.y_nodes
    Boundary{T}(undef, Nx, Ny)
end

"""
    Gauge(grid::SpecCartGrid3D)

Create a `Gauge` struct with arrays of `size = (1,grid.x_nodes,grid.y_nodes)`
"""
function Gauge(grid::SpecCartGrid3D{T}) where {T}
    Nx = grid.x_nodes
    Ny = grid.y_nodes
    Gauge{T}(undef, Nx, Ny)
end

"""
    BulkEvolvedPartition(grid::SpecCartGrid3D)

Returns `BulkPartition`, of `length = 1 + grid.u_outer_domains`, with
elements of type `BulkEvolved`. The first `BulkEvolved` has arrays of `size =
(grid.u_inner_nodes, grid.x_nodes, grid.y_nodes)`, and the remaining ones have
`size = (grid.u_outer_nodes, grid.x_nodes, grid.y_nodes)`
"""
function BulkEvolvedPartition(grid::SpecCartGrid3D{T}) where {T}
    Nx = grid.x_nodes
    Ny = grid.y_nodes
    Nu_in  = grid.u_inner_nodes
    Nu_out = grid.u_outer_nodes
    N_outer_sys = grid.u_outer_domains

    bulk_in   = BulkEvolved{T}(undef, Nu_in, Nx, Ny)
    bulk_out  = [BulkEvolved{T}(undef, Nu_out, Nx, Ny) for i in 1:N_outer_sys]

    BulkPartition((bulk_in, bulk_out...))
end

"""
    BulkConstrainedPartition(grid::SpecCartGrid3D)

Returns `BulkPartition`, of `length = 1 + grid.u_outer_domains`, with
elements of type `BulkConstrained`. The first `BulkConstrained` has arrays of `size =
(grid.u_inner_nodes, grid.x_nodes, grid.y_nodes)`, and the remaining ones have
`size = (grid.u_outer_nodes, grid.x_nodes, grid.y_nodes)`
"""
function BulkConstrainedPartition(grid::SpecCartGrid3D{T}) where {T}
    Nx = grid.x_nodes
    Ny = grid.y_nodes
    Nu_in  = grid.u_inner_nodes
    Nu_out = grid.u_outer_nodes
    N_outer_sys = grid.u_outer_domains

    bulk_in   = BulkConstrained{T}(undef, Nu_in, Nx, Ny)
    bulk_out  = [BulkConstrained{T}(undef, Nu_out, Nx, Ny) for i in 1:N_outer_sys]

    BulkPartition((bulk_in, bulk_out...))
end
"""
    BulkConstrainedPartition(systems::SystemPartition)

Returns `BulkPartition`, of `length = length(systems)`
"""
function BulkConstrainedPartition(systems::SystemPartition)
    T     = Jecco.coord_eltype(systems[1].ucoord)
    bulks = [BulkConstrained{T}(undef, size(sys)...) for sys in systems]

    BulkPartition(bulks...)
end

function BulkDerivPartition(grid::SpecCartGrid3D{T}) where {T}
    Nx = grid.x_nodes
    Ny = grid.y_nodes
    Nu_in  = grid.u_inner_nodes
    Nu_out = grid.u_outer_nodes
    N_outer_sys = grid.u_outer_domains

    bulk_in   = BulkDeriv{T}(undef, Nu_in, Nx, Ny)
    bulk_out  = [BulkDeriv{T}(undef, Nu_out, Nx, Ny) for i in 1:N_outer_sys]

    BulkPartition((bulk_in, bulk_out...))
end
function BulkDerivPartition(systems::SystemPartition)
    T     = Jecco.coord_eltype(systems[1].ucoord)
    bulks = [BulkDeriv{T}(undef, size(sys)...) for sys in systems]

    BulkPartition(bulks...)
end

function HorizonCache(sys::System, ord::Int)
    _, Nx, Ny = size(sys)
    hx = Jecco.delta(sys.xcoord)
    hy = Jecco.delta(sys.ycoord)
    T1 = Jecco.coord_eltype(sys.ucoord)
    M  = Nx * Ny

    # there is no method for lu when acting on SparseMatrixCSC{BigFloat,Int64}, so
    # let's use instead type Float64 for such cases
    if T1 == BigFloat
        T2 = Float64
    else
        T2 = T1
    end

    bulkhorizon = BulkHorizon{T1}(Nx, Ny)

    axx    = Vector{T2}(undef, M)
    ayy    = Vector{T2}(undef, M)
    axy    = Vector{T2}(undef, M)
    bx     = Vector{T2}(undef, M)
    by     = Vector{T2}(undef, M)
    cc     = Vector{T2}(undef, M)
    b_vec  = Vector{T2}(undef, M)

    Dx_    = CenteredDiff{1}(1, ord, hx, Nx)
    Dxx_   = CenteredDiff{1}(2, ord, hx, Nx)
    Dy_    = CenteredDiff{2}(1, ord, hy, Ny)
    Dyy_   = CenteredDiff{2}(2, ord, hy, Ny)

    #=
    use the Kronecker product (kron) to build 2-dimensional derivation matrices
    from 1-dimensional ones. see for instance:

    https://en.wikipedia.org/wiki/Kronecker_product
    https://arxiv.org/pdf/1801.01483.pdf (section 5)
    =#
    Dx_2D  = T2.(kron(I(Ny), SparseMatrixCSC(Dx_)))
    Dxx_2D = T2.(kron(I(Ny), SparseMatrixCSC(Dxx_)))
    Dy_2D  = T2.(kron(SparseMatrixCSC(Dy_), I(Nx)))
    Dyy_2D = T2.(kron(SparseMatrixCSC(Dyy_), I(Nx)))
    Dxy_2D = T2.(Dx_2D * Dy_2D)

    _Dx_2D  = copy(Dx_2D)
    _Dxx_2D = copy(Dxx_2D)
    _Dy_2D  = copy(Dy_2D)
    _Dyy_2D = copy(Dyy_2D)
    _Dxy_2D = copy(Dxy_2D)

    HorizonCache{T1,T2,typeof(Dx_2D)}(bulkhorizon, axx, ayy, axy, bx, by, cc, b_vec,
                                  Dx_2D,  Dy_2D, Dxx_2D,  Dyy_2D,  Dxy_2D,
                                  _Dx_2D, _Dy_2D, _Dxx_2D, _Dyy_2D, _Dxy_2D)
end

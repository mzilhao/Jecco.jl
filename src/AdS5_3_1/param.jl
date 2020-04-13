

# TODO: add phi0 param here

@with_kw struct ParamBase
    which_potential :: String
end


#= Grid

We use a 3D grid with the following configuration. x and y are periodic
coordinates, uniformly spaced. The u coordinate is split into an "inner" and
"outer" portion; we therefore will have an inner and an outer grid. The outer
grid is further decomposed into several u-domains, and in each of these domains
the u-direction uses Gauss-Lobatto points.

Since the boundary conditions for the nested system are always specified at u=0,
the inner grid necessarily starts at u=0 and finishes at u=u_outer_min.

=#

@with_kw struct ParamGrid
    x_min            :: Float64
    x_max            :: Float64
    x_nodes          :: Int
    y_min            :: Float64
    y_max            :: Float64
    y_nodes          :: Int
    u_outer_min      :: Float64
    u_outer_max      :: Float64
    u_outer_domains  :: Int     = 1
    u_outer_nodes    :: Int # number of points per domain
    u_inner_nodes    :: Int
end

@with_kw struct ParamID
    ID_type     :: String
    A0x         :: Float64  = 0.0
    A0y         :: Float64  = 0.0
    Lx          :: Float64  = 1.0
    Ly          :: Float64  = 1.0
end


# TODO: add kappa param here

@with_kw struct ParamEvol
    dt          :: Float64
    tmax        :: Float64
    ODE_method  :: String   = "RK4"
end

@with_kw struct ParamIO
    out_every   :: Int
    folder      :: String  = "./data"
    prefix      :: String  = "phi_"
    overwrite   :: Bool    = false
end


Base.@kwdef struct ParamBase
    which_potential :: String
end

Base.@kwdef struct ParamGrid
    xmin        :: Float64
    xmax        :: Float64
    xnodes      :: Int
    ymin        :: Float64
    ymax        :: Float64
    ynodes      :: Int
    umin        :: Float64
    umax        :: Float64
    udomains    :: Int     = 1
    unodes      :: Int # number of points per domain
end

Base.@kwdef struct ParamID
    ID_type     :: String
    A0x         :: Float64  = 0.0
    A0y         :: Float64  = 0.0
    Lx          :: Float64  = 1.0
    Ly          :: Float64  = 1.0
end

Base.@kwdef struct ParamEvol
    dt          :: Float64
    tmax        :: Float64
    ODE_method  :: String   = "RK4"
end

Base.@kwdef struct ParamIO
    out_every   :: Int
    folder      :: String  = "./data"
    prefix      :: String  = "phi_"
    overwrite   :: Bool    = false
end

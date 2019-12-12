
@with_kw struct ParamBase
    which_potential :: String
end

@with_kw struct ParamGrid
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

@with_kw struct ParamID
    ID_type      :: String
    A0x          :: Float64
    A0y          :: Float64
    Lx           :: Float64
    Ly           :: Float64
end

@with_kw struct ParamEvol
    dt          :: Float64
    tmax        :: Float64
end

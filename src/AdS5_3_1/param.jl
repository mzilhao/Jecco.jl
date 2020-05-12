
# @with_kw struct ParamID
#     ID_type     :: String
#     A0x         :: Float64  = 0.0
#     A0y         :: Float64  = 0.0
#     Lx          :: Float64  = 1.0
#     Ly          :: Float64  = 1.0
# end


@with_kw struct ParamEvol
    dt          :: Float64
    tmax        :: Float64
    ODE_method  :: String   = "RK4"
    kappa       :: Float64  = 1.0
end

@with_kw struct ParamIO
    out_every   :: Int
    folder      :: String  = "./data"
    prefix      :: String  = "phi_"
    overwrite   :: Bool    = false
end

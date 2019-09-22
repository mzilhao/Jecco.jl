module KG_3_1

using Jecco
using Vivi
using Parameters

export Param
export System

@with_kw struct Param
    A0x         :: Float64
    A0y         :: Float64

    tmax        :: Float64
    out_every_t :: Float64

    xmin        :: Float64
    xmax        :: Float64
    xnodes      :: Int
    ymin        :: Float64
    ymax        :: Float64
    ynodes      :: Int
    umin        :: Float64
    umax        :: Float64
    unodes      :: Int

    dtfac       :: Float64    = 0.5

    folder      :: String  = "./data"
end


struct System{C,D} <: Vivi.System
    coords :: C
    derivs :: D
    # _dt    :: Float64
    # param  :: Param
end

function System(p::Param)
    ucoord  = Vivi.SpectralCoord("u", p.umin, p.umax, p.unodes)

    xcoord  = Vivi.CartCoord("x", p.xmin, p.xmax, p.xnodes, endpoint=false)
    ycoord  = Vivi.CartCoord("y", p.ymin, p.ymax, p.ynodes, endpoint=false)

    coords = Vivi.CoordSystem{Float64}("uxy", [ucoord, xcoord, ycoord])

    # FIXME
    ord    = 4
    BC     = :periodic

    # dx     = xcoord.delta
    # dy     = ycoord.delta
    # dt0    = p.dtfac * min(dx, dy)

    derivs = Vivi.Deriv(coords, (nothing, ord, ord), (nothing, BC, BC))

    # System{typeof(coords), typeof(derivs)}(coords, derivs, dt0, p)
    System{typeof(coords), typeof(derivs)}(coords, derivs)
end

include("initial_data.jl")

end

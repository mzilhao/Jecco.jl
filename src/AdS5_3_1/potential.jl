
abstract type PotentialType end
struct ZeroPotential   <: PotentialType end
struct SquarePotential <: PotentialType end

struct Potential{T<:PotentialType,fType,fpType}
    _type :: T
    f     :: fType
    fp    :: fpType
end

(ff::Potential)(phi)  = ff.f(phi)
∂(ff::Potential)      = ff.fp

Potential(::ZeroPotential) = Potential(ZeroPotential(), phi -> -3, phi -> 0)

Potential(::SquarePotential) = Potential(SquarePotential(),
                                         phi -> -3 + phi*phi / 2, phi -> phi)

function Potential(par_base::ParamBase)
    if par_base.which_potential == :zero
        Potential(ZeroPotential())
    elseif par_base.which_potential == :square
        Potential(SquarePotential())
    else
        error("Unknown potential.")
    end
end

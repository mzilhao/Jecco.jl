
abstract type AbstractPotential{fType,fpType} end

(ff::AbstractPotential)(phi)   = ff.f(phi)
∂(ff::AbstractPotential)       = ff.fp

struct ZeroPotential{fType,fpType} <: AbstractPotential{fType,fpType}
    f  :: fType
    fp :: fpType
end
ZeroPotential() = ZeroPotential(phi -> -3, phi -> 0)

struct SquarePotential{fType,fpType} <: AbstractPotential{fType,fpType}
    f  :: fType
    fp :: fpType
end
SquarePotential() = SquarePotential(phi -> -3 + phi*phi / 2, phi -> phi)

function Potential(par_base::ParamBase)
    if par_base.which_potential == :zero
        ZeroPotential()
    elseif par_base.which_potential == :square
        SquarePotential()
    else
        error("Unknown potential.")
    end
end

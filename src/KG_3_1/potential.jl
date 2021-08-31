
struct ConstPotential{T} <: Potential end
ConstPotential() = ConstPotential{Float64}()

UU(phi,  ::ConstPotential) = phi
UUp(phi, ::ConstPotential{T}) where {T} = one(T)

parameters(potential::ConstPotential) = ()


struct SquarePotential{T} <: Potential end
SquarePotential() = SquarePotential{Float64}()

UU(phi,  ::SquarePotential) = -1 + phi*phi / 2
UUp(phi, ::SquarePotential) = phi

parameters(potential::SquarePotential) = ()

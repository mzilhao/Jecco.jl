
struct ZeroPotential <: Potential end

UU(phi,  ::ZeroPotential) = 0
UUp(phi, ::ZeroPotential) = 0

parameters(potential::ZeroPotential) = ()

"""
* `oophiM2` : 1/ϕM^2
* `oophiQ`  : 1/ϕQ
"""
Base.@kwdef struct Phi8Potential{T} <: Potential
    oophiM2 :: T = 0.0
    oophiQ  :: T = 0.0
end

parameters(potential::Phi8Potential) = (oophiM2=potential.oophiM2, oophiQ=potential.oophiQ)

@doc raw"""
```math
U = -\frac{1}{3} + \frac{ϕ^2}{3 ϕM^2} + \frac{ϕ^2}{2 ϕM^4} - \frac{2 ϕ^2}{ϕQ}
 + \frac{4 ϕ^4}{3 ϕQ} + \frac{6 ϕ^4}{ϕM^2 ϕQ} - \frac{ϕ^4}{12 ϕM^4}
 + \frac{18 ϕ^6}{ϕQ^2} - \frac{2 ϕ^6}{3 ϕM^2 ϕQ}
 - \frac{4ϕ^8}{3 ϕQ^2}
```
"""
function UU(phi, potential::Phi8Potential)
    phi2 = phi  * phi
    phi4 = phi2 * phi2
    phi6 = phi2 * phi4
    phi8 = phi4 * phi4

    oophiM2 = potential.oophiM2
    oophiM4 = oophiM2 * oophiM2
    oophiQ  = potential.oophiQ
    oophiQ2 = oophiQ * oophiQ

    -1/3 + phi2 * (oophiM4 / 2 + oophiM2 / 3 - 2 * oophiQ) +
        phi4 * ( -oophiM4 / 12 + 4 * oophiQ / 3 + 6 * oophiM2 * oophiQ) +
        phi6 * (18 * oophiQ2 - 2 * oophiM2 * oophiQ / 3) -
        phi8 * 4 * oophiQ2 / 3
end

@doc raw"""
```math
Up = \frac{dU}{dϕ}
```
"""
function UUp(phi, potential::Phi8Potential)
    phi2 = phi  * phi
    phi3 = phi  * phi2
    phi5 = phi3 * phi2
    phi7 = phi5 * phi2

    oophiM2 = potential.oophiM2
    oophiM4 = oophiM2 * oophiM2
    oophiQ  = potential.oophiQ
    oophiQ2 = oophiQ * oophiQ

    phi * (oophiM4 + 2 * oophiM2 / 3 - 4 * oophiQ) +
        phi3 * ( -oophiM4 / 3 + 16 * oophiQ / 3 + 24 * oophiM2 * oophiQ) +
        phi5 * (108 * oophiQ2 - 4 * oophiM2 * oophiQ) -
        phi7 * 32 * oophiQ2 / 3
end


#New potential coming from a superpotential with extra term alpha*phi^8

Base.@kwdef struct PhiAlphaPotential{T} <: Potential
    oophiM2 :: T = 0.0
    oophiQ  :: T = 0.0
    alpha   :: T = 0.0
end

parameters(potential::PhiAlphaPotential) = (oophiM2=potential.oophiM2, oophiQ=potential.oophiQ, alpha=potential.alpha)

function UU(phi, potential::PhiAlphaPotential)
    phi2  = phi  * phi
    phi4  = phi2 * phi2
    phi6  = phi2 * phi4
    phi8  = phi4 * phi4
    phi10 = phi8 * phi2
    phi12 = phi6 * phi6

    oophiM2 = potential.oophiM2
    oophiM4 = oophiM2 * oophiM2
    oophiQ  = potential.oophiQ
    oophiQ2 = oophiQ * oophiQ
    alpha   = potential.alpha
    alpha2  = alpha * alpha

    return -1/3 + (1/2 * oophiM4 + 1/3 * oophiM2 - 2 * oophiQ) * phi2 + (-4 * alpha - 1/12 * oophiM4 +
    4/3 * oophiQ + 6 * oophiM2 * oophiQ) * phi4 + (alpha * (4/3 + 8 * oophiM2) + 18 * oophiQ2 -
    2/3 * oophiM2 * oophiQ) * phi6 + (-2/3 * alpha * oophiM2 - 4/3 * oophiQ2 + 48 * alpha * oophiQ) * phi8 +
    (32 * alpha2 - 8/3 * alpha * oophiQ) * phi10 - 4/3 * alpha2 * phi12
end

@doc raw"""
```math
Up = \frac{dU}{dϕ}
```
"""
function UUp(phi, potential::PhiAlphaPotential)
    phi2  = phi  * phi
    phi3  = phi  * phi2
    phi5  = phi3 * phi2
    phi7  = phi5 * phi2
    phi9  = phi7 * phi2
    phi11 = phi9 * phi2

    oophiM2 = potential.oophiM2
    oophiM4 = oophiM2 * oophiM2
    oophiQ  = potential.oophiQ
    oophiQ2 = oophiQ * oophiQ
    alpha   = potential.alpha
    alpha2  = alpha * alpha

    return (oophiM4 + 2/3 * oophiM2 - 4 * oophiQ) * phi + (-16 * alpha - 1/3 * oophiM4 + 16/3 * oophiQ
    + 24 * oophiM2 * oophiQ) * phi3 + (8 * alpha + 48 * alpha * oophiM2 + 108 * oophiQ2
    - 4 * oophiM2 * oophiQ) * phi5 + 16/3 * (-alpha * oophiM2 -2 * oophiQ2 + 72 * alpha * oophiQ) * phi7
    + (320 * alpha2 - 80/3 * alpha * oophiQ) * phi9 - 16 * alpha2 * phi11

end

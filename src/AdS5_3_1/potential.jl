
struct ZeroPotential <: Potential end

UU(phi,  ::ZeroPotential) = 0
UUp(phi, ::ZeroPotential) = 0

parameters(potential::ZeroPotential) = ()

"""
* `oophiM2` : 1/ϕM^2
* `oophiQ`  : 1/ϕQ

Note that, in the notation used in the paper XXXX.XXXX,
```math
λ₄ = 1/(4 ϕM^2), λ₆ = 1/ϕQ
```
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

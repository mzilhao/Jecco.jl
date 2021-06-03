
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


#New potential coming from a superpotential with extra term alpha*phi^8

Base.@kwdef struct PhiAlphaBetaPotential{T} <: Potential
    oophiM2 :: T   = 0.0
    oophiQ  :: T   = 0.0
    alpha   :: T   = 0.0
    beta    :: Int = 1
end

parameters(potential::PhiAlphaBetaPotential) = (oophiM2=potential.oophiM2, oophiQ=potential.oophiQ, alpha=potential.alpha)

function UU(phi, potential::PhiAlphaBetaPotential)
    phi2         = phi  * phi
    phi4         = phi2 * phi2
    phi6         = phi2 * phi4
    phi8         = phi4 * phi4

    oophiM2 = potential.oophiM2
    oophiM4 = oophiM2 * oophiM2
    oophiQ  = potential.oophiQ
    oophiQ2 = oophiQ * oophiQ
    alpha   = potential.alpha
    alpha2  = alpha * alpha
    beta    = potential.beta
    beta2   = beta * beta

    phi4beta     = phi^(beta-4)
    phi2beta     = phi4beta * phi2
    phiplus2beta = phi2beta * phi4
    phibeta      = phi4beta * phi4
    phi62beta    = phi^(2*beta-6)

    return -1/3+4*alpha*phi4beta-alpha*beta*phi4beta+4/3*alpha*phi2beta+
            1/2*alpha2*beta2*phi62beta-4/3*alpha2*phi2beta^2+1/2*phi2*oophiM4-
            1/12*oophiM4*phi4+1/3*phi2*oophiM2+alpha*beta*oophiM2*phi2beta-
            2/3*alpha*phibeta*oophiM2+18*oophiQ2*phi6-4/3*oophiQ2*phi8-2*oophiQ*phi2+
            4/3*oophiQ*phi4+6*alpha*beta*oophiQ*phibeta-8/3*alpha*oophiQ*phiplus2beta+
            6*oophiM2*oophiQ*phi4-2/3*oophiM2*oophiQ*phi6

end

@doc raw"""
```math
Up = \frac{dU}{dϕ}
```
"""
function UUp(phi, potential::PhiAlphaBetaPotential)
    phi2  = phi  * phi
    phi3  = phi  * phi2
    phi4  = phi  * phi3
    phi5  = phi3 * phi2
    phi7  = phi5 * phi2

    oophiM2 = potential.oophiM2
    oophiM4 = oophiM2 * oophiM2
    oophiQ  = potential.oophiQ
    oophiQ2 = oophiQ * oophiQ
    alpha   = potential.alpha
    alpha2  = alpha * alpha
    beta    = potential.beta
    beta2   = beta * beta
    beta3   = beta2 * beta

    phi5beta     = phi^(beta-5)
    phi1beta     = phi5beta * phi4
    phi3beta     = phi5beta * phi2
    phiplus1beta = phi1beta * phi2
    phi72beta    = phi^(2*beta-7)
    phi52beta    = phi72beta * phi2

    return -16*alpha*phi5beta+8*alpha*beta*phi5beta-alpha*beta2*phi5beta-
            8/3*alpha*phi3beta+4/3*alpha*beta*phi3beta-3*alpha2*beta2*phi72beta+
            alpha2*beta3*phi72beta+16/3*alpha2*phi52beta-8/3*alpha2*beta*phi52beta+
            oophiM4*phi-1/3*oophiM4*phi3+2/3*oophiM2*phi-2*alpha*beta*oophiM2*phi3beta+
            alpha*beta2*oophiM2*phi3beta-
            2/3*alpha*beta*oophiM2*phi1beta+108*oophiQ2*phi5-32/3*oophiQ2*phi7-
            4*oophiQ*phi+16/3*oophiQ*phi3+6*alpha*beta2*oophiQ*phi1beta-
            16/3*alpha*oophiQ*phiplus1beta-8/3*alpha*beta*oophiQ*phiplus1beta+
            24*oophiM2*oophiQ*phi3-4*oophiM2*oophiQ*phi5

end

Base.@kwdef struct PhiPoli{T} <: Potential
    alpha :: T   = 0.0
    beta  :: T   = 0.0
    gamma :: T   = 0.0
end

parameters(potential::PhiPoli) = (alpha=potential.alpha, beta=potential.beta, gamma=potential.gamma)


function UU(phi, potential::PhiPoli)
    phi2 = phi  * phi
    phi4 = phi2 * phi2
    phi6 = phi2 * phi4

    alpha = potential.alpha
    beta  = potential.beta
    gamma = potential.gamma

    -1/3 + alpha * phi2 + beta * phi4 + gamma * phi6
end

@doc raw"""
```math
Up = \frac{dU}{dϕ}
```
"""
function UUp(phi, potential::PhiPoli)
    phi2 = phi  * phi
    phi3 = phi  * phi2
    phi5 = phi2 * phi3

    alpha = potential.alpha
    beta  = potential.beta
    gamma = potential.gamma

    2 * alpha * phi + 4 * beta * phi3 + 6 * gamma * phi5
end


# Vf(phi)  = VV(phi)
# Vfp(phi) = ∂(VV)(phi)

# FIXME
#V(phi) = -3.0
#Vp(phi) = 0.0


# assuming
# (A d_uu + B d_u + C Id) f = -S

function S_eq_coeff!(ABCS::Vector, vars::AllVars{Inner})
    u    = vars.u
    phi0 = vars.phi0

    xi   = vars.xi

    B1   = vars.B1
    B1p  = vars.B1p

    B2   = vars.B2
    B2p  = vars.B2p

    G   = vars.G
    Gp  = vars.Gp

    phi  = vars.phi
    phip = vars.phip

    coshGu4sq = cosh(*(G, u ^ 4))


    ABCS[1] = *(6, u ^ 7)

    ABCS[2] = *(48, u ^ 6)

    ABCS[3] = *(u ^ 5, 72 + *(Gp ^ 2, u ^ 6) + *(3, B2p ^ 2, u ^ 6) + *(16, G ^ 2, u ^ 8) + *(48, B2 ^ 2, u ^ 8) + *(u ^ 6, (B1p + *(-4, B1, u)) ^ 2, coshGu4sq) + *(-24, B2, B2p, u ^ 7) + *(-8, G, Gp, u ^ 7) + *(4, phi0 ^ 2, u ^ 2, (1 + *(-2, u, xi)) ^ 2) + *(4, phi0 ^ 6, u ^ 4, (phip + *(-3, phi, u)) ^ 2) + *(8, phi0 ^ 4, u ^ 3, -1 + *(2, u, xi), phip + *(-3, phi, u)))

    ABCS[4] = *(1/3, u ^ 4, *(phi0 ^ 2, u ^ 2, *(48, xi ^ 3) + *(-1, Gp ^ 2, u ^ 3) + *(-48, B2 ^ 2, u ^ 5) + *(-16, G ^ 2, u ^ 5) + *(xi, Gp ^ 2, u ^ 4) + *(3, B2p ^ 2, u ^ 3, -1 + *(u, xi)) + *(8, G, Gp, u ^ 4) + *(16, xi, G ^ 2, u ^ 6) + *(48, xi, B2 ^ 2, u ^ 6) + *(u ^ 3, (B1p + *(-4, B1, u)) ^ 2, coshGu4sq, -1 + *(u, xi)) + *(-24, B2, B2p, u ^ 4, -1 + *(u, xi)) + *(-8, G, Gp, xi, u ^ 5)) + *(3, u ^ 3, 1 + *(u, xi), Gp ^ 2 + *(3, B2p ^ 2) + *((B1p + *(-4, B1, u)) ^ 2, coshGu4sq) + *(16, G ^ 2, u ^ 2) + *(48, B2 ^ 2, u ^ 2) + *(-24, B2, B2p, u) + *(-8, G, Gp, u)) + *(4, phi0 ^ 4, -1 + *(2, u, xi), u + *(-18, phi, u) + *(6, phip, 1 + *(u, xi)) + *(xi, u ^ 2, -3 + *(-18, phi) + *(2, u, xi))) + *(4, u, phi0 ^ 6, phip + *(-3, phi, u), *(u, 2 + *(-9, phi)) + *(3, phip, 1 + *(u, xi)) + *(xi, u ^ 2, -6 + *(-9, phi) + *(4, u, xi))) + *(4, phi0 ^ 8, u ^ 3, (phip + *(-3, phi, u)) ^ 2, -1 + *(u, xi)))

    nothing
end


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

    coshGu4sq = cosh(*(G, u ^ 4)) ^ 2


    ABCS[1] = *(6, u ^ 7)

    ABCS[2] = *(48, u ^ 6)

    ABCS[3] = *(u ^ 5, 72 + *(Gp ^ 2, u ^ 6) + *(3, B2p ^ 2, u ^ 6) + *(16, G ^ 2, u ^ 8) + *(48, B2 ^ 2, u ^ 8) + *(u ^ 6, (B1p + *(-4, B1, u)) ^ 2, coshGu4sq) + *(-24, B2, B2p, u ^ 7) + *(-8, G, Gp, u ^ 7) + *(4, phi0 ^ 2, u ^ 2, (1 + *(-2, u, xi)) ^ 2) + *(4, phi0 ^ 6, u ^ 4, (phip + *(-3, phi, u)) ^ 2) + *(8, phi0 ^ 4, u ^ 3, -1 + *(2, u, xi), phip + *(-3, phi, u)))

    ABCS[4] = *(1/3, u ^ 4, *(phi0 ^ 2, u ^ 2, *(48, xi ^ 3) + *(-1, Gp ^ 2, u ^ 3) + *(-48, B2 ^ 2, u ^ 5) + *(-16, G ^ 2, u ^ 5) + *(xi, Gp ^ 2, u ^ 4) + *(3, B2p ^ 2, u ^ 3, -1 + *(u, xi)) + *(8, G, Gp, u ^ 4) + *(16, xi, G ^ 2, u ^ 6) + *(48, xi, B2 ^ 2, u ^ 6) + *(u ^ 3, (B1p + *(-4, B1, u)) ^ 2, coshGu4sq, -1 + *(u, xi)) + *(-24, B2, B2p, u ^ 4, -1 + *(u, xi)) + *(-8, G, Gp, xi, u ^ 5)) + *(3, u ^ 3, 1 + *(u, xi), Gp ^ 2 + *(3, B2p ^ 2) + *((B1p + *(-4, B1, u)) ^ 2, coshGu4sq) + *(16, G ^ 2, u ^ 2) + *(48, B2 ^ 2, u ^ 2) + *(-24, B2, B2p, u) + *(-8, G, Gp, u)) + *(4, phi0 ^ 4, -1 + *(2, u, xi), u + *(-18, phi, u) + *(6, phip, 1 + *(u, xi)) + *(xi, u ^ 2, -3 + *(-18, phi) + *(2, u, xi))) + *(4, u, phi0 ^ 6, phip + *(-3, phi, u), *(u, 2 + *(-9, phi)) + *(3, phip, 1 + *(u, xi)) + *(xi, u ^ 2, -6 + *(-9, phi) + *(4, u, xi))) + *(4, phi0 ^ 8, u ^ 3, (phip + *(-3, phi, u)) ^ 2, -1 + *(u, xi)))

    nothing
end


# this is a coupled equation for Fx and Fy. the notation used is
#
# ( A11 d_uu Fx + A12 d_uu Fy + B11 d_u Fx + B12 d_u Fy + C11 Fx + C12 Fy ) = -S1
# ( A21 d_uu Fx + A22 d_uu Fy + B21 d_u Fx + B22 d_u Fy + C21 Fx + C22 Fy ) = -S2

function Fxy_eq_coeff!(AA::Matrix, BB::Matrix, CC::Matrix, SS::Vector, vars::FxyVars{Inner})
    u    = vars.u
    phi0 = vars.phi0

    xi    = vars.xi
    xi_x  = vars.xi_x
    xi_y  = vars.xi_y

    B1    = vars.B1
    B1p   = vars.B1p
    # B1_x  = vars.B1_x
    # B1_y  = vars.B1_y
    # B1pp  = vars.B1pp
    # B1p_x = vars.B1p_x
    # B1p_y = vars.B1p_y

    B2    = vars.B2
    B2p   = vars.B2p
    # B2_x  = vars.B2_x
    # B2_y  = vars.B2_y
    # B2pp  = vars.B2pp
    # B2p_x = vars.B2p_x
    # B2p_y = vars.B2p_y

    G     = vars.G
    Gp    = vars.Gp
    # G_x   = vars.G_x
    # G_y   = vars.G_y
    # Gpp   = vars.Gpp
    # Gp_x  = vars.Gp_x
    # Gp_y  = vars.Gp_y

    phi   = vars.phi
    phip  = vars.phip
    # phi_x = vars.phi_x
    # phi_y = vars.phi_y

    S     = vars.S
    Sp    = vars.Sp
    # S_x   = vars.S_x
    # S_y   = vars.S_y
    # Spp   = vars.Spp
    # Sp_x  = vars.Sp_x
    # Sp_y  = vars.Sp_y

    expB1u4  = exp(*(B1, u ^ 4))
    cosh2Gu4 = cosh(*(2, G, u ^ 4))
    sinh2Gu4 = sinh(*(2, G, u ^ 4))

    coshGu4sq = cosh(*(G, u ^ 4)) ^ 2

    AA[1,1] = *(2/9, u ^ 4, (3 + *(3, S, u ^ 4) + *(3, u, xi) + *(phi0 ^ 2, u ^ 2, -1 + *(u, xi))) ^ 2, expB1u4)
    AA[1,2] = 0
    AA[2,1] = 0
    AA[2,2] = *(2/9, u ^ 4, (3 + *(3, S, u ^ 4) + *(3, u, xi) + *(phi0 ^ 2, u ^ 2, -1 + *(u, xi))) ^ 2)


    BB[1,1] = *(2/9, u ^ 3, 3 + *(3, S, u ^ 4) + *(3, u, xi) + *(phi0 ^ 2, u ^ 2, -1 + *(u, xi)), 15 + *(-1, phi0 ^ 2, u ^ 2) + *(-1, u ^ 3, *(3, Sp) + *(-2, xi, phi0 ^ 2)) + *(18, u, xi) + *(27, S, u ^ 4) + *(-1, u ^ 3, B2p + *(coshGu4sq, B1p + *(-4, B1, u)) + *(-4, B2, u), 3 + *(3, S, u ^ 4) + *(3, u, xi) + *(phi0 ^ 2, u ^ 2, -1 + *(u, xi))) + *(6, phi0 ^ 2, u ^ 2, -1 + *(u, xi)), expB1u4)

    BB[1,2] = *(-1/9, u ^ 6, (3 + *(3, S, u ^ 4) + *(3, u, xi) + *(phi0 ^ 2, u ^ 2, -1 + *(u, xi))) ^ 2, *(-2, Gp) + *(-1, B1p + *(-4, B1, u), sinh2Gu4) + *(8, G, u))

    BB[2,1] = *(-1/9, u ^ 6, (3 + *(3, S, u ^ 4) + *(3, u, xi) + *(phi0 ^ 2, u ^ 2, -1 + *(u, xi))) ^ 2, *(-2, Gp) + *(-1, B1p + *(-4, B1, u), sinh2Gu4) + *(8, G, u))

    BB[2,2] = *(2/9, u ^ 3, 3 + *(3, S, u ^ 4) + *(3, u, xi) + *(phi0 ^ 2, u ^ 2, -1 + *(u, xi)), 15 + *(-1, phi0 ^ 2, u ^ 2) + *(-1, u ^ 3, *(3, Sp) + *(-2, xi, phi0 ^ 2)) + *(18, u, xi) + *(27, S, u ^ 4) + *(-1, u ^ 3, B2p + *(-1, coshGu4sq, B1p + *(-4, B1, u)) + *(-4, B2, u), 3 + *(3, S, u ^ 4) + *(3, u, xi) + *(phi0 ^ 2, u ^ 2, -1 + *(u, xi))) + *(6, phi0 ^ 2, u ^ 2, -1 + *(u, xi)))


    CC[1,1] = *(2/9, *(6, (*(3, u, 1 + *(S, u ^ 4) + *(u, xi)) + *(phi0 ^ 2, u ^ 3, -1 + *(u, xi))) ^ 2) + *(u ^ 2, 3 + *(3, S, u ^ 4) + *(3, u, xi) + *(phi0 ^ 2, u ^ 2, -1 + *(u, xi)), -6 + *(-42, Sp, u ^ 3) + *(-36, B2, u ^ 4) + *(9, B2p, u ^ 3) + *(162, S, u ^ 4) + *(phi0 ^ 2, u ^ 2, -10 + *(28, u, xi) + *(3, B2p, u ^ 3, 1 + *(-2, u, xi)) + *(12, B2, u ^ 4, -1 + *(2, u, xi))) + *(-27, B2p, S, u ^ 7) + *(108, B2, S, u ^ 8) + *(-3, u ^ 3, coshGu4sq, B1p + *(-4, B1, u), -3 + *(9, S, u ^ 4) + *(phi0 ^ 2, u ^ 2, -1 + *(2, u, xi)))) + *(-4, u ^ 2, -3 + *(9, S, u ^ 4) + *(phi0 ^ 2, u ^ 2, -1 + *(2, u, xi)), -3 + *(-3, Sp, u ^ 3) + *(9, S, u ^ 4) + *(phi0 ^ 2, u ^ 2, -1 + *(2, u, xi))) + *(2, u ^ 4, (3 + *(3, S, u ^ 4) + *(3, u, xi) + *(phi0 ^ 2, u ^ 2, -1 + *(u, xi))) ^ 2, *(u, *(-3, B2p) + *(coshGu4sq, *(-1, B1p, 3 + *(2, B1, u ^ 4)) + *(2, B1, u, 7 + *(4, B1, u ^ 4))) + *(8, G ^ 2, u ^ 5) + *(14, B2, u) + *(24, B2 ^ 2, u ^ 5) + *(-6, B2, B2p, u ^ 4) + *(-2, G, Gp, u ^ 4) + *(-2, G, u ^ 4, B1p + *(-4, B1, u), sinh2Gu4)) + *(2, phi0 ^ 2, (1 + *(-2, u, xi)) ^ 2) + *(2, u, phi0 ^ 4, -1 + *(2, u, xi), phip + *(-6, phi, u)) + *(6, phi, phi0 ^ 6, u ^ 3, *(-1, phip) + *(3, phi, u))), expB1u4)

    CC[1,2] = *(1/9, u ^ 5, 3 + *(3, S, u ^ 4) + *(3, u, xi) + *(phi0 ^ 2, u ^ 2, -1 + *(u, xi)), *(*(B1p, 9 + *(18, u, xi) + *(45, S, u ^ 4) + *(phi0 ^ 2, u ^ 2, -9 + *(12, u, xi) + *(-4, B1, u ^ 4, -1 + *(u, xi))) + *(-12, B1, u ^ 4, 1 + *(S, u ^ 4) + *(u, xi))) + *(4, B1, u, -12 + *(-48, S, u ^ 4) + *(-21, u, xi) + *(phi0 ^ 2, u ^ 2, 10 + *(-13, u, xi) + *(4, B1, u ^ 4, -1 + *(u, xi))) + *(12, B1, u ^ 4, 1 + *(S, u ^ 4) + *(u, xi))), sinh2Gu4) + *(2, Gp, 9 + *(18, u, xi) + *(45, S, u ^ 4) + *(phi0 ^ 2, u ^ 2, -9 + *(12, u, xi) + *(-4, B1, u ^ 4, -1 + *(u, xi))) + *(-12, B1, u ^ 4, 1 + *(S, u ^ 4) + *(u, xi))) + *(8, G, u, -12 + *(-48, S, u ^ 4) + *(-21, u, xi) + *(phi0 ^ 2, u ^ 2, 10 + *(-13, u, xi) + *(4, B1, u ^ 4, -1 + *(u, xi))) + *(12, B1, u ^ 4, 1 + *(S, u ^ 4) + *(u, xi))) + *(8, G, u ^ 4, B1p + *(-4, B1, u), 3 + *(3, S, u ^ 4) + *(3, u, xi) + *(phi0 ^ 2, u ^ 2, -1 + *(u, xi)), cosh2Gu4))

    CC[2,1] = *(1/9, u ^ 5, 3 + *(3, S, u ^ 4) + *(3, u, xi) + *(phi0 ^ 2, u ^ 2, -1 + *(u, xi)), *(*(-1, B1p, 9 + *(18, u, xi) + *(45, S, u ^ 4) + *(phi0 ^ 2, u ^ 2, -9 + *(12, u, xi) + *(4, B1, u ^ 4, -1 + *(u, xi))) + *(12, B1, u ^ 4, 1 + *(S, u ^ 4) + *(u, xi))) + *(4, B1, u, 12 + *(21, u, xi) + *(48, S, u ^ 4) + *(phi0 ^ 2, u ^ 2, -10 + *(13, u, xi) + *(4, B1, u ^ 4, -1 + *(u, xi))) + *(12, B1, u ^ 4, 1 + *(S, u ^ 4) + *(u, xi))), sinh2Gu4) + *(2, Gp, 9 + *(18, u, xi) + *(45, S, u ^ 4) + *(phi0 ^ 2, u ^ 2, -9 + *(12, u, xi) + *(4, B1, u ^ 4, -1 + *(u, xi))) + *(12, B1, u ^ 4, 1 + *(S, u ^ 4) + *(u, xi))) + *(-8, G, u, 12 + *(21, u, xi) + *(48, S, u ^ 4) + *(phi0 ^ 2, u ^ 2, -10 + *(13, u, xi) + *(4, B1, u ^ 4, -1 + *(u, xi))) + *(12, B1, u ^ 4, 1 + *(S, u ^ 4) + *(u, xi))) + *(-8, G, u ^ 4, B1p + *(-4, B1, u), 3 + *(3, S, u ^ 4) + *(3, u, xi) + *(phi0 ^ 2, u ^ 2, -1 + *(u, xi)), cosh2Gu4), expB1u4)

    CC[2,2] = *(4/3, (*(3, u, 1 + *(S, u ^ 4) + *(u, xi)) + *(phi0 ^ 2, u ^ 3, -1 + *(u, xi))) ^ 2) + *(-8/9, u ^ 2, -3 + *(9, S, u ^ 4) + *(phi0 ^ 2, u ^ 2, -1 + *(2, u, xi)), -3 + *(-3, Sp, u ^ 3) + *(9, S, u ^ 4) + *(phi0 ^ 2, u ^ 2, -1 + *(2, u, xi))) + *(-2/9, u ^ 2, 3 + *(3, S, u ^ 4) + *(3, u, xi) + *(phi0 ^ 2, u ^ 2, -1 + *(u, xi)), 6 + *(-162, S, u ^ 4) + *(-9, B2p, u ^ 3) + *(36, B2, u ^ 4) + *(42, Sp, u ^ 3) + *(phi0 ^ 2, u ^ 2, 10 + *(-28, u, xi) + *(3, B2p, u ^ 3, -1 + *(2, u, xi)) + *(12, B2, u ^ 4, 1 + *(-2, u, xi))) + *(-108, B2, S, u ^ 8) + *(27, B2p, S, u ^ 7) + *(-3, u ^ 3, coshGu4sq, B1p + *(-4, B1, u), -3 + *(9, S, u ^ 4) + *(phi0 ^ 2, u ^ 2, -1 + *(2, u, xi)))) + *(4/9, u ^ 4, (3 + *(3, S, u ^ 4) + *(3, u, xi) + *(phi0 ^ 2, u ^ 2, -1 + *(u, xi))) ^ 2, *(u, *(-3, B2p) + *(coshGu4sq, *(B1p, 3 + *(-2, B1, u ^ 4)) + *(2, B1, u, -7 + *(4, B1, u ^ 4))) + *(8, G ^ 2, u ^ 5) + *(14, B2, u) + *(24, B2 ^ 2, u ^ 5) + *(-6, B2, B2p, u ^ 4) + *(-2, G, Gp, u ^ 4) + *(2, G, u ^ 4, B1p + *(-4, B1, u), sinh2Gu4)) + *(2, phi0 ^ 2, (1 + *(-2, u, xi)) ^ 2) + *(2, u, phi0 ^ 4, -1 + *(2, u, xi), phip + *(-6, phi, u)) + *(6, phi, phi0 ^ 6, u ^ 3, *(-1, phip) + *(3, phi, u)))

    SS[1]   = *(-1/9, (*(3, u, 1 + *(S, u ^ 4) + *(u, xi)) + *(phi0 ^ 2, u ^ 3, -1 + *(u, xi))) ^ 2, *(-2, Ghp) + *(*(-1, B1hp) + *(B1h, u, 4 + *(B1p, u ^ 3) + *(-4, B1, u ^ 4)) + *(4, u, xi_y, *(B1p, -1 + *(B1, u ^ 4)) + *(B1, u, 5 + *(-4, B1, u ^ 4))), sinh2Gu4) + *(8, Gh, u) + *(-8, B1h, G, u ^ 5) + *(-8, Gp, u, xi_y) + *(2, B1h, Gp, u ^ 4) + *(40, G, xi_y, u ^ 2) + *(-32, B1, G, xi_y, u ^ 6) + *(-2, u ^ 4, B1p + *(-4, B1, u), Gh + *(4, G, u, xi_y), cosh2Gu4) + *(8, B1, Gp, xi_y, u ^ 5)) + *(2/9, *((*(3, u, 1 + *(S, u ^ 4) + *(u, xi)) + *(phi0 ^ 2, u ^ 3, -1 + *(u, xi))) ^ 2, *(-1, B2tp) + *(-1, coshGu4sq, B1tp + *(B1t, u, -4 + *(B1p, u ^ 3) + *(-4, B1, u ^ 4)) + *(4, u, xi_x, B1p + *(B1, B1p, u ^ 4) + *(-1, B1, u, 5 + *(4, B1, u ^ 4)))) + *(4, B2t, u) + *(-1, Gp, Gt, u ^ 4) + *(-4, B2p, u, xi_x) + *(-3, B2p, B2t, u ^ 4) + *(4, G, Gt, u ^ 5) + *(4, phit, u, phi0 ^ 4) + *(12, B2, B2t, u ^ 5) + *(16, xi_x, G ^ 2, u ^ 6) + *(20, B2, xi_x, u ^ 2) + *(48, xi_x, B2 ^ 2, u ^ 6) + *(u ^ 4, Gt + *(4, G, u, xi_x), *(-1, B1p) + *(4, B1, u), sinh2Gu4) + *(-12, B2, B2p, xi_x, u ^ 5) + *(-8, phit, xi, phi0 ^ 4, u ^ 2) + *(-8, u, xi, xi_x, phi0 ^ 2) + *(-4, G, Gp, xi_x, u ^ 5) + *(-4, phip, phit, phi0 ^ 6, u ^ 2) + *(12, phi, phit, phi0 ^ 6, u ^ 3) + *(12, phi, xi_x, phi0 ^ 4, u ^ 2) + *(16, xi_x, phi0 ^ 2, u ^ 2, xi ^ 2) + *(36, xi_x, phi ^ 2, phi0 ^ 6, u ^ 4) + *(-48, phi, xi, xi_x, phi0 ^ 4, u ^ 3) + *(-12, phi, phip, xi_x, phi0 ^ 6, u ^ 3) + *(8, phip, xi, xi_x, phi0 ^ 4, u ^ 2)) + *(-4, u ^ 3, *(3, St) + *(2, xi, xi_x, phi0 ^ 2) + *(9, S, u, xi_x), -3 + *(-3, Sp, u ^ 3) + *(9, S, u ^ 4) + *(phi0 ^ 2, u ^ 2, -1 + *(2, u, xi))) + *(3, u ^ 2, 3 + *(3, S, u ^ 4) + *(3, u, xi) + *(phi0 ^ 2, u ^ 2, -1 + *(u, xi)), *(-4, Stp) + *(u, xi_x, *(-12, Sp) + *(S, *(48, u) + *(-9, B2p, u ^ 4) + *(36, B2, u ^ 5)) + *(2, xi, phi0 ^ 2, 4 + *(-1, B2p, u ^ 3) + *(4, B2, u ^ 4))) + *(3, St, u, 4 + *(-1, B2p, u ^ 3) + *(4, B2, u ^ 4)) + *(u ^ 4, coshGu4sq, *(-1, B1p) + *(4, B1, u), *(3, St) + *(2, xi, xi_x, phi0 ^ 2) + *(9, S, u, xi_x))), expB1u4) + *(-1/3, u ^ 6, *(-2, Gp) + *(-1, B1p + *(-4, B1, u), sinh2Gu4) + *(8, G, u), *(3, Sh) + *(2, xi, xi_y, phi0 ^ 2) + *(9, S, u, xi_y), 3 + *(3, S, u ^ 4) + *(3, u, xi) + *(phi0 ^ 2, u ^ 2, -1 + *(u, xi)))

    SS[2]   = *(2/9, (*(3, u, 1 + *(S, u ^ 4) + *(u, xi)) + *(phi0 ^ 2, u ^ 3, -1 + *(u, xi))) ^ 2, *(-1, B2hp) + *(B1hp, coshGu4sq) + *(u, *(4, B2h) + *(-1, coshGu4sq, *(B1h, 4 + *(B1p, u ^ 3) + *(-4, B1, u ^ 4)) + *(xi_y, *(4, B1p, -1 + *(B1, u ^ 4)) + *(4, B1, u, 5 + *(-4, B1, u ^ 4)))) + *(-4, B2p, xi_y) + *(4, phih, phi0 ^ 4) + *(-1, Gh, Gp, u ^ 3) + *(-8, xi, xi_y, phi0 ^ 2) + *(-3, B2h, B2p, u ^ 3) + *(4, G, Gh, u ^ 4) + *(12, B2, B2h, u ^ 4) + *(16, xi_y, G ^ 2, u ^ 5) + *(20, B2, u, xi_y) + *(48, xi_y, B2 ^ 2, u ^ 5) + *(u ^ 3, B1p + *(-4, B1, u), Gh + *(4, G, u, xi_y), sinh2Gu4) + *(-12, B2, B2p, xi_y, u ^ 4) + *(-8, phih, u, xi, phi0 ^ 4) + *(-4, G, Gp, xi_y, u ^ 4) + *(-4, phih, phip, u, phi0 ^ 6) + *(12, phi, phih, phi0 ^ 6, u ^ 2) + *(12, phi, u, xi_y, phi0 ^ 4) + *(16, u, xi_y, phi0 ^ 2, xi ^ 2) + *(36, xi_y, phi ^ 2, phi0 ^ 6, u ^ 3) + *(-48, phi, xi, xi_y, phi0 ^ 4, u ^ 2) + *(-12, phi, phip, xi_y, phi0 ^ 6, u ^ 2) + *(8, phip, u, xi, xi_y, phi0 ^ 4))) + *(-8/9, u ^ 3, *(3, Sh) + *(2, xi, xi_y, phi0 ^ 2) + *(9, S, u, xi_y), -3 + *(-3, Sp, u ^ 3) + *(9, S, u ^ 4) + *(phi0 ^ 2, u ^ 2, -1 + *(2, u, xi))) + *(-2/3, u ^ 2, 3 + *(3, S, u ^ 4) + *(3, u, xi) + *(phi0 ^ 2, u ^ 2, -1 + *(u, xi)), *(4, Shp) + *(u, xi_y, *(12, Sp) + *(S, *(-48, u) + *(-36, B2, u ^ 5) + *(9, B2p, u ^ 4)) + *(2, xi, phi0 ^ 2, -4 + *(B2p, u ^ 3) + *(-4, B2, u ^ 4))) + *(-3, Sh, u, 4 + *(-1, B2p, u ^ 3) + *(4, B2, u ^ 4)) + *(u ^ 4, coshGu4sq, *(-1, B1p) + *(4, B1, u), *(3, Sh) + *(2, xi, xi_y, phi0 ^ 2) + *(9, S, u, xi_y))) + *(1/3, u ^ 2, *(3 + *(3, S, u ^ 4) + *(3, u, xi) + *(phi0 ^ 2, u ^ 2, -1 + *(u, xi)), *(2, Gtp) + *(-1, B1tp + *(B1t, u, -4 + *(B1p, u ^ 3) + *(-4, B1, u ^ 4)) + *(4, u, xi_x, B1p + *(B1, B1p, u ^ 4) + *(-1, B1, u, 5 + *(4, B1, u ^ 4))), sinh2Gu4) + *(-8, Gt, u) + *(2, B1t, u ^ 4, Gp + *(-4, G, u)) + *(8, u, xi_x, Gp + *(B1, Gp, u ^ 4) + *(-1, G, u, 5 + *(4, B1, u ^ 4))) + *(-2, u ^ 4, B1p + *(-4, B1, u), Gt + *(4, G, u, xi_x), cosh2Gu4)) + *(-3, u ^ 4, *(-2, Gp) + *(B1p + *(-4, B1, u), sinh2Gu4) + *(8, G, u), *(3, St) + *(2, xi, xi_x, phi0 ^ 2) + *(9, S, u, xi_x)), 1 + *(S, u ^ 4) + *(u, xi) + *(1/3, phi0 ^ 2, u ^ 2, -1 + *(u, xi)), expB1u4)

    nothing
end

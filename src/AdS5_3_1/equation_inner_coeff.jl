
# Vf(phi)  = VV(phi)
# Vfp(phi) = ∂(VV)(phi)

# FIXME
#V(phi) = -3.0
#Vp(phi) = 0.0


# assuming
# (A d_uu + B d_u + C Id) f = -S

function S_eq_coeff!(ABCS::Vector, vars::SVars{Inner})
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
    B1_x  = vars.B1_x
    B1_y  = vars.B1_y
    B1pp  = vars.B1pp
    B1p_x = vars.B1p_x
    B1p_y = vars.B1p_y

    B2    = vars.B2
    B2p   = vars.B2p
    B2_x  = vars.B2_x
    B2_y  = vars.B2_y
    B2pp  = vars.B2pp
    B2p_x = vars.B2p_x
    B2p_y = vars.B2p_y

    G     = vars.G
    Gp    = vars.Gp
    G_x   = vars.G_x
    G_y   = vars.G_y
    Gpp   = vars.Gpp
    Gp_x  = vars.Gp_x
    Gp_y  = vars.Gp_y

    phi   = vars.phi
    phip  = vars.phip
    phi_x = vars.phi_x
    phi_y = vars.phi_y

    S     = vars.S
    Sp    = vars.Sp
    S_x   = vars.S_x
    S_y   = vars.S_y
    Spp   = vars.Spp
    Sp_x  = vars.Sp_x
    Sp_y  = vars.Sp_y

    u2 = u*u
    u3 = u*u2
    u4 = u2*u2
    u5 = u4*u
    u5 = u3*u2
    u6 = u3*u3
    u8 = u4*u4

    expB1u4  = exp(*(B1, u4))
    cosh2Gu4 = cosh(*(2, G, u4))
    sinh2Gu4 = sinh(*(2, G, u4))

    coshGu4sq = cosh(*(G, u4)) ^ 2
    sinhGu4sq = sinh(*(G, u4)) ^ 2


    AA[1,1] = *(2/9, u4, (3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi))) ^ 2, expB1u4)

    AA[1,2] = 0


    BB[1,1] = *(u2, *(4/9, u, (3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi))) ^ 2, expB1u4) + *(2/9, u, 3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi)), 9 + *(-3, Sp, u3) + *(12, B2, u4) + *(12, u, xi) + *(21, S, u4) + *(phi0 ^ 2, u2, -5 + *(B2p, u3 + *(-1, xi, u4)) + *(-4, B2, u4) + *(6, u, xi) + *(-4, B1, u4, coshGu4sq) + *(4, B2, xi, u5) + *(B1p, u3, coshGu4sq, 1 + *(-1, u, xi)) + *(4, B1, xi, u5, coshGu4sq)) + *(-3, B2p, u3, 1 + *(S, u4) + *(u, xi)) + *(12, B1, u4, coshGu4sq) + *(12, B2, S, u8) + *(12, B2, xi, u5) + *(-3, B1p, u3, coshGu4sq, 1 + *(S, u4) + *(u, xi)) + *(12, B1, S, u8, coshGu4sq) + *(12, B1, xi, u5, coshGu4sq), expB1u4))

    BB[1,2] = *(1/9, u6, (3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi))) ^ 2, *(2, Gp) + *(B1p, sinh2Gu4) + *(-4, u, *(2, G) + *(B1, sinh2Gu4)))


    CC[1,1] = *(2/9, *(6, (*(3, u, 1 + *(S, u4) + *(u, xi)) + *(phi0 ^ 2, u3, -1 + *(u, xi))) ^ 2) + *(u2, *(-4, (-3 + *(-3, Sp, u3) + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) ^ 2) + *((*(3, u, 1 + *(S, u4) + *(u, xi)) + *(phi0 ^ 2, u3, -1 + *(u, xi))) ^ 2, B2pp + *(4, (phi0 + *(u, phi0 ^ 3, *(-1, phip) + *(3, phi, u)) + *(-2, phi0, u, xi)) ^ 2) + *(u4, (Gp + *(-4, G, u)) ^ 2) + *(coshGu4sq, B1pp + *(4, u, *(-2, B1p) + *(5, B1, u))) + *(-8, B2p, u) + *(3, u4, (B2p + *(-4, B2, u)) ^ 2) + *(20, B2, u2) + *(u4, (B1p + *(-4, B1, u)) ^ 2, coshGu4sq) + *(u4, B1p + *(-4, B1, u), Gp + *(-4, G, u), sinh2Gu4)) + *(u2, 3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi)), *(12, Spp) + *(4, phi0 ^ 2, -2 + *(6, u, xi)) + *(72, u, *(-1, Sp) + *(2, S, u)) + *(-3, u, B2p + *(-4, B2, u), -3 + *(-3, Sp, u3) + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(-3, u, coshGu4sq, B1p + *(-4, B1, u), -3 + *(-3, Sp, u3) + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))))) + *(2, u2, 3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi)), -3 + *(-3, Sp, u3) + *(9, S, u4) + *(12, B2, u4) + *(phi0 ^ 2, u2, -1 + *(B2p, u3 + *(-1, xi, u4)) + *(-4, B2, u4) + *(2, u, xi) + *(-4, B1, u4, coshGu4sq) + *(4, B2, xi, u5) + *(B1p, u3, coshGu4sq, 1 + *(-1, u, xi)) + *(4, B1, xi, u5, coshGu4sq)) + *(-3, B2p, u3, 1 + *(S, u4) + *(u, xi)) + *(12, B1, u4, coshGu4sq) + *(12, B2, S, u8) + *(12, B2, xi, u5) + *(-3, B1p, u3, coshGu4sq, 1 + *(S, u4) + *(u, xi)) + *(12, B1, S, u8, coshGu4sq) + *(12, B1, xi, u5, coshGu4sq)), expB1u4)

    CC[1,2] = *(1/9, u4, *(u, *(B1p, *(21 + *(-9, Sp, u3) + *(30, u, xi) + *(57, S, u4) + *(phi0 ^ 2, u2, -13 + *(16, u, xi) + *(-8, B1, u4, -1 + *(u, xi))) + *(-24, B1, u4, 1 + *(S, u4) + *(u, xi)), sinh2Gu4) + *(16, G, u4, sinhGu4sq, 3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi)))) + *(-4, u, *(6, G, 4 + *(-3, Sp, u3) + *(7, u, xi) + *(16, S, u4) + *(8, B1, u4, sinhGu4sq, 1 + *(S, u4) + *(u, xi))) + *(-3, B1, -4 + *(-16, S, u4) + *(-7, u, xi) + *(3, Sp, u3) + *(4, B1, u4, 1 + *(S, u4) + *(u, xi)), sinh2Gu4) + *(2, G, phi0 ^ 2, u2, -10 + *(13, u, xi) + *(8, B1, u4, sinhGu4sq, -1 + *(u, xi))) + *(B1, phi0 ^ 2, u2, -10 + *(13, u, xi) + *(-4, B1, u4, -1 + *(u, xi)), sinh2Gu4)) + *(2, Gp, 21 + *(-9, Sp, u3) + *(30, u, xi) + *(57, S, u4) + *(phi0 ^ 2, u2, -13 + *(16, u, xi) + *(-2, B1p, u3, sinhGu4sq, -1 + *(u, xi)) + *(8, B1, u4, sinhGu4sq, -1 + *(u, xi))) + *(24, B1, u4, sinhGu4sq) + *(-6, B1p, u3, sinhGu4sq, 1 + *(S, u4) + *(u, xi)) + *(24, B1, S, u8, sinhGu4sq) + *(24, B1, xi, u5, sinhGu4sq)) + *(B1p ^ 2, u3, 3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi)), sinh2Gu4)) + *(-2, Gpp, 3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi))) + *(-1, B1pp, 3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi)), sinh2Gu4), 3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi)))


    SS[1] = *(1/9, (*(3, u, 1 + *(S, u4) + *(u, xi)) + *(phi0 ^ 2, u3, -1 + *(u, xi))) ^ 2, *(2, Gp_y) + *(B1p_y, sinh2Gu4) + *(-8, G_y, u) + *(-2, B2p_x, expB1u4) + *(-8, xi_x, phi0 ^ 2, expB1u4) + *(-4, B1_y, u, sinh2Gu4) + *(-2, B1_y, Gp, u4) + *(-2, B1p_x, coshGu4sq, expB1u4) + *(8, B1_y, G, u5) + *(8, B2_x, u, expB1u4) + *(-1, B1_y, B1p, u4, sinh2Gu4) + *(-8, B1, G_y, u5, cosh2Gu4) + *(-6, B2_x, B2p, u4, expB1u4) + *(-2, G_x, Gp, u4, expB1u4) + *(2, B1p, G_y, u4, cosh2Gu4) + *(4, B1, B1_y, u5, sinh2Gu4) + *(8, B1_x, u, coshGu4sq, expB1u4) + *(8, G, G_x, u5, expB1u4) + *(8, phi_x, u, phi0 ^ 4, expB1u4) + *(24, B2, B2_x, u5, expB1u4) + *(-24, phi, xi_x, phi0 ^ 4, u2, expB1u4) + *(-16, phi_x, xi, phi0 ^ 4, u2, expB1u4) + *(-8, phi_x, phip, phi0 ^ 6, u2, expB1u4) + *(-2, B1_x, B1p, u4, coshGu4sq, expB1u4) + *(-2, B1p, G_x, u4, expB1u4, sinh2Gu4) + *(8, B1, B1_x, u5, coshGu4sq, expB1u4) + *(8, B1, G_x, u5, expB1u4, sinh2Gu4) + *(8, phip, u, xi_x, phi0 ^ 4, expB1u4) + *(16, u, xi, xi_x, phi0 ^ 2, expB1u4) + *(24, phi, phi_x, phi0 ^ 6, u3, expB1u4)) + *(1/9, u2, *(-2, *(12, Sp_x) + *(xi_x, phi0 ^ 2, -8 + *(-12, B2, u4) + *(3, B2p, u3) + *(-12, B1, u4, coshGu4sq) + *(3, B1p, u3, coshGu4sq)) + *(-9, S_x, u, 4 + *(-1, B2p, u3) + *(4, B2, u4) + *(-1, B1p, u3, coshGu4sq) + *(4, B1, u4, coshGu4sq)) + *(9, u, xi_x, B2p + *(B1p, coshGu4sq) + *(-4, u, B2 + *(B1, coshGu4sq))), expB1u4) + *(-3, u, *(xi_y, 3 + *(phi0 ^ 2, u2)) + *(3, S_y, u3), *(-2, Gp) + *(-1, B1p, sinh2Gu4) + *(8, G, u) + *(4, B1, u, sinh2Gu4)), 3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi))) + *(2/9, xi_x, *(-4, (-3 + *(-3, Sp, u3) + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) ^ 2) + *((*(3, u, 1 + *(S, u4) + *(u, xi)) + *(phi0 ^ 2, u3, -1 + *(u, xi))) ^ 2, B2pp + *(4, (phi0 + *(u, phi0 ^ 3, *(-1, phip) + *(3, phi, u)) + *(-2, phi0, u, xi)) ^ 2) + *(u4, (Gp + *(-4, G, u)) ^ 2) + *(coshGu4sq, B1pp + *(4, u, *(-2, B1p) + *(5, B1, u))) + *(-8, B2p, u) + *(3, u4, (B2p + *(-4, B2, u)) ^ 2) + *(20, B2, u2) + *(u4, (B1p + *(-4, B1, u)) ^ 2, coshGu4sq) + *(u4, B1p + *(-4, B1, u), Gp + *(-4, G, u), sinh2Gu4)) + *(u2, 3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi)), *(12, Spp) + *(4, phi0 ^ 2, -2 + *(6, u, xi)) + *(72, u, *(-1, Sp) + *(2, S, u)) + *(-3, u, B2p + *(-4, B2, u), -3 + *(-3, Sp, u3) + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(-3, u, coshGu4sq, B1p + *(-4, B1, u), -3 + *(-3, Sp, u3) + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi)))), expB1u4) + *(8/9, *(xi_x, 3 + *(phi0 ^ 2, u2)) + *(3, S_x, u3), 3 + *(phi0 ^ 2, u2 + *(-2, xi, u3)) + *(-9, S, u4) + *(3, Sp, u3), expB1u4) + *(-1/9, xi_y, u2, *(-1, u, *(B1p, *(15 + *(-9, Sp, u3) + *(24, u, xi) + *(51, S, u4) + *(phi0 ^ 2, u2, -11 + *(14, u, xi) + *(-8, B1, u4, -1 + *(u, xi))) + *(-24, B1, u4, 1 + *(S, u4) + *(u, xi)), sinh2Gu4) + *(16, G, u4, sinhGu4sq, 3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi)))) + *(-4, u, *(6, G, 2 + *(-3, Sp, u3) + *(5, u, xi) + *(14, S, u4) + *(8, B1, u4, sinhGu4sq, 1 + *(S, u4) + *(u, xi))) + *(-3, B1, -2 + *(-14, S, u4) + *(-5, u, xi) + *(3, Sp, u3) + *(4, B1, u4, 1 + *(S, u4) + *(u, xi)), sinh2Gu4) + *(2, G, phi0 ^ 2, u2, -8 + *(11, u, xi) + *(8, B1, u4, sinhGu4sq, -1 + *(u, xi))) + *(B1, phi0 ^ 2, u2, -8 + *(11, u, xi) + *(-4, B1, u4, -1 + *(u, xi)), sinh2Gu4)) + *(2, Gp, 15 + *(-9, Sp, u3) + *(24, u, xi) + *(51, S, u4) + *(phi0 ^ 2, u2, -11 + *(14, u, xi) + *(-2, B1p, u3, sinhGu4sq, -1 + *(u, xi)) + *(8, B1, u4, sinhGu4sq, -1 + *(u, xi))) + *(24, B1, u4, sinhGu4sq) + *(-6, B1p, u3, sinhGu4sq, 1 + *(S, u4) + *(u, xi)) + *(24, B1, S, u8, sinhGu4sq) + *(24, B1, xi, u5, sinhGu4sq)) + *(B1p ^ 2, u3, 3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi)), sinh2Gu4)) + *(2, Gpp, 3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi))) + *(B1pp, 3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi)), sinh2Gu4), 3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi)))


    AA[2,1] = 0

    AA[2,2] = *(2/9, u4, (3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi))) ^ 2)


    BB[2,1] = *(1/9, u6, (3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi))) ^ 2, *(2, Gp) + *(-1, B1p + *(-4, B1, u), sinh2Gu4) + *(-8, G, u), expB1u4)

    BB[2,2] = *(2/9, u3, 3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi)), 15 + *(-3, Sp, u3) + *(12, B2, u4) + *(18, u, xi) + *(27, S, u4) + *(phi0 ^ 2, u2, -7 + *(B2p, u3 + *(-1, xi, u4)) + *(8, u, xi) + *(4, B2, u4, -1 + *(u, xi))) + *(-3, B2p, u3, 1 + *(S, u4) + *(u, xi)) + *(12, B2, S, u8) + *(12, B2, xi, u5) + *(u3, coshGu4sq, B1p + *(-4, B1, u), 3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi))))


    CC[2,1] = *(1/9, u4, 3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi)), *(*(3, B1pp, 1 + *(S, u4) + *(u, xi)) + *(B1p ^ 2, u4, 3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi))) + *(-1, B1p, u, 21 + *(-9, Sp, u3) + *(30, u, xi) + *(57, S, u4) + *(phi0 ^ 2, u2, -13 + *(16, u, xi) + *(8, B1, u4, -1 + *(u, xi))) + *(24, B1, u4, 1 + *(S, u4) + *(u, xi))) + *(4, B1, u2, 12 + *(-9, Sp, u3) + *(21, u, xi) + *(48, S, u4) + *(phi0 ^ 2, u2, -10 + *(13, u, xi) + *(4, B1, u4, -1 + *(u, xi))) + *(12, B1, u4, 1 + *(S, u4) + *(u, xi))) + *(B1pp, phi0 ^ 2, u2, -1 + *(u, xi)), sinh2Gu4) + *(-6, Gpp, 1 + *(S, u4) + *(u, xi)) + *(-8, G, u2, 12 + *(-9, Sp, u3) + *(21, u, xi) + *(48, S, u4) + *(phi0 ^ 2, u2, -10 + *(13, u, xi))) + *(2, Gp, u, 21 + *(-9, Sp, u3) + *(30, u, xi) + *(57, S, u4) + *(phi0 ^ 2, u2, -13 + *(16, u, xi))) + *(-2, Gpp, phi0 ^ 2, u2, -1 + *(u, xi)) + *(4, u4, sinhGu4sq, B1p + *(-4, B1, u), Gp + *(-4, G, u), 3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi))), expB1u4)

    CC[2,2] = *(4/3, (*(3, u, 1 + *(S, u4) + *(u, xi)) + *(phi0 ^ 2, u3, -1 + *(u, xi))) ^ 2) + *(2/9, u2, *(-4, (-3 + *(-3, Sp, u3) + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) ^ 2) + *((*(3, u, 1 + *(S, u4) + *(u, xi)) + *(phi0 ^ 2, u3, -1 + *(u, xi))) ^ 2, B2pp + *(4, (phi0 + *(u, phi0 ^ 3, *(-1, phip) + *(3, phi, u)) + *(-2, phi0, u, xi)) ^ 2) + *(u4, (Gp + *(-4, G, u)) ^ 2) + *(-1, coshGu4sq, B1pp + *(4, u, *(-2, B1p) + *(5, B1, u))) + *(-8, B2p, u) + *(3, u4, (B2p + *(-4, B2, u)) ^ 2) + *(20, B2, u2) + *(u4, (B1p + *(-4, B1, u)) ^ 2, coshGu4sq) + *(u4, B1p + *(-4, B1, u), *(-1, Gp) + *(4, G, u), sinh2Gu4)) + *(u2, 3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi)), *(12, Spp) + *(4, phi0 ^ 2, -2 + *(6, u, xi)) + *(72, u, *(-1, Sp) + *(2, S, u)) + *(-3, u, B2p + *(-4, B2, u), -3 + *(-3, Sp, u3) + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(3, u, coshGu4sq, B1p + *(-4, B1, u), -3 + *(-3, Sp, u3) + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))))) + *(-4/9, u2, 3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi)), 3 + *(phi0 ^ 2, u2 + *(u6, *(4, B2) + *(B2p, xi)) + *(-1, B2p, u5) + *(-2, xi, u3) + *(-4, B2, xi, u ^ 7)) + *(-12, B2, u4) + *(-9, S, u4) + *(3, Sp, u3) + *(-12, B2, S, u8) + *(-12, B2, xi, u5) + *(3, B2p, u3, 1 + *(S, u4) + *(u, xi)) + *(u3, coshGu4sq, *(-1, B1p) + *(4, B1, u), 3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi))))


    SS[2] = *(1/9, (*(3, u, 1 + *(S, u4) + *(u, xi)) + *(phi0 ^ 2, u3, -1 + *(u, xi))) ^ 2, *(-2, B2p_y) + *(-8, xi_y, phi0 ^ 2) + *(2, B1p_y, coshGu4sq) + *(2, Gp_x, expB1u4) + *(8, B2_y, u) + *(-1, B1p_x, expB1u4, sinh2Gu4) + *(-8, B1_y, u, coshGu4sq) + *(-8, G_x, u, expB1u4) + *(-6, B2_y, B2p, u4) + *(-2, G_y, Gp, u4) + *(8, G, G_y, u5) + *(8, phi_y, u, phi0 ^ 4) + *(24, B2, B2_y, u5) + *(-24, phi, xi_y, phi0 ^ 4, u2) + *(-16, phi_y, xi, phi0 ^ 4, u2) + *(-8, B1, G_y, u5, sinh2Gu4) + *(-8, B1_x, G, u5, expB1u4) + *(-8, phi_y, phip, phi0 ^ 6, u2) + *(-2, B1_y, B1p, u4, coshGu4sq) + *(2, B1_x, Gp, u4, expB1u4) + *(2, B1p, G_y, u4, sinh2Gu4) + *(4, B1_x, u, expB1u4, sinh2Gu4) + *(8, B1, B1_y, u5, coshGu4sq) + *(8, phip, u, xi_y, phi0 ^ 4) + *(16, u, xi, xi_y, phi0 ^ 2) + *(24, phi, phi_y, phi0 ^ 6, u3) + *(-1, B1_x, B1p, u4, expB1u4, sinh2Gu4) + *(-2, B1p, G_x, u4, cosh2Gu4, expB1u4) + *(4, B1, B1_x, u5, expB1u4, sinh2Gu4) + *(8, B1, G_x, u5, cosh2Gu4, expB1u4)) + *(1/9, *(xi_y, 3 + *(phi0 ^ 2, u2)) + *(3, S_y, u3), 24 + *(-72, S, u4) + *(8, phi0 ^ 2, u2 + *(-2, xi, u3)) + *(24, Sp, u3)) + *(2/9, xi_y, *(-4, (-3 + *(-3, Sp, u3) + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) ^ 2) + *((*(3, u, 1 + *(S, u4) + *(u, xi)) + *(phi0 ^ 2, u3, -1 + *(u, xi))) ^ 2, B2pp + *(4, (phi0 + *(u, phi0 ^ 3, *(-1, phip) + *(3, phi, u)) + *(-2, phi0, u, xi)) ^ 2) + *(u4, (Gp + *(-4, G, u)) ^ 2) + *(-1, coshGu4sq, B1pp + *(4, u, *(-2, B1p) + *(5, B1, u))) + *(-8, B2p, u) + *(3, u4, (B2p + *(-4, B2, u)) ^ 2) + *(20, B2, u2) + *(u4, (B1p + *(-4, B1, u)) ^ 2, coshGu4sq) + *(u4, B1p + *(-4, B1, u), *(-1, Gp) + *(4, G, u), sinh2Gu4)) + *(u2, 3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi)), *(12, Spp) + *(4, phi0 ^ 2, -2 + *(6, u, xi)) + *(72, u, *(-1, Sp) + *(2, S, u)) + *(-3, u, B2p + *(-4, B2, u), -3 + *(-3, Sp, u3) + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(3, u, coshGu4sq, B1p + *(-4, B1, u), -3 + *(-3, Sp, u3) + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))))) + *(1/3, u2, 1 + *(S, u4) + *(u, xi) + *(1/3, phi0 ^ 2, u2, -1 + *(u, xi)), *(-24, Sp_y) + *(16, xi_y, phi0 ^ 2) + *(72, S_y, u) + *(-6, u, B2p + *(-4, B2, u), *(xi_y, 3 + *(phi0 ^ 2, u2)) + *(3, S_y, u3)) + *(6, u, coshGu4sq, B1p + *(-4, B1, u), *(xi_y, 3 + *(phi0 ^ 2, u2)) + *(3, S_y, u3)) + *(6, u, Gp + *(-4, G, u), *(xi_x, 3 + *(phi0 ^ 2, u2)) + *(3, S_x, u3), expB1u4) + *(3, u, *(-1, B1p) + *(4, B1, u), *(xi_x, 3 + *(phi0 ^ 2, u2)) + *(3, S_x, u3), expB1u4, sinh2Gu4)) + *(1/9, xi_x, u2, *(u, *(-1, B1p, *(15 + *(-9, Sp, u3) + *(24, u, xi) + *(51, S, u4) + *(phi0 ^ 2, u2, -11 + *(14, u, xi) + *(8, B1, u4, -1 + *(u, xi))) + *(24, B1, u4, 1 + *(S, u4) + *(u, xi)), sinh2Gu4) + *(16, G, u4, sinhGu4sq, 3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi)))) + *(4, u, *(6, G, -2 + *(-14, S, u4) + *(-5, u, xi) + *(3, Sp, u3) + *(8, B1, u4, sinhGu4sq, 1 + *(S, u4) + *(u, xi))) + *(3, B1, 2 + *(-3, Sp, u3) + *(5, u, xi) + *(14, S, u4) + *(4, B1, u4, 1 + *(S, u4) + *(u, xi)), sinh2Gu4) + *(2, G, phi0 ^ 2, u2, 8 + *(-11, u, xi) + *(8, B1, u4, sinhGu4sq, -1 + *(u, xi))) + *(B1, phi0 ^ 2, u2, -8 + *(11, u, xi) + *(4, B1, u4, -1 + *(u, xi)), sinh2Gu4)) + *(6, Gp, 5 + *(-3, Sp, u3) + *(8, u, xi) + *(17, S, u4) + *(-8, B1, u4, sinhGu4sq) + *(-8, B1, S, u8, sinhGu4sq) + *(-8, B1, xi, u5, sinhGu4sq) + *(2, B1p, u3, sinhGu4sq, 1 + *(S, u4) + *(u, xi))) + *(B1p ^ 2, u3, 3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi)), sinh2Gu4) + *(-2, Gp, phi0 ^ 2, u2, 11 + *(-14, u, xi) + *(-2, B1p, u3, sinhGu4sq, -1 + *(u, xi)) + *(8, B1, u4, sinhGu4sq, -1 + *(u, xi)))) + *(-2, Gpp, 3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi))) + *(B1pp, 3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi)), sinh2Gu4), 3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi)), expB1u4)

    nothing
end


function Sd_eq_coeff!(ABCS::Vector, vars::AllVars{Inner})
    u      = vars.u
    phi0   = vars.phi0

    xi     = vars.xi
    xi_x   = vars.xi_x
    xi_y   = vars.xi_y
    xi_xx  = vars.xi_xx
    xi_xy  = vars.xi_xy
    xi_yy  = vars.xi_yy

    B1     = vars.B1
    B1p    = vars.B1p
    B1t    = vars.B1t
    B1h    = vars.B1h
    B1b    = vars.B1b
    B1s    = vars.B1s
    B1pt   = vars.B1pt
    B1ph   = vars.B1ph

    B2     = vars.B2
    B2p    = vars.B2p
    B2t    = vars.B2t
    B2h    = vars.B2h
    B2b    = vars.B2b
    B2s    = vars.B2s
    B2pt   = vars.B2pt
    B2ph   = vars.B2ph

    G      = vars.G
    Gp     = vars.Gp
    Gt     = vars.Gt
    Gh     = vars.Gh
    Gb     = vars.Gb
    Gs     = vars.Gs
    Gpt    = vars.Gpt
    Gph    = vars.Gph

    phi    = vars.phi
    phip   = vars.phip
    phit   = vars.phit
    phih   = vars.phih
    phib   = vars.phib
    phis   = vars.phis
    phipt  = vars.phipt
    phiph  = vars.phiph

    S      = vars.S
    Sp     = vars.Sp
    St     = vars.St
    Sh     = vars.Sh
    Sb     = vars.Sb
    Ss     = vars.Ss
    Spt    = vars.Spt
    Sph    = vars.Sph

    Fx     = vars.Fx
    Fxp    = vars.Fxp
    Fxt    = vars.Fxt
    Fxh    = vars.Fxh
    # Fxb   = vars.Fxb
    # Fxs   = vars.Fxs
    Fxpt   = vars.Fxpt
    Fxph   = vars.Fxph

    Fy     = vars.Fy
    Fyp    = vars.Fyp
    Fyt    = vars.Fyt
    Fyh    = vars.Fyh
    # Fyb   = vars.Fyb
    # Fys   = vars.Fys
    Fypt   = vars.Fypt
    Fyph   = vars.Fyph

    B2c   = vars.B2c
    Gc    = vars.Gc
    Sc    = vars.Sc

    Spp   = vars.Spp

    u2 = u*u
    u3 = u*u2
    u4 = u2*u2
    u5 = u4*u
    u5 = u3*u2
    u6 = u3*u3
    u8 = u4*u4

    coshGu4 = cosh(*(G, u4))
    sinhGu4 = sinh(*(G, u4))

    expB1u4 = exp(*(B1, u ^ 4))
    expB2u4 = exp(*(B2, u ^ 4))


    ABCS[1] = 0

    ABCS[2] = *(-4/9, u2, (3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi))) ^ 3, expB1u4)

    ABCS[3] = *(-8/9, (*(3, u, 1 + *(S, u4) + *(u, xi)) + *(phi0 ^ 2, u3, -1 + *(u, xi))) ^ 2, *(3, xi) + *(-3, Sp, u2) + *(12, S, u3) + *(u, phi0 ^ 2, -2 + *(3, u, xi)), expB1u4)

    ABCS[4] = *(u, *(-4, (xi_x + *(u3, St + *(Sp, xi_x + *(Fx, u2))) + *(-1, xi_x + *(Fx, u2), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3)) + *(1/3, xi_x, phi0 ^ 2, u2)) ^ 2, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-4, (xi_y + *(u3, Sh + *(Sp, xi_y + *(Fy, u2))) + *(-1, xi_y + *(Fy, u2), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3)) + *(1/3, xi_y, phi0 ^ 2, u2)) ^ 2, coshGu4, expB2u4) + *(8, xi_x + *(u3, St + *(Sp, xi_x + *(Fx, u2))) + *(-1, xi_x + *(Fx, u2), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3)) + *(1/3, xi_x, phi0 ^ 2, u2), xi_y + *(u3, Sh + *(Sp, xi_y + *(Fy, u2))) + *(-1, xi_y + *(Fy, u2), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3)) + *(1/3, xi_y, phi0 ^ 2, u2), exp(*(B1, u4) + *(B2, u4)), sinhGu4)) + *(1 + *(S, u4) + *(u, xi) + *(-1/3, phi0 ^ 2, u2) + *(1/3, xi, phi0 ^ 2, u3), *(-16, xi_xy + *(u3, Sc + *(Spt, xi_y + *(Fy, u2)) + *(Sph + *(Spp, xi_y + *(Fy, u2)), xi_x + *(Fx, u2))) + *(-1, xi_x + *(Fx, u2), *(u3, Sph + *(Spp, xi_y + *(Fy, u2))) + *(-3, u4, Sh + *(Sp, xi_y + *(Fy, u2))) + *(-2/3, xi_y, phi0 ^ 2, u3)) + *(-1, xi_y + *(Fy, u2), *(u3, Spt + *(Spp, xi_x + *(Fx, u2))) + *(-1, xi_x + *(Fx, u2), *(Spp, u3) + *(-6, Sp, u4) + *(12, S, u5) + *(-2/3, phi0 ^ 2, u3) + *(2, xi, phi0 ^ 2, u4)) + *(-3, u4, St + *(Sp, xi_x + *(Fx, u2))) + *(-2/3, xi_x, phi0 ^ 2, u3)) + *(1/3, xi_xy, phi0 ^ 2, u2), exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(8, xi_xx + *(u3, Sb + *(Spt, xi_x + *(Fx, u2)) + *(Spt + *(Spp, xi_x + *(Fx, u2)), xi_x + *(Fx, u2))) + *(-1, xi_x + *(Fx, u2), *(u3, Spt + *(Spp, xi_x + *(Fx, u2))) + *(-3, u4, St + *(Sp, xi_x + *(Fx, u2))) + *(-2/3, xi_x, phi0 ^ 2, u3)) + *(-1, xi_x + *(Fx, u2), *(u3, Spt + *(Spp, xi_x + *(Fx, u2))) + *(-1, xi_x + *(Fx, u2), *(Spp, u3) + *(-6, Sp, u4) + *(12, S, u5) + *(-2/3, phi0 ^ 2, u3) + *(2, xi, phi0 ^ 2, u4)) + *(-3, u4, St + *(Sp, xi_x + *(Fx, u2))) + *(-2/3, xi_x, phi0 ^ 2, u3)) + *(1/3, xi_xx, phi0 ^ 2, u2), coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(8, xi_yy + *(u3, Ss + *(Sph, xi_y + *(Fy, u2)) + *(Sph + *(Spp, xi_y + *(Fy, u2)), xi_y + *(Fy, u2))) + *(-1, xi_y + *(Fy, u2), *(u3, Sph + *(Spp, xi_y + *(Fy, u2))) + *(-3, u4, Sh + *(Sp, xi_y + *(Fy, u2))) + *(-2/3, xi_y, phi0 ^ 2, u3)) + *(-1, xi_y + *(Fy, u2), *(u3, Sph + *(Spp, xi_y + *(Fy, u2))) + *(-1, xi_y + *(Fy, u2), *(Spp, u3) + *(-6, Sp, u4) + *(12, S, u5) + *(-2/3, phi0 ^ 2, u3) + *(2, xi, phi0 ^ 2, u4)) + *(-3, u4, Sh + *(Sp, xi_y + *(Fy, u2))) + *(-2/3, xi_y, phi0 ^ 2, u3)) + *(1/3, xi_yy, phi0 ^ 2, u2), coshGu4, expB2u4) + *(-8, *(u4, B1h + *(B1p, xi_y + *(Fy, u2))) + *(-1, xi_y + *(Fy, u2), *(B1p, u4) + *(-4, B1, u5)), xi_y + *(u3, Sh + *(Sp, xi_y + *(Fy, u2))) + *(-1, xi_y + *(Fy, u2), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3)) + *(1/3, xi_y, phi0 ^ 2, u2), coshGu4, expB2u4) + *(-8, *(u4, B2h + *(B2p, xi_y + *(Fy, u2))) + *(-1, xi_y + *(Fy, u2), *(B2p, u4) + *(-4, B2, u5)), xi_x + *(u3, St + *(Sp, xi_x + *(Fx, u2))) + *(-1, xi_x + *(Fx, u2), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3)) + *(1/3, xi_x, phi0 ^ 2, u2), exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-8, *(u4, B2t + *(B2p, xi_x + *(Fx, u2))) + *(-1, xi_x + *(Fx, u2), *(B2p, u4) + *(-4, B2, u5)), xi_y + *(u3, Sh + *(Sp, xi_y + *(Fy, u2))) + *(-1, xi_y + *(Fy, u2), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3)) + *(1/3, xi_y, phi0 ^ 2, u2), exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-8, *(u4, Gh + *(Gp, xi_y + *(Fy, u2))) + *(-1, xi_y + *(Fy, u2), *(Gp, u4) + *(-4, G, u5)), xi_x + *(u3, St + *(Sp, xi_x + *(Fx, u2))) + *(-1, xi_x + *(Fx, u2), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3)) + *(1/3, xi_x, phi0 ^ 2, u2), coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-8, *(u4, Gt + *(Gp, xi_x + *(Fx, u2))) + *(-1, xi_x + *(Fx, u2), *(Gp, u4) + *(-4, G, u5)), xi_y + *(u3, Sh + *(Sp, xi_y + *(Fy, u2))) + *(-1, xi_y + *(Fy, u2), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3)) + *(1/3, xi_y, phi0 ^ 2, u2), coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-8, xi_xx + *(u2, Fxt + *(Fxp, xi_x + *(Fx, u2))) + *(-1, xi_x + *(Fx, u2), *(Fxp, u2) + *(-2, Fx, u3)), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3), coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-8, xi_yy + *(u2, Fyh + *(Fyp, xi_y + *(Fy, u2))) + *(-1, xi_y + *(Fy, u2), *(Fyp, u2) + *(-2, Fy, u3)), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3), coshGu4, expB2u4) + *(-2, *(Fxp, u2) + *(-2, Fx, u3), xi_x + *(u3, St + *(Sp, xi_x + *(Fx, u2))) + *(-1, xi_x + *(Fx, u2), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3)) + *(1/3, xi_x, phi0 ^ 2, u2), coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-2, *(Fyp, u2) + *(-2, Fy, u3), xi_y + *(u3, Sh + *(Sp, xi_y + *(Fy, u2))) + *(-1, xi_y + *(Fy, u2), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3)) + *(1/3, xi_y, phi0 ^ 2, u2), coshGu4, expB2u4) + *(2, *(Fxp, u2) + *(-2, Fx, u3), xi_y + *(u3, Sh + *(Sp, xi_y + *(Fy, u2))) + *(-1, xi_y + *(Fy, u2), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3)) + *(1/3, xi_y, phi0 ^ 2, u2), exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(2, *(Fyp, u2) + *(-2, Fy, u3), xi_x + *(u3, St + *(Sp, xi_x + *(Fx, u2))) + *(-1, xi_x + *(Fx, u2), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3)) + *(1/3, xi_x, phi0 ^ 2, u2), exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(8, *(u4, B1t + *(B1p, xi_x + *(Fx, u2))) + *(-1, xi_x + *(Fx, u2), *(B1p, u4) + *(-4, B1, u5)), xi_x + *(u3, St + *(Sp, xi_x + *(Fx, u2))) + *(-1, xi_x + *(Fx, u2), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3)) + *(1/3, xi_x, phi0 ^ 2, u2), coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(8, *(u4, B2h + *(B2p, xi_y + *(Fy, u2))) + *(-1, xi_y + *(Fy, u2), *(B2p, u4) + *(-4, B2, u5)), xi_y + *(u3, Sh + *(Sp, xi_y + *(Fy, u2))) + *(-1, xi_y + *(Fy, u2), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3)) + *(1/3, xi_y, phi0 ^ 2, u2), coshGu4, expB2u4) + *(8, *(u4, B2t + *(B2p, xi_x + *(Fx, u2))) + *(-1, xi_x + *(Fx, u2), *(B2p, u4) + *(-4, B2, u5)), xi_x + *(u3, St + *(Sp, xi_x + *(Fx, u2))) + *(-1, xi_x + *(Fx, u2), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3)) + *(1/3, xi_x, phi0 ^ 2, u2), coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(8, *(u4, Gh + *(Gp, xi_y + *(Fy, u2))) + *(-1, xi_y + *(Fy, u2), *(Gp, u4) + *(-4, G, u5)), xi_y + *(u3, Sh + *(Sp, xi_y + *(Fy, u2))) + *(-1, xi_y + *(Fy, u2), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3)) + *(1/3, xi_y, phi0 ^ 2, u2), expB2u4, sinhGu4) + *(8, *(u4, Gt + *(Gp, xi_x + *(Fx, u2))) + *(-1, xi_x + *(Fx, u2), *(Gp, u4) + *(-4, G, u5)), xi_x + *(u3, St + *(Sp, xi_x + *(Fx, u2))) + *(-1, xi_x + *(Fx, u2), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3)) + *(1/3, xi_x, phi0 ^ 2, u2), exp(*(B2, u4) + *(2, B1, u4)), sinhGu4) + *(8, xi_xy + *(u2, Fxh + *(Fxp, xi_y + *(Fy, u2))) + *(-1, xi_y + *(Fy, u2), *(Fxp, u2) + *(-2, Fx, u3)), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3), exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(8, xi_xy + *(u2, Fyt + *(Fyp, xi_x + *(Fx, u2))) + *(-1, xi_x + *(Fx, u2), *(Fyp, u2) + *(-2, Fy, u3)), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3), exp(*(B1, u4) + *(B2, u4)), sinhGu4)) + *(-4/27, *(-324, xi ^ 3) + *(phi0 ^ 8, u5) + *(-81, u, xi ^ 4) + *(-6, phi0 ^ 6, u3) + *(108, xi, phi0 ^ 2) + *(-81, S ^ 3, u ^ 9, *(-3, (1 + *(u, xi)) ^ 2) + *(phi0 ^ 2, u2)) + *(-63, phi0 ^ 4, u3, xi ^ 2) + *(-54, phi0 ^ 4, u4, xi ^ 3) + *(-24, phi0 ^ 6, u6, xi ^ 3) + *(-9, S ^ 2, u5, *(-9, (1 + *(u, xi)) ^ 2, 5 + *(6, u, xi)) + *(phi0 ^ 4, u4, -7 + *(8, u, xi)) + *(-3, phi0 ^ 2, u2, -12 + *(-12, u, xi) + *(8, u3, xi ^ 3) + *(9, u2, xi ^ 2))) + *(-4, xi, phi0 ^ 8, u6) + *(-3, S, u, *(-27, (1 + *(u, xi)) ^ 3, 1 + *(3, u, xi)) + *(phi0 ^ 6, u6, 5 + *(-12, u, xi) + *(7, u2, xi ^ 2)) + *(-9, phi0 ^ 2, u2, -7 + *(u2, xi ^ 2) + *(-16, u, xi) + *(10, u4, xi ^ 4) + *(20, u3, xi ^ 3)) + *(-3, phi0 ^ 4, u4, 11 + *(-22, u2, xi ^ 2) + *(-2, u, xi) + *(2, u3, xi ^ 3) + *(7, u4, xi ^ 4))) + *(-3, phi0 ^ 6, u ^ 7, xi ^ 4) + *(-2, phi0 ^ 8, u8, xi ^ 3) + *(3, Sp, (3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi))) ^ 2, *(-3, (1 + *(u, xi)) ^ 2) + *(phi0 ^ 2, u2)) + *(5, phi0 ^ 8, u ^ 7, xi ^ 2) + *(6, phi0 ^ 6, u8, xi ^ 5) + *(12, xi, phi0 ^ 6, u4) + *(15, phi0 ^ 6, u5, xi ^ 2) + *(36, phi0 ^ 4, u6, xi ^ 5) + *(45, phi0 ^ 4, u5, xi ^ 4) + *(54, phi0 ^ 2, u4, xi ^ 5) + *(81, u, phi0 ^ 2, xi ^ 2) + *(108, phi0 ^ 2, u2, xi ^ 3) + *(135, phi0 ^ 2, u3, xi ^ 4), expB1u4) + *(4/9, *(27, phi0 ^ 2, *(xi, -2 + *(u3, xi ^ 3) + *(2, u2, xi ^ 2)) + *(S ^ 2, u ^ 7, -1 + *(u2, xi ^ 2)) + *(2, S, u3, (1 + *(u, xi)) ^ 2, -1 + *(u, xi))) + *(27, xi ^ 3, 4 + *(u, xi)) + *(81, u5, (S + *(S, u, xi)) ^ 2) + *(27, S ^ 3, u ^ 9, 1 + *(u, xi)) + *(81, S, u, (1 + *(u, xi)) ^ 3) + *(phi0 ^ 6, u3, (-1 + *(u, xi)) ^ 3, 1 + *(u, xi)) + *(9, u, phi0 ^ 4, (-1 + *(u, xi)) ^ 2, 1 + *(u, xi), 1 + *(S, u4) + *(u, xi)), expB1u4) + *(4/81, *(-1944, xi ^ 3) + *(27, phi0 ^ 2, *(30, xi) + *(u3, *(-8, U) + *(39, xi ^ 4)) + *(-2, u4, *(-9, xi ^ 5) + *(8, U, xi)) + *(18, u, xi ^ 2) + *(24, u2, xi ^ 3) + *(8, U, u ^ 7, xi ^ 4) + *(16, U, u6, xi ^ 3)) + *(162, u, U + *(-3, xi ^ 4)) + *(81, S ^ 4, u ^ 13, -6 + *(2, U, u4) + *(3, phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(162, U, u5, xi ^ 4) + *(648, U, xi, u2) + *(648, U, u4, xi ^ 3) + *(972, U, u3, xi ^ 2) + *(2, phi0 ^ 8, u5, (-1 + *(u, xi)) ^ 3, -15 + *(-1, U, u4) + *(15, u, xi) + *(36, u2, xi ^ 2) + *(U, xi, u5)) + *(3, phi0 ^ 10, u ^ 7, (-1 + *(u, xi)) ^ 4, -1 + *(2, u, xi)) + *(12, S, u, (3 + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi))) ^ 3, -6 + *(2, U, u4) + *(3, phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(54, S ^ 2, u5, (3 + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi))) ^ 2, -6 + *(2, U, u4) + *(3, phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(108, S ^ 3, u ^ 9, -6 + *(2, U, u4) + *(3, phi0 ^ 2, u2, -1 + *(2, u, xi)), 3 + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi))) + *(6, phi0 ^ 6, u3, (-1 + *(u, xi)) ^ 2, 1 + *(u, xi), -15 + *(-4, U, u4) + *(15, u, xi) + *(54, u2, xi ^ 2) + *(4, U, xi, u5)) + *(108, phi0 ^ 4, u3, (1 + *(u, xi)) ^ 2, -1 + *(u, xi), *(6, xi ^ 2) + *(-1, U, u2) + *(U, xi, u3)), expB1u4) + *(u, (1 + *(S, u4) + *(u, xi) + *(-1/3, phi0 ^ 2, u2) + *(1/3, xi, phi0 ^ 2, u3)) ^ 2, *(u2, *(Fxp ^ 2, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(Fyp ^ 2, coshGu4, expB2u4) + *(-4, B2c, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-4, Gc, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-2, B1s, coshGu4, expB2u4) + *(2, B1b, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(2, B2s, coshGu4, expB2u4) + *(2, B2b, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(2, Gs, expB2u4, sinhGu4) + *(2, Gb, exp(*(B2, u4) + *(2, B1, u4)), sinhGu4) + *(-12, Fx, xi_y, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-12, Fy, xi_x, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-2, B1p, xi_xx, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-2, B2p, xi_xx, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-2, B2p, xi_yy, coshGu4, expB2u4) + *(-2, Fxp, Fyp, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-2, Gp, xi_xx, exp(*(B2, u4) + *(2, B1, u4)), sinhGu4) + *(-2, Gp, xi_yy, expB2u4, sinhGu4) + *(2, B1p, xi_yy, coshGu4, expB2u4) + *(4, B2p, xi_xy, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(4, Gp, xi_xy, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(12, Fx, xi_x, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(12, Fy, xi_y, coshGu4, expB2u4)) + *(-8, u3, *(*(*(B1, xi_yy) + *(Fy, Fyp) + *(-1, B2, xi_yy) + *(-2, B2h, xi_y) + *(2, B1h, xi_y), expB2u4) + *(*(Fx, Fxp) + *(-1, B1, xi_xx) + *(-1, B2, xi_xx) + *(-2, B1t, xi_x) + *(-2, B2t, xi_x), exp(*(B2, u4) + *(2, B1, u4))) + *(2, *(G, xi_xy) + *(Gh, xi_x) + *(Gt, xi_y), exp(*(B1, u4) + *(B2, u4))), coshGu4) + *(-1, *(*(Fx, Fyp) + *(Fxp, Fy) + *(-2, B2, xi_xy) + *(-2, B2h, xi_x) + *(-2, B2t, xi_y), exp(*(B1, u4) + *(B2, u4))) + *(G, xi_xx, exp(*(B2, u4) + *(2, B1, u4))) + *(G, xi_yy, expB2u4) + *(2, Gh, xi_y, expB2u4) + *(2, Gt, xi_x, exp(*(B2, u4) + *(2, B1, u4))), sinhGu4)) + *(-4, u ^ 7, *(B2p, *(Fx ^ 2, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(Fy ^ 2, coshGu4, expB2u4) + *(-2, Fx, Fy, exp(*(B1, u4) + *(B2, u4)), sinhGu4)) + *(B1p, *(Fx ^ 2, exp(*(B2, u4) + *(2, B1, u4))) + *(-1, Fy ^ 2, expB2u4), coshGu4) + *(Gp, Fx ^ 2, exp(*(B2, u4) + *(2, B1, u4)), sinhGu4) + *(Gp, Fy ^ 2, expB2u4, sinhGu4) + *(-8, B2, B2h, xi_y, coshGu4, expB2u4) + *(-8, B2, B2t, xi_x, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-4, B1, B1h, xi_y, coshGu4, expB2u4) + *(-4, B1, B1t, xi_x, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-4, B1, Gt, xi_x, exp(*(B2, u4) + *(2, B1, u4)), sinhGu4) + *(-4, B1t, G, xi_x, exp(*(B2, u4) + *(2, B1, u4)), sinhGu4) + *(-4, G, Gh, xi_y, coshGu4, expB2u4) + *(-4, G, Gt, xi_x, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-2, B1, B2t, xi_x, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-2, B1, Fy, Fyp, coshGu4, expB2u4) + *(-2, B1, Gt, xi_y, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-2, B1h, G, xi_x, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-2, B1t, B2, xi_x, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-2, B2, Fx, Fyp, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-2, B2, Fxp, Fy, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-2, B2, Gh, xi_y, expB2u4, sinhGu4) + *(-2, B2, Gt, xi_x, exp(*(B2, u4) + *(2, B1, u4)), sinhGu4) + *(-2, B2h, G, xi_y, expB2u4, sinhGu4) + *(-2, B2t, G, xi_x, exp(*(B2, u4) + *(2, B1, u4)), sinhGu4) + *(-2, Fx, Fy, Gp, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-2, Fx, Fyp, G, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-2, Fxp, Fy, G, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(2, B1, B2h, xi_y, coshGu4, expB2u4) + *(2, B1, Fx, Fxp, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(2, B1, Gh, xi_x, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(2, B1h, B2, xi_y, coshGu4, expB2u4) + *(2, B1t, G, xi_y, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(2, B2, Fx, Fxp, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(2, B2, Fy, Fyp, coshGu4, expB2u4) + *(2, B2, Gh, xi_x, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(2, B2, Gt, xi_y, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(2, B2h, G, xi_x, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(2, B2t, G, xi_y, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(2, Fx, Fxp, G, exp(*(B2, u4) + *(2, B1, u4)), sinhGu4) + *(2, Fy, Fyp, G, expB2u4, sinhGu4) + *(4, B1, Gh, xi_y, expB2u4, sinhGu4) + *(4, B1h, G, xi_y, expB2u4, sinhGu4) + *(4, G, Gh, xi_x, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(4, G, Gt, xi_y, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(4, xi, Fx ^ 2, phi0 ^ 2, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(4, xi, Fy ^ 2, phi0 ^ 2, coshGu4, expB2u4) + *(8, B2, B2h, xi_x, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(8, B2, B2t, xi_y, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-8, Fx, Fy, xi, phi0 ^ 2, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-6, Fx, phi, phit, phi0 ^ 6, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-6, Fy, phi, phih, phi0 ^ 6, coshGu4, expB2u4) + *(6, Fx, phi, phih, phi0 ^ 6, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(6, Fy, phi, phit, phi0 ^ 6, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-24, Fx, phi, xi, xi_y, phi0 ^ 4, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-24, Fy, phi, xi, xi_x, phi0 ^ 4, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(24, Fx, phi, xi, xi_x, phi0 ^ 4, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(24, Fy, phi, xi, xi_y, phi0 ^ 4, coshGu4, expB2u4)) + *(2, u4, *(*(*(Fxh, Gp) + *(Fxp, Gh) + *(Fyp, Gt) + *(Fyt, Gp) + *(-40, G, xi_x, xi_y), exp(*(B1, u4) + *(B2, u4))) + *(*(8, Fx ^ 2) + *(-1, B1p, Fxt) + *(-1, B1t, Fxp) + *(-1, B2p, Fxt) + *(-1, B2t, Fxp) + *(2, phi0 ^ 6, phit ^ 2) + *(20, B1, xi_x ^ 2) + *(20, B2, xi_x ^ 2) + *(8, phi0 ^ 2, xi ^ 2, xi_x ^ 2) + *(-8, phit, xi, xi_x, phi0 ^ 4), exp(*(B2, u4) + *(2, B1, u4))) + *(*(8, Fy ^ 2) + *(B1h, Fyp) + *(B1p, Fyh) + *(-1, B2h, Fyp) + *(-1, B2p, Fyh) + *(-20, B1, xi_y ^ 2) + *(2, phi0 ^ 6, phih ^ 2) + *(20, B2, xi_y ^ 2) + *(8, phi0 ^ 2, xi ^ 2, xi_y ^ 2) + *(-8, phih, xi, xi_y, phi0 ^ 4), expB2u4), coshGu4) + *(*(B2h, Fxp, exp(*(B1, u4) + *(B2, u4))) + *(B2p, Fxh + Fyt, exp(*(B1, u4) + *(B2, u4))) + *(B2t, Fyp, exp(*(B1, u4) + *(B2, u4))) + *(-1, Fxp, Gt, exp(*(B2, u4) + *(2, B1, u4))) + *(-1, Fxt, Gp, exp(*(B2, u4) + *(2, B1, u4))) + *(-1, Fyh, Gp, expB2u4) + *(-1, Fyp, Gh, expB2u4) + *(-16, Fx, Fy, exp(*(B1, u4) + *(B2, u4))) + *(20, G, xi_x ^ 2, exp(*(B2, u4) + *(2, B1, u4))) + *(20, G, xi_y ^ 2, expB2u4) + *(-40, B2, xi_x, xi_y, exp(*(B1, u4) + *(B2, u4))) + *(-4, phih, phit, phi0 ^ 6, exp(*(B1, u4) + *(B2, u4))) + *(-16, xi_x, xi_y, phi0 ^ 2, xi ^ 2, exp(*(B1, u4) + *(B2, u4))) + *(8, phih, xi, xi_x, phi0 ^ 4, exp(*(B1, u4) + *(B2, u4))) + *(8, phit, xi, xi_y, phi0 ^ 4, exp(*(B1, u4) + *(B2, u4))), sinhGu4)) + *(2, u6, *(B1t, *(B2t, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-1, Gh, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(2, Gt, exp(*(B2, u4) + *(2, B1, u4)), sinhGu4)) + *(B2t, *(Gt, exp(*(B2, u4) + *(2, B1, u4)), sinhGu4) + *(-1, Gh, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-4, B2h, exp(*(B1, u4) + *(B2, u4)), sinhGu4)) + *(B1h ^ 2, coshGu4, expB2u4) + *(B1t ^ 2, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(Gh ^ 2, coshGu4, expB2u4) + *(Gt ^ 2, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(2, B2h ^ 2, coshGu4, expB2u4) + *(2, B2t ^ 2, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(B1h, Gt, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(B2h, Gh, expB2u4, sinhGu4) + *(-1, B1h, B2h, coshGu4, expB2u4) + *(-1, B2h, Gt, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-2, B1h, Gh, expB2u4, sinhGu4) + *(-2, Gh, Gt, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(2, Fx ^ 2, phi0 ^ 2, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(2, Fy ^ 2, phi0 ^ 2, coshGu4, expB2u4) + *(-56, B1, Fy, xi_y, coshGu4, expB2u4) + *(-56, B2, Fx, xi_y, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-56, B2, Fy, xi_x, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-56, Fx, G, xi_y, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-56, Fy, G, xi_x, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-4, Fx, Fy, phi0 ^ 2, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(18, phi ^ 2, phi0 ^ 6, xi_x ^ 2, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(18, phi ^ 2, phi0 ^ 6, xi_y ^ 2, coshGu4, expB2u4) + *(56, B1, Fx, xi_x, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(56, B2, Fx, xi_x, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(56, B2, Fy, xi_y, coshGu4, expB2u4) + *(56, Fx, G, xi_x, exp(*(B2, u4) + *(2, B1, u4)), sinhGu4) + *(56, Fy, G, xi_y, expB2u4, sinhGu4) + *(-36, xi_x, xi_y, phi ^ 2, phi0 ^ 6, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-16, Fx, xi_y, phi0 ^ 2, xi ^ 2, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-16, Fy, xi_x, phi0 ^ 2, xi ^ 2, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-12, Fx, phi, xi_y, phi0 ^ 4, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-12, Fy, phi, xi_x, phi0 ^ 4, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-8, Fx, phit, xi, phi0 ^ 4, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-8, Fy, phih, xi, phi0 ^ 4, coshGu4, expB2u4) + *(8, Fx, phih, xi, phi0 ^ 4, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(8, Fy, phit, xi, phi0 ^ 4, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(12, Fx, phi, xi_x, phi0 ^ 4, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(12, Fy, phi, xi_y, phi0 ^ 4, coshGu4, expB2u4) + *(16, Fx, xi_x, phi0 ^ 2, xi ^ 2, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(16, Fy, xi_y, phi0 ^ 2, xi ^ 2, coshGu4, expB2u4)) + *(4, u, *(Fxt + *(-1, Fxp, xi_x), coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(Fyh + *(-1, Fyp, xi_y), coshGu4, expB2u4) + *(-1, Fxh + Fyt + *(-1, Fxp, xi_y) + *(-1, Fyp, xi_x), exp(*(B1, u4) + *(B2, u4)), sinhGu4)) + *(4, u5, *(-5, B1h, Fy, coshGu4, expB2u4) + *(-5, B2h, Fx, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-5, B2t, Fy, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-5, Fx, Gh, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-5, Fy, Gt, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-2, B1, Fyh, coshGu4, expB2u4) + *(-2, B2, Fxh, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-2, B2, Fyt, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-2, Fxh, G, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-2, Fyt, G, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(2, B1, Fxt, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(2, B2, Fxt, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(2, B2, Fyh, coshGu4, expB2u4) + *(2, Fxt, G, exp(*(B2, u4) + *(2, B1, u4)), sinhGu4) + *(2, Fyh, G, expB2u4, sinhGu4) + *(5, B1t, Fx, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(5, B2h, Fy, coshGu4, expB2u4) + *(5, B2t, Fx, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(5, Fx, Gt, exp(*(B2, u4) + *(2, B1, u4)), sinhGu4) + *(5, Fy, Gh, expB2u4, sinhGu4) + *(B1p, Fy, xi_y, coshGu4, expB2u4) + *(B2p, Fx, xi_y, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(B2p, Fy, xi_x, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(Fx, Gp, xi_y, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(Fy, Gp, xi_x, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-1, B1p, Fx, xi_x, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-1, B2p, Fx, xi_x, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-1, B2p, Fy, xi_y, coshGu4, expB2u4) + *(-1, Fx, Gp, xi_x, exp(*(B2, u4) + *(2, B1, u4)), sinhGu4) + *(-1, Fy, Gp, xi_y, expB2u4, sinhGu4) + *(-2, B1, Fxp, xi_x, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-2, B2, Fxp, xi_x, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-2, B2, Fyp, xi_y, coshGu4, expB2u4) + *(-2, Fx, phih, phi0 ^ 4, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-2, Fxp, G, xi_x, exp(*(B2, u4) + *(2, B1, u4)), sinhGu4) + *(-2, Fy, phit, phi0 ^ 4, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-2, Fyp, G, xi_y, expB2u4, sinhGu4) + *(2, B1, Fyp, xi_y, coshGu4, expB2u4) + *(2, B2, Fxp, xi_y, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(2, B2, Fyp, xi_x, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(2, Fx, phit, phi0 ^ 4, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(2, Fxp, G, xi_y, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(2, Fy, phih, phi0 ^ 4, coshGu4, expB2u4) + *(2, Fyp, G, xi_x, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-12, phi, xi, phi0 ^ 4, xi_x ^ 2, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-12, phi, xi, phi0 ^ 4, xi_y ^ 2, coshGu4, expB2u4) + *(-6, phi, phih, xi_x, phi0 ^ 6, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-6, phi, phit, xi_y, phi0 ^ 6, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-4, Fx, xi, xi_x, phi0 ^ 2, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-4, Fy, xi, xi_y, phi0 ^ 2, coshGu4, expB2u4) + *(4, Fx, xi, xi_y, phi0 ^ 2, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(4, Fy, xi, xi_x, phi0 ^ 2, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(6, phi, phih, xi_y, phi0 ^ 6, coshGu4, expB2u4) + *(6, phi, phit, xi_x, phi0 ^ 6, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(24, phi, xi, xi_x, xi_y, phi0 ^ 4, exp(*(B1, u4) + *(B2, u4)), sinhGu4)) + *(4, u ^ 10, *(*(Fx, *(xi_x, *(16, B1 ^ 2) + *(16, G ^ 2) + *(32, B2 ^ 2) + *(16, B1, B2)) + *(9, Fx, phi ^ 2, phi0 ^ 6), exp(*(B2, u4) + *(2, B1, u4))) + *(Fy, *(xi_y, *(16, B1 ^ 2) + *(16, G ^ 2) + *(32, B2 ^ 2) + *(-16, B1, B2)) + *(9, Fy, phi ^ 2, phi0 ^ 6), expB2u4) + *(-16, B2, G, *(Fx, xi_y) + *(Fy, xi_x), exp(*(B1, u4) + *(B2, u4))), coshGu4) + *(-2, *(8, Fx, xi_y, G ^ 2 + *(2, B2 ^ 2)) + *(8, Fy, xi_x, G ^ 2 + *(2, B2 ^ 2)) + *(9, Fx, Fy, phi ^ 2, phi0 ^ 6), exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(16, G, *(B2, Fx, xi_x, exp(*(B2, u4) + *(2, B1, u4))) + *(B2, Fy, xi_y, expB2u4) + *(-2, B1, Fy, xi_y, expB2u4) + *(2, B1, Fx, xi_x, exp(*(B2, u4) + *(2, B1, u4))), sinhGu4)) + *(8, u8, *(B2, *(*(9, Fx ^ 2, exp(*(B2, u4) + *(2, B1, u4))) + *(9, Fy ^ 2, expB2u4) + *(-8, G, xi_x, xi_y, exp(*(B1, u4) + *(B2, u4))), coshGu4) + *(4, G, *(xi_x ^ 2, exp(*(B2, u4) + *(2, B1, u4))) + *(xi_y ^ 2, expB2u4), sinhGu4) + *(-18, Fx, Fy, exp(*(B1, u4) + *(B2, u4)), sinhGu4)) + *(8, B2 ^ 2, *(xi_x ^ 2, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(xi_y ^ 2, coshGu4, expB2u4) + *(-2, xi_x, xi_y, exp(*(B1, u4) + *(B2, u4)), sinhGu4)) + *(B1, *(-9, Fy ^ 2, expB2u4) + *(9, Fx ^ 2, exp(*(B2, u4) + *(2, B1, u4))) + *(-4, B2, xi_y ^ 2, expB2u4) + *(4, B2, xi_x ^ 2, exp(*(B2, u4) + *(2, B1, u4))), coshGu4) + *(4, B1 ^ 2, *(xi_x ^ 2, exp(*(B2, u4) + *(2, B1, u4))) + *(xi_y ^ 2, expB2u4), coshGu4) + *(4, G ^ 2, xi_x ^ 2, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(4, G ^ 2, xi_y ^ 2, coshGu4, expB2u4) + *(8, B1, G, *(xi_x ^ 2, exp(*(B2, u4) + *(2, B1, u4))) + *(-1, xi_y ^ 2, expB2u4), sinhGu4) + *(9, G, Fx ^ 2, exp(*(B2, u4) + *(2, B1, u4)), sinhGu4) + *(9, G, Fy ^ 2, expB2u4, sinhGu4) + *(-18, Fx, Fy, G, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-8, xi_x, xi_y, G ^ 2, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(2, Fx ^ 2, phi0 ^ 2, xi ^ 2, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(2, Fy ^ 2, phi0 ^ 2, xi ^ 2, coshGu4, expB2u4) + *(3, phi, Fx ^ 2, phi0 ^ 4, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(3, phi, Fy ^ 2, phi0 ^ 4, coshGu4, expB2u4) + *(-9, Fx, xi_y, phi ^ 2, phi0 ^ 6, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-9, Fy, xi_x, phi ^ 2, phi0 ^ 6, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-6, Fx, Fy, phi, phi0 ^ 4, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-4, Fx, Fy, phi0 ^ 2, xi ^ 2, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(9, Fx, xi_x, phi ^ 2, phi0 ^ 6, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(9, Fy, xi_y, phi ^ 2, phi0 ^ 6, coshGu4, expB2u4)) + *(8, u ^ 9, *(B1, *(B2t, Fx, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(Fy, Gt, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-1, B2h, Fy, coshGu4, expB2u4) + *(-1, Fx, Gh, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-2, Fy, Gh, expB2u4, sinhGu4) + *(2, B1h, Fy, coshGu4, expB2u4) + *(2, B1t, Fx, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(2, Fx, Gt, exp(*(B2, u4) + *(2, B1, u4)), sinhGu4)) + *(B1t, *(B2, Fx, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-1, Fy, G, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(2, Fx, G, exp(*(B2, u4) + *(2, B1, u4)), sinhGu4)) + *(B1h, Fx, G, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(B2, Fx, Gt, exp(*(B2, u4) + *(2, B1, u4)), sinhGu4) + *(B2, Fy, Gh, expB2u4, sinhGu4) + *(B2h, Fy, G, expB2u4, sinhGu4) + *(B2t, Fx, G, exp(*(B2, u4) + *(2, B1, u4)), sinhGu4) + *(-1, B1h, B2, Fy, coshGu4, expB2u4) + *(-1, B2, Fx, Gh, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-1, B2, Fy, Gt, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-1, B2h, Fx, G, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-1, B2t, Fy, G, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-4, B2, B2h, Fx, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-4, B2, B2t, Fy, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-2, B1h, Fy, G, expB2u4, sinhGu4) + *(-2, Fx, G, Gh, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-2, Fy, G, Gt, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(2, Fx, G, Gt, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(2, Fy, G, Gh, coshGu4, expB2u4) + *(4, B2, B2h, Fy, coshGu4, expB2u4) + *(4, B2, B2t, Fx, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-6, phi, xi, Fx ^ 2, phi0 ^ 4, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-6, phi, xi, Fy ^ 2, phi0 ^ 4, coshGu4, expB2u4) + *(12, Fx, Fy, phi, xi, phi0 ^ 4, exp(*(B1, u4) + *(B2, u4)), sinhGu4)) + *(32, u ^ 12, *(G ^ 2, *(Fx ^ 2, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(Fy ^ 2, coshGu4, expB2u4) + *(-2, Fx, Fy, exp(*(B1, u4) + *(B2, u4)), sinhGu4)) + *(2, B2 ^ 2, *(Fx ^ 2, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(Fy ^ 2, coshGu4, expB2u4) + *(-2, Fx, Fy, exp(*(B1, u4) + *(B2, u4)), sinhGu4)) + *(B1, *(B2, coshGu4) + *(2, G, sinhGu4), *(Fx ^ 2, exp(*(B2, u4) + *(2, B1, u4))) + *(-1, Fy ^ 2, expB2u4)) + *(B2, G, *(Fx ^ 2, exp(*(B2, u4) + *(2, B1, u4)), sinhGu4) + *(Fy ^ 2, expB2u4, sinhGu4) + *(-2, Fx, Fy, coshGu4, exp(*(B1, u4) + *(B2, u4)))) + *(B1 ^ 2, *(Fx ^ 2, exp(*(B2, u4) + *(2, B1, u4))) + *(Fy ^ 2, expB2u4), coshGu4)) + *(-2, Fxpt, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-2, Fyph, coshGu4, expB2u4) + *(2, Fxph, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(2, Fypt, exp(*(B1, u4) + *(B2, u4)), sinhGu4))

    nothing
end

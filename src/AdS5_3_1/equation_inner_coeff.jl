
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

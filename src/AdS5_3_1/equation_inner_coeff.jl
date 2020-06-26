
#= tilde, hat, etc, definitions

We use these macros as shorthand notation. For instance

  @tilde_inner("B1")
  @hat_inner("B2")

should expand to

  B1t = B1_x - (Fx * u * u + xi_x) * B1p
  B2h = B2_y - (Fy * u * u + xi_y) * B2p

etc.

=#
macro tilde_inner(fname::String)
    ft    = Symbol(fname, "t")
    f_x   = Symbol(fname, "_x")
    fp    = Symbol(fname, "p")
    return esc( :($ft = $f_x - (Fx * u * u + xi_x) * $fp) )
end
macro hat_inner(fname::String)
    fh    = Symbol(fname, "h")
    f_y   = Symbol(fname, "_y")
    fp    = Symbol(fname, "p")
    return esc( :($fh = $f_y - (Fy * u * u + xi_y) * $fp) )
end

macro bar_inner(fname::String)
    fb    = Symbol(fname, "b")
    f_xx  = Symbol(fname, "_xx")
    fpp   = Symbol(fname, "pp")
    fp_x  = Symbol(fname, "p_x")
    return esc( :($fb = $f_xx + (Fx * u * u + xi_x) * ( -2*($fp_x) + (Fx * u * u + xi_x) * ($fpp) )) )
end

macro star_inner(fname::String)
    fs    = Symbol(fname, "s")
    f_yy  = Symbol(fname, "_yy")
    fpp   = Symbol(fname, "pp")
    fp_y  = Symbol(fname, "p_y")
    return esc( :($fs = $f_yy + (Fy * u * u + xi_y) * ( -2*($fp_y) + (Fy * u * u + xi_y) * ($fpp) )) )
end

macro cross_inner(fname::String)
    fc    = Symbol(fname, "c")
    f_xy  = Symbol(fname, "_xy")
    fpp   = Symbol(fname, "pp")
    fp_x  = Symbol(fname, "p_x")
    fp_y  = Symbol(fname, "p_y")
    return esc( :($fc = $f_xy  - (Fx * u * u + xi_x) * ($fp_y) -
                  (Fy * u * u + xi_y) * ( $fp_x - (Fx * u * u + xi_x) * ($fpp) ) ) )
end



# assuming
# (A d_uu + B d_u + C Id) f = -S

function S_eq_coeff!(ABCS::Vector, vars::Tuple, ::Inner)
    (phi0, u, xi, B1, B1p, B2, B2p, G, Gp, phi, phip) = vars

    x0 = u^7
    x1 = u^6
    x2 = u^5
    x3 = Gp^2
    x4 = x1*x3
    x5 = 24*B2*B2p
    x6 = 8*x0
    x7 = G*Gp
    x8 = x6*x7
    x9 = u^8
    x10 = 48*B2^2
    x11 = x10*x9
    x12 = 3*B2p^2
    x13 = 16*G^2
    x14 = x13*x9
    x15 = u^4
    x16 = cosh(G*x15)^2
    x17 = B1p^2*x16
    x18 = B1*B1p*x16
    x19 = 16*B1^2*x16
    x20 = x19*x9
    x21 = phi0^2
    x22 = u^2
    x23 = u*xi
    x24 = 2*x23
    x25 = 3*phi
    x26 = phip - u*x25
    x27 = x26^2
    x28 = 4*phi0^6
    x29 = u^3
    x30 = phi0^4*(x24 - 1)
    x31 = x23 - 1
    x32 = x23 + 1
    x33 = phi*x32
    x34 = x22*xi^2
    x35 = -3*x23 + 2*x34 + 1
    x36 = 3*x32
    x37 = 9*x33
    x38 = 8*x7
    x39 = 8*x18
    x40 = x1*x31
    x41 = x2*x31

    ABCS[1] = 6*x0

    ABCS[2] = 48*x1

    ABCS[3] = x2*(-x0*x5 + x1*x12 + x1*x17 + x11 + x14 + x15*x27*x28 - x18*x6 + x20 + 4*x21*x22*(1 - x24)^2 + 8*x26*x29*x30 + x4 - x8 + 72)

    ABCS[4] = x15*(4*phi0^8*x27*x29*x31 + u*x28*(phip^2*x36 + 2*phip*u*(x35 - x37) + x22*x25*(6*x23 - 4*x34 + x37 - 2)) + x21*(-x0*x10 - x0*x13 - x0*x19 + x1*x38 + x11*xi + x12*x41 + x14*xi + x17*x41 - x2*x3 + x20*xi + 48*x22*xi^3 - x39*x40 + x4*xi - x40*x5 - x8*xi) + x29*x36*(-u*x38 - u*x39 - u*x5 + x10*x22 + x12 + x13*x22 + x17 + x19*x22 + x3) + 4*x30*(6*phip*x32 + u*(-18*x33 + x35)))/3

    nothing
end


# this is a coupled equation for Fx and Fy. the notation used is
#
# ( A11 d_uu Fx + A12 d_uu Fy + B11 d_u Fx + B12 d_u Fy + C11 Fx + C12 Fy ) = -S1
# ( A21 d_uu Fx + A22 d_uu Fy + B21 d_u Fx + B22 d_u Fy + C21 Fx + C22 Fy ) = -S2

function Fxy_eq_coeff!(AA::Matrix, BB::Matrix, CC::Matrix, SS::Vector, vars::Tuple, ::Inner)
    (
        phi0, u, xi, xi_x, xi_y,
        B1     ,    B2     ,    G      ,    phi    ,    S      ,
        B1p    ,    B2p    ,    Gp     ,    phip   ,    Sp     ,
        B1pp   ,    B2pp   ,    Gpp    ,                Spp    ,
        B1_x   ,    B2_x   ,    G_x    ,    phi_x  ,    S_x    ,
        B1_y   ,    B2_y   ,    G_y    ,    phi_y  ,    S_y    ,
        B1p_x  ,    B2p_x  ,    Gp_x   ,                Sp_x   ,
        B1p_y  ,    B2p_y  ,    Gp_y   ,                Sp_y   ,
    ) = vars

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


function Sd_eq_coeff!(ABCS::Vector, vars::Tuple, ::Inner)
    (
        potential, phi0, u, xi, xi_x, xi_y, xi_xx, xi_yy, xi_xy,
        B1     ,    B2     ,    G      ,    phi    ,    S      ,    Fx     ,    Fy     ,
        B1p    ,    B2p    ,    Gp     ,    phip   ,    Sp     ,    Fxp    ,    Fyp    ,
        B1pp   ,    B2pp   ,    Gpp    ,    phipp  ,    Spp    ,    Fxpp   ,    Fypp   ,
        B1_x   ,    B2_x   ,    G_x    ,    phi_x  ,    S_x    ,    Fx_x   ,    Fy_x   ,
        B1_y   ,    B2_y   ,    G_y    ,    phi_y  ,    S_y    ,    Fx_y   ,    Fy_y   ,
        B1p_x  ,    B2p_x  ,    Gp_x   ,    phip_x ,    Sp_x   ,    Fxp_x  ,    Fyp_x  ,
        B1p_y  ,    B2p_y  ,    Gp_y   ,    phip_y ,    Sp_y   ,    Fxp_y  ,    Fyp_y  ,
        B1_xx  ,    B2_xx  ,    G_xx   ,    phi_xx ,    S_xx   ,
        B1_yy  ,    B2_yy  ,    G_yy   ,    phi_yy ,    S_yy   ,
                    B2_xy  ,    G_xy   ,                S_xy
    ) = vars

    @tilde_inner("B1")
    @tilde_inner("B2")
    @tilde_inner("G")
    @tilde_inner("phi")
    @tilde_inner("S")
    @tilde_inner("Fx")
    @tilde_inner("Fy")

    @hat_inner("B1")
    @hat_inner("B2")
    @hat_inner("G")
    @hat_inner("phi")
    @hat_inner("S")
    @hat_inner("Fx")
    @hat_inner("Fy")

    @bar_inner("B1")
    @bar_inner("B2")
    @bar_inner("G")
    @bar_inner("phi")
    @bar_inner("S")

    @star_inner("B1")
    @star_inner("B2")
    @star_inner("G")
    @star_inner("phi")
    @star_inner("S")

    @tilde_inner("Sp")
    @tilde_inner("Fxp")
    @tilde_inner("Fyp")

    @hat_inner("Sp")
    @hat_inner("Fxp")
    @hat_inner("Fyp")

    @cross_inner("B2")
    @cross_inner("G")
    @cross_inner("S")

    phiouter = u * phi0 - u^2 * phi0 * xi + u^3 * phi0^3 * phi

    x0 = u^2
    x1 = u^4
    x2 = B1*x1
    x3 = exp(x2)
    x4 = S*x1
    x5 = 3*x4
    x6 = u*xi
    x7 = 3*x6
    x8 = x6 - 1
    x9 = phi0^2
    x10 = x0*x9
    x11 = x10*x8 + x7 + 3
    x12 = x11 + x5
    x13 = u^3
    x14 = 12*x13
    x15 = 3*Sp
    x16 = u*x9
    x17 = x13*x9
    x18 = x6 + 1
    x19 = x18 + x4
    x20 = 3*u
    x21 = xi^3
    x22 = x18^3
    x23 = S*u
    x24 = u^5
    x25 = u^9
    x26 = S^3*x25
    x27 = phi0^6
    x28 = x13*x27
    x29 = 4*x28
    x30 = phi0^4
    x31 = x8^2
    x32 = S^2
    x33 = u^7
    x34 = xi^2
    x35 = x0*x34
    x36 = x18^2
    x37 = 2*x13
    x38 = x13*x21
    x39 = 3*Sh
    x40 = 9*S
    x41 = u*x40
    x42 = x9*xi
    x43 = 2*xi_y
    x44 = u*(x39 + x41*xi_y + x42*x43)
    x45 = 9*x4
    x46 = 2*x6
    x47 = x46 - 1
    x48 = x10*x47 + x45 - 3
    x49 = Fy*x48
    x50 = x44 + x49
    x51 = B2*x1
    x52 = exp(x51)
    x53 = G*x1
    x54 = cosh(x53)
    x55 = x52*x54
    x56 = 3*St
    x57 = 2*xi_x
    x58 = u*(x41*xi_x + x42*x57 + x56)
    x59 = Fx*x48
    x60 = x58 + x59
    x61 = exp(2*x2 + x51)
    x62 = x54*x61
    x63 = exp(x2 + x51)
    x64 = sinh(x53)
    x65 = x63*x64
    x66 = 1296*x21
    x67 = xi^4
    x68 = u*x67
    x69 = phi0^8
    x70 = x24*x69
    x71 = x16*x34
    x72 = x1*x27
    x73 = 48*xi
    x74 = u^6
    x75 = 16*x74
    x76 = x69*xi
    x77 = x17*x67
    x78 = xi^5
    x79 = x1*x9
    x80 = x78*x79
    x81 = x13*x30
    x82 = x34*x81
    x83 = x1*x30
    x84 = 180*x67
    x85 = x24*x30
    x86 = 144*x74
    x87 = x24*x27
    x88 = x34*x87
    x89 = x27*x74
    x90 = x21*x89
    x91 = 12*x33
    x92 = x27*x67
    x93 = u^8
    x94 = x27*x78*x93
    x95 = x33*x34
    x96 = x69*x95
    x97 = 8*x93
    x98 = x21*x69
    x99 = x10 - 3*x36
    x100 = x12^2
    x101 = 12*Sp
    x102 = -12*x6
    x103 = 3*x10
    x104 = x24*x32
    x105 = x1*x67
    x106 = -x46
    x107 = UU(phiouter,potential)*x30
    x108 = phi0^10
    x109 = 864*UU(phiouter,potential)
    x110 = 432*UU(phiouter,potential)
    x111 = UU(phiouter,potential)*x108
    x112 = 96*x111
    x113 = UU(phiouter,potential)*x25
    x114 = phi0^12
    x115 = 8*x114
    x116 = x113*x115
    x117 = x72*xi
    x118 = x74*x76
    x119 = 72*x108
    x120 = x93*xi
    x121 = 1728*UU(phiouter,potential)
    x122 = 576*x111
    x123 = u^10
    x124 = UU(phiouter,potential)*x114
    x125 = 64*x124
    x126 = xi^6
    x127 = x24*x9
    x128 = x67*x85
    x129 = x30*x33
    x130 = x126*x129
    x131 = 648*x25
    x132 = x126*x27
    x133 = x93*x98
    x134 = x25*x67*x69
    x135 = x123*x69*x78
    x136 = u^11
    x137 = x126*x136*x69
    x138 = x108*x25*x34
    x139 = x123*x21
    x140 = u^12
    x141 = x140*x78
    x142 = u^13
    x143 = x126*x142
    x144 = 2592*UU(phiouter,potential)
    x145 = xi^8
    x146 = 5184*UU(phiouter,potential)
    x147 = xi^7
    x148 = x121*x147
    x149 = 1344*x111
    x150 = u^14
    x151 = u^15
    x152 = 224*x124
    x153 = 448*x124
    x154 = x12^4
    x155 = 2*UU(phiouter,potential)*x8^4*x83 - x103*x31 - 6
    x156 = 4*UU(phiouter,potential)*x10*x31
    x157 = 8*u
    x158 = phi^2
    x159 = x158*x27
    x160 = Fy*x0
    x161 = 4*G
    x162 = u*x161
    x163 = Gh + x162*(x160 + xi_y)
    x164 = Fx*x0
    x165 = 4*x1
    x166 = x42*xi_x
    x167 = x10*(4*x6 - 1) + 18*x4 - 3
    x168 = 2*x0
    x169 = -4*x13*x15 + 4*x48
    x170 = 36*x23
    x171 = 72*u
    x172 = x171*xi_y
    x173 = 12*x1
    x174 = Sh*x173
    x175 = 8*x42
    x176 = xi_y^2
    x177 = 8*x9
    x178 = B1*Sh
    x179 = 48*xi_y
    x180 = x24*xi_y
    x181 = 36*S
    x182 = x180*x181
    x183 = B2*x24
    x184 = x179*x183
    x185 = Fyp*x13
    x186 = x0*x176
    x187 = 144*S
    x188 = B1*x176
    x189 = S*x86
    x190 = B2*x189
    x191 = x6*x9
    x192 = 24*x191
    x193 = B1h*xi_y
    x194 = x1*x42
    x195 = 8*x194
    x196 = B2h*xi_y
    x197 = x10*xi
    x198 = 32*x42
    x199 = x198*x24
    x200 = x183*x198
    x201 = Fy^2
    x202 = 24*x51
    x203 = 117*x4
    x204 = Sp*x14
    x205 = B2*x93
    x206 = S*x205
    x207 = 72*x206
    x208 = 24*x2*(x5 - 1)
    x209 = 22*x6
    x210 = 8*x2
    x211 = x210*x47
    x212 = 8*x51
    x213 = 24*xi_y
    x214 = 78*x13
    x215 = Sh*x214
    x216 = 48*x33
    x217 = B2*x216
    x218 = Sh*x217
    x219 = 378*x4
    x220 = Sp*x13
    x221 = B1*x93
    x222 = S*x221
    x223 = 288*xi_y
    x224 = 8*x10
    x225 = B1*x75
    x226 = x9*xi_y
    x227 = B2*x75
    x228 = 68*x17
    x229 = xi*xi_y
    x230 = B1*x33
    x231 = 64*x42
    x232 = x231*xi_y
    x233 = B2*x33
    x234 = 4*x13
    x235 = x234*x48
    x236 = x24*x40
    x237 = 2*x194
    x238 = x17 + x20 - x236 - x237
    x239 = x171*xi_x
    x240 = St*x173
    x241 = x0*x56
    x242 = xi_x^2
    x243 = 48*xi_x
    x244 = B1*St
    x245 = B1t*xi_x
    x246 = x24*x245
    x247 = x183*x243
    x248 = B2t*xi_x
    x249 = x181*x24
    x250 = Fxp*x13
    x251 = x40*xi_x
    x252 = x0*x242
    x253 = B1*x242
    x254 = Fx^2
    x255 = 24*xi_x
    x256 = St*x214
    x257 = St*x217
    x258 = 288*xi_x
    x259 = x9*xi_x
    x260 = xi*xi_x
    x261 = 64*x166
    x262 = Fx*x213
    x263 = Fy*x255
    x264 = Fx*Fyp
    x265 = B2h*Fx
    x266 = B2t*Fy
    x267 = xi_x*xi_y
    x268 = 16*x9
    x269 = Fx*Fy
    x270 = B2*x74
    x271 = Fx*xi_y
    x272 = 48*x51
    x273 = Fy*xi_x
    x274 = x181*x33
    x275 = B2h*xi_x
    x276 = B2t*xi_y
    x277 = x269*x74
    x278 = x24*x269
    x279 = x0*x267
    x280 = 288*S
    x281 = B2*x269
    x282 = 288*x206
    x283 = 4*x24
    x284 = B2h*x283
    x285 = x283*x9
    x286 = 32*x205
    x287 = x227*x9
    x288 = x175*x74
    x289 = x228*xi
    x290 = Fyp*x57
    x291 = x231*x233
    x292 = Fxp*u
    x293 = Gh*xi_y
    x294 = Gh*x1
    x295 = 2*B1h
    x296 = 10*x13
    x297 = Fy*Gh
    x298 = x13*x161
    x299 = Gp*x0
    x300 = 8*B1
    x301 = x297*x33
    x302 = Fy*x33
    x303 = 8*G
    x304 = B1h*x303
    x305 = 4*B2
    x306 = B2*x283
    x307 = B2h*x161
    x308 = Fy*xi_y
    x309 = 56*x53
    x310 = Gp*x13
    x311 = x201*x74
    x312 = 36*G
    x313 = 2*x24
    x314 = Gp*x313
    x315 = 20*G
    x316 = G*x308
    x317 = 64*x221
    x318 = G*x123
    x319 = B1*x201
    x320 = 32*x319
    x321 = G*x74
    x322 = 32*x321
    x323 = B2*x318
    x324 = 16*x323
    x325 = G*x227
    x326 = x0*x163
    x327 = x157*xi_x
    x328 = Gt*x1
    x329 = 2*B1t
    x330 = Fx*x296
    x331 = Gt*x0
    x332 = Fx*Gt
    x333 = 8*xi_x
    x334 = B1*x24
    x335 = B1t*Fx
    x336 = x305*x33
    x337 = Gt*xi_x
    x338 = B2t*Fx
    x339 = x161*x33
    x340 = x161*x24
    x341 = Fx*xi_x
    x342 = x161*xi_x
    x343 = x254*x74
    x344 = G*x341
    x345 = B1*x254
    x346 = 2*u
    x347 = B2p*x346
    x348 = 6*u
    x349 = 8*x0
    x350 = B2*x349
    x351 = 10*x1
    x352 = 16*x13
    x353 = 4*Fyp
    x354 = 4*Fxp
    x355 = 16*x205
    x356 = 72*x33
    x357 = B2*x356
    x358 = Fx*x270
    x359 = 56*x271
    x360 = Fy*x270
    x361 = 4*x51
    x362 = 56*x273
    x363 = 40*x13
    x364 = B2*x363
    x365 = B2p*x1
    x366 = Fx*x43
    x367 = Fy*x57
    x368 = G*x97
    x369 = Fx*Gh
    x370 = Fy*Gt
    x371 = Gt*xi_y
    x372 = B2^2
    x373 = 64*x372
    x374 = x136*x269
    x375 = x25*x373
    x376 = x267*x33
    x377 = G^2
    x378 = 32*x377
    x379 = x25*x378
    x380 = 4*x83
    x381 = Fx*phih
    x382 = Fy*phit
    x383 = 24*x129
    x384 = phi*x383
    x385 = x42*x75
    x386 = 12*phi
    x387 = x386*x85
    x388 = x386*x89
    x389 = x85*xi
    x390 = 8*x389
    x391 = x386*x72
    x392 = x30*x73
    x393 = x392*x93
    x394 = phi*x393
    x395 = phi*x74
    x396 = x392*x395
    x397 = phi*x73*x83
    x398 = 36*x158
    x399 = x27*x398
    x400 = x25*x399
    x401 = x268*x95
    x402 = x33*x399
    x403 = 16*x34
    x404 = x127*x403
    x405 = x398*x87
    x406 = x17*x403
    x407 = B2p*x13
    x408 = x361 - x407 + 2
    x409 = 4*x230
    x410 = B1*x283
    x411 = Gh*xi_x
    x412 = B1t*x161
    x413 = G*x286
    x414 = 32*x270
    x415 = B1p*x346
    x416 = 12*u
    x417 = B1*x349
    x418 = 20*x1
    x419 = Fy*x418
    x420 = 16*x0
    x421 = x420*xi_y
    x422 = 8*Fyp
    x423 = 16*x93
    x424 = B1*x423
    x425 = B1h*Fy
    x426 = B1*x97
    x427 = B2h*Fy
    x428 = x300*x74
    x429 = B1*x74
    x430 = 112*x308
    x431 = Fyp*xi_y
    x432 = B2*x97
    x433 = 8*x270
    x434 = x165*x308
    x435 = G*x423
    x436 = G*x75
    x437 = 4*x311
    x438 = x25*x308
    x439 = 64*x438
    x440 = B1*B2
    x441 = 32*B2
    x442 = B1^2
    x443 = 128*x372
    x444 = 8*x83
    x445 = Fy*phih
    x446 = 32*x442
    x447 = x136*x201
    x448 = x176*x33
    x449 = phi*x213
    x450 = 24*phi*x89
    x451 = 16*x194
    x452 = 16*x389
    x453 = x30*x352
    x454 = 96*x30*x395*xi
    x455 = x159*x356
    x456 = 32*x127*x34
    x457 = B1p*x13
    x458 = 8*Fxp
    x459 = 112*x341
    x460 = Fxp*xi_x
    x461 = x165*x341
    x462 = 4*x343
    x463 = x25*x341
    x464 = 64*x463
    x465 = x136*x254
    x466 = B1*x441
    x467 = x242*x33
    x468 = Fx*phit
    x469 = phi*x255
    x470 = phi*x254

    ABCS[1] = 0

    ABCS[2] = -4*x0*x12^3*x3/9

    ABCS[3] = -8*x3*(x17*x8 + x19*x20)^2*(S*x14 - x0*x15 + x16*(x7 - 2) + 3*xi)/9

    ABCS[4] = u*x100*(x54*(-2*Fyph*x52 + u*x52*(B1h^2*x313 - B1h*B2h*x313 - B1h*x419 - B1h*x421 + B1p*x434 + B1p*x437 - B1s*x346 - B2*x136*x320 + B2h^2*x283 - 2*B2h*x185 + B2h*x419 + B2h*x421 - B2p*x434 - B2p*x437 + B2s*x346 + Fy*x422*x429 + Fy*x449*x85 + 2*Fyh*(-4*x2 + x408 + x457) + Fyp^2*u + Gh^2*x313 + phih^2*x29 - phih*x229*x453 + phih*x449*x72 - x160*x422 + x176*x364 - x176*x397 + x176*x405 + x176*x406 + x183*x430 + x185*x295 - x188*x33*x441 - x188*x363 + x193*x225 - x193*x433 + x196*x414 - x196*x428 + x201*x285 + x201*x352 + x201*x357 + x201*x384 - x201*x385 - x201*x394 + x201*x400 + x201*x401 + x210*x431 - x212*x431 + x286*x427 + x293*x436 + x297*x435 + x308*x416 - x308*x451 - x308*x454 + x308*x455 + x308*x456 - x319*x356 - x334*x430 - x347*xi_yy + x350*xi_yy - x353*xi_y - x360*x422 + x373*x447 + x373*x448 + x377*x439 + x378*x447 + x378*x448 + x415*xi_yy - x417*xi_yy + x424*x425 - x425*x432 - x426*x427 + x438*x443 - x439*x440 + x439*x442 + x444*x445 + x445*x450 - x445*x452 + x446*x447 + x446*x448) - x168*x63*(-B1h*Fx*x339 - B1h*x328 - B1h*x340*xi_x + B1t*x294 + B2h*x328 + B2t*x294 + Fxh*x298 - Fxh*x299 - Fxp*x326 - Fyp*x331 + Fyt*x298 - Fyt*x299 + G*x157*xi_xy + G*x267*x414 + 72*G*x277 + 40*G*x279 + 2*Gc + Gh*x327 + Gh*x330 - Gp*x269*x283 - 2*Gp*xi_xy + x157*x371 + x180*x412 - x185*x342 - x264*x340 + x265*x339 + x266*x339 + 32*x269*x323 + x271*x413 + x273*x413 + x275*x340 + x276*x340 + x296*x370 + x302*x412 + x306*x371 + x306*x411 - x310*x366 - x310*x367 + x336*x369 + x336*x370 + x359*x53 + x362*x53 + x369*x409 - x370*x409 - x371*x410 + x410*x411) + x61*(-2*Fxpt + u*(B1b*x346 - B1p*x461 - B1p*x462 + B1t^2*x313 + B1t*B2t*x313 + B2b*x346 - B2p*x461 - B2p*x462 + B2t^2*x283 - 2*B2t*x250 - Fx*x429*x458 + Fx*x469*x85 + Fxp^2*u + Fxt*(x210 + x212 - 2*x407 - 2*x457 + 4) + Gt^2*x313 + phit^2*x29 - phit*x260*x453 + phit*x469*x72 - x164*x458 + x183*x459 - x210*x460 - x212*x460 + x225*x245 + x242*x364 - x242*x397 + x242*x405 + x242*x406 + x245*x420 + x245*x433 + x248*x414 + x248*x420 + x248*x428 - x250*x329 + x253*x363 + x254*x285 + x254*x352 + x254*x357 - x254*x385 + x254*x400 + x254*x401 + x286*x338 + x332*x435 + x334*x459 + x335*x418 + x335*x424 + x335*x432 + x337*x436 + x338*x418 + x338*x426 + x341*x416 - x341*x451 - x341*x454 + x341*x455 + x341*x456 + x345*x356 - x347*xi_xx + x350*xi_xx - x354*xi_x - x358*x458 + x373*x465 + x373*x467 + x377*x464 + x378*x465 + x378*x467 + x383*x470 - x393*x470 - x415*xi_xx + x417*xi_xx + x440*x464 + x442*x464 + x443*x463 + x444*x468 + x446*x465 + x446*x467 + x450*x468 - x452*x468 + x465*x466 + x466*x467))) + x64*(x168*(x52*(B2h*x294 - Fy*x310*x43 + Fyh*x298 - Fyh*x299 - Fyp*x326 - Gp*xi_yy + Gs + x157*x293 + x162*xi_yy + x176*x325 - x180*x304 + x180*x307 + x186*x315 - x188*x322 - x201*x314 + x201*x324 - x24*x293*x300 + x286*x316 + x293*x306 - x294*x295 + x296*x297 - x300*x301 + x301*x305 - x302*x304 + x302*x307 + x308*x309 + x311*x312 - x316*x317 - x318*x320) + x61*(B2t*x328 - Fx*Fxp*x340 - Fx*x310*x57 - Fxp*x331 + Fxt*x0*(-Gp + x162) + Gb - Gp*xi_xx + Gt*x327 + Gt*x330 + Gt*x333*x334 + x162*xi_xx + x242*x325 + x246*x303 + x248*x340 - x250*x342 + x252*x315 + x253*x322 - x254*x314 + x254*x324 + x286*x344 + x300*x33*x332 + x303*x33*x335 + x306*x337 + x309*x341 + x312*x343 + x317*x344 + 32*x318*x345 + x328*x329 + x332*x336 + x338*x339)) + x63*(2*Fxph + 2*Fypt - x346*(B2c*x346 - B2h*x250 - 4*B2p*x277 - B2t*x185 + B2t*x284 + Fxh*x408 - Fxp*x361*xi_y - Fxp*x43 + Fyp*x292 - Fyp*x361*xi_x + Fyt*x408 + Gh*Gt*x313 + Gh*x321*x333 + phih*phit*x29 - phih*x333*x81*xi + phih*x391*xi_x - 8*phit*x229*x81 + phit*x391*xi_y - x160*x354 - x164*x353 + x183*x359 + x183*x362 - x195*x271 - x195*x273 + x227*x275 + x227*x276 + x265*x351 + x265*x355 + x266*x351 + x266*x355 + x267*x364 - x267*x397 + x267*x405 + x267*x406 + x269*x285 + x269*x352 + x269*x357 + x269*x384 - x269*x385 - x269*x394 + x269*x400 + x269*x401 + x271*x348 + x271*x375 + x271*x379 + x271*x387 - x271*x396 + x271*x402 + x271*x404 + x273*x348 + x273*x375 + x273*x379 + x273*x387 - x273*x396 + x273*x402 + x273*x404 + x275*x349 + x276*x349 - x290 + 8*x321*x371 - x347*xi_xy + x350*xi_xy - x353*x358 - x354*x360 - x365*x366 - x365*x367 + x368*x369 + x368*x370 + x373*x374 + x373*x376 + x374*x378 + x376*x378 + x380*x381 + x380*x382 + x381*x388 - x381*x390 + x382*x388 - x382*x390))))/9 + 2*x0*x12*(-x165*x54*x63*(Fx*x162*(x13*x39 + x167*xi_y + x168*x49) + Fy*x162*(x13*x56 + x167*xi_x) + Gh*x58 + Gh*x59 + Gt*x44 + Gt*x49 + x0*x161*(x39*xi_x + xi_y*(4*x166 + 18*x23*xi_x + x56))) + x165*x64*(x163*x50*x52 + x60*x61*(Gt + x162*(x164 + xi_x))) + x55*(Fyh*x169 + u*(-B1h*x174 - B1h*x182 + B2h*x174 + B2h*x182 + Fy*(-B1h*x235 - B2h*x234*(-x45 + x9*(x0 - x37*xi) + 3) + Fyp*x238 - x178*x216 + x179*x2 - x179*x51 + x206*x223 - x213*x220 - x213 + x215 + x218 + x219*xi_y - x222*x223 - x224*xi_y + x225*x226 - x226*x227 + x228*x229 - x230*x232 + x232*x233) - Fyp*x0*x39 - Fyp*x197*x43 + Sh*x172 + Sh*x184 + 12*Ss - x101*xi_yy - x168*x201*(x10*(-x209 + x211 + x212*(x106 + 1) + 9) + x202 - x203 + x204 - x207 + x208 + 15) + x170*xi_yy + x175*xi_yy + x176*x177 + x176*x190 + x176*x192 + x176*x200 - x178*x179*x24 - x185*x40*xi_y + x186*x187 - x188*x189 - x188*x199 - x193*x195 + x195*x196)) + x62*(Fxt*x169 + u*(B1t*x240 + B2t*x240 + Fx*(B1t*x235 + B2t*x235 + Fxp*x238 - x2*x243 + x206*x258 + x216*x244 + x219*xi_x - x220*x255 + x222*x258 - x224*xi_x - x225*x259 - x227*x259 + x228*x260 + x230*x261 + x233*x261 - x243*x51 - x255 + x256 + x257) - Fxp*x197*x57 - Fxp*x241 + 12*Sb + St*x239 + St*x247 - x101*xi_xx + x168*x254*(x10*(x209 + x211 + x212*x47 - 9) - x202 + x203 - x204 + x207 + x208 - 15) + x170*xi_xx + x175*xi_xx + x177*x242 + x181*x246 + x187*x252 + x189*x253 + x190*x242 + x192*x242 + x195*x245 + x195*x248 + x199*x253 + x200*x242 + x24*x243*x244 + x248*x249 - x250*x251)) - x65*(Fxh*x169 + Fyt*x169 + u*(B2h*x240 + B2t*x174 - 60*Fx*x160 + Fx*x215 + Fx*x218 - Fx*x284*x9 + Fy*x256 + Fy*x257 - Fyp*x241 + 468*S*x277 + 24*Sc + Sh*x239 + Sh*x247 - 48*Sp*x278 - 24*Sp*xi_xy + St*x172 + St*x184 + x123*x280*x281 - x14*x265 - x14*x266 + x17*x264 + x183*x231*x267 - x185*x251 + x191*x243*xi_y + x195*x275 + x195*x276 - x197*x290 + x20*x264 + x219*x271 + x219*x273 - x220*x262 - x220*x263 - x224*x271 - x224*x273 + 72*x23*xi_xy + x231*x25*x281 - x236*x264 - x237*x264 + x249*x275 + x249*x276 - x262 - x263 + x265*x274 + x265*x288 + x266*x274 - x266*x285 + x266*x288 + x267*x268 + x267*x270*x280 - 96*x269*x270 - x269*x286*x9 - 36*x269*x79 - x271*x272 + x271*x282 - x271*x287 + x271*x289 + x271*x291 - x272*x273 + x273*x282 - x273*x287 + x273*x289 + x273*x291 + 88*x278*x42 + x279*x280 - x292*x50 + 16*x42*xi_xy)))/9 + 4*x24*(-x50^2*x55 + x50*x65*(2*x58 + 2*x59) - x60^2*x62)/9 + x3*(36*u*x18*x19*x30*x31 + 108*x18*x26 + x18*x29*x8^3 + 108*x21*(x6 + 4) + 324*x22*x23 + 324*x24*(S*x6 + S)^2 + 108*x9*(S*x36*x37*x8 + x32*x33*(x35 - 1) + xi*(2*x35 + x38 - 2)))/9 - x3*(432*x10*x21 + x100*x101*x99 - 36*x104*(-x103*(x102 + 9*x35 + 8*x38 - 12) - 9*x36*(6*x6 + 5) + x83*(8*x6 - 7)) - 216*x21*x83 - 12*x23*(-9*x10*(10*x105 + x35 + 20*x38 - 16*x6 - 7) - 27*x22*(x7 + 1) - 3*x83*(7*x105 + x106 - 22*x35 + 2*x38 + 11) + x89*(x102 + 7*x35 + 5)) - 324*x26*x99 - 24*x28 + x30*x78*x86 + 432*x42 - x66 - 324*x68 + 4*x70 + 324*x71 + x72*x73 - x75*x76 + 540*x77 + 216*x80 - 252*x82 + x84*x85 + 60*x88 - 96*x90 - x91*x92 + 24*x94 + 20*x96 - x97*x98)/27 + x3*(324*S^4*x142*x155 + UU(phiouter,potential)*u^17*x115*x145 + 3888*UU(phiouter,potential)*x128 - 4320*UU(phiouter,potential)*x134 - 1344*UU(phiouter,potential)*x138 + phi^4*x116*x154 - 32*phi^3*x111*x154*x33*x8 - phi*x157*x8*(x156 - 3)*(phi0^3*x0*x8 + 3*phi0*x19)^4 - u^16*x125*x147 + 648*u*x107 - x10*x66 + 216*x104*x11^2*x155 + x107*x131*x145 - x108*x136*x84 + 240*x108*x139 - 12*x108*x143 - x108*x91 + x109*x136*x145*x27 - x109*x28 + 48*x11^3*x155*x23 + 432*x11*x155*x26 + x110*x142*x145*x69 + x110*x70 + x112*x145*x151 - x112*x33 - 1728*x113*x132 + x116 + x117*x121 + 720*x117 - x118*x121 - 480*x118 + x119*x120 + x119*x141 + x120*x122 + x121*x133 + x121*x135 + x121*x137 + x121*x88 + x121*x96 - x122*x147*x150 - x123*x125*xi - x123*x148*x27 + 560*x124*x142*x67 - 972*x126*x127 + x126*x151*x152 + 2592*x128 - x130*x144 - 1296*x130 - x131*x132 + 96*x133 - 744*x134 + 576*x135 + x136*x152*x34 - 144*x137 - 180*x138 + x139*x149 + x14*x154*x159*(x156 - 1) - x140*x148*x69 - x140*x153*x21 - x141*x149 + x143*x149 - x144*x82 - x146*x90 + x146*x94 - x150*x153*x78 - 7776*x21 - 360*x28 + 360*x33*x92 + 3240*x42 - 1944*x68 + 120*x70 + 972*x71 - 1620*x77 - 1944*x80 - 1296*x82 + 648*x88 - 2016*x90 + 1296*x94 + 576*x96)/81

    nothing
end


function B2d_eq_coeff!(ABCS::Vector, vars::Tuple, ::Inner)
    (
        phi0, u, xi, xi_x, xi_y, xi_xx, xi_yy, xi_xy,
        B1     ,    B2     ,    G      ,    phi    ,    S      ,    Fx     ,    Fy     ,  Sd,
        B1p    ,    B2p    ,    Gp     ,    phip   ,    Sp     ,    Fxp    ,    Fyp    ,
        B1pp   ,    B2pp   ,    Gpp    ,    phipp  ,    Spp    ,    Fxpp   ,    Fypp   ,
        B1_x   ,    B2_x   ,    G_x    ,    phi_x  ,    S_x    ,    Fx_x   ,    Fy_x   ,
        B1_y   ,    B2_y   ,    G_y    ,    phi_y  ,    S_y    ,    Fx_y   ,    Fy_y   ,
        B1p_x  ,    B2p_x  ,    Gp_x   ,    phip_x ,    Sp_x   ,    Fxp_x  ,    Fyp_x  ,
        B1p_y  ,    B2p_y  ,    Gp_y   ,    phip_y ,    Sp_y   ,    Fxp_y  ,    Fyp_y  ,
        B1_xx  ,    B2_xx  ,    G_xx   ,    phi_xx ,    S_xx   ,
        B1_yy  ,    B2_yy  ,    G_yy   ,    phi_yy ,    S_yy   ,
                    B2_xy  ,    G_xy   ,                S_xy
    ) = vars

    @tilde_inner("B1")
    @tilde_inner("B2")
    @tilde_inner("G")
    @tilde_inner("phi")
    @tilde_inner("S")
    @tilde_inner("Fx")
    @tilde_inner("Fy")

    @hat_inner("B1")
    @hat_inner("B2")
    @hat_inner("G")
    @hat_inner("phi")
    @hat_inner("S")
    @hat_inner("Fx")
    @hat_inner("Fy")

    @bar_inner("B1")
    @bar_inner("B2")
    @bar_inner("G")
    @bar_inner("phi")
    @bar_inner("S")

    @star_inner("B1")
    @star_inner("B2")
    @star_inner("G")
    @star_inner("phi")
    @star_inner("S")

    @tilde_inner("Fxp")
    @tilde_inner("Fyp")

    @hat_inner("Fxp")
    @hat_inner("Fyp")

    @cross_inner("B2")
    @cross_inner("G")
    @cross_inner("S")

    x0 = u^2
    x1 = u^4
    x2 = B1*x1
    x3 = exp(x2)
    x4 = u*xi
    x5 = S*x1
    x6 = phi0^2
    x7 = x0*x6
    x8 = 3*x4 + 3*x5 + x7*(x4 - 1) + 3
    x9 = x3*x8^3
    x10 = 3*Sp
    x11 = u^3
    x12 = x10*x11
    x13 = -x12
    x14 = x13 + 3
    x15 = 6*x4
    x16 = 4*x4
    x17 = 4*u
    x18 = 6*x1
    x19 = u^5
    x20 = B2*x1
    x21 = exp(x20)
    x22 = G*x1
    x23 = cosh(x22)
    x24 = x21*x23
    x25 = 3*Sh
    x26 = 9*u*xi_y
    x27 = 2*x6*xi
    x28 = x27*xi_y
    x29 = u*(S*x26 + x25 + x28)
    x30 = 9*x5
    x31 = 2*x4 - 1
    x32 = x30 + x31*x7 - 3
    x33 = Fy*x32
    x34 = x29 + x33
    x35 = 2*x2
    x36 = exp(x20 + x35)
    x37 = x23*x36
    x38 = 3*St
    x39 = S*u*xi_x
    x40 = 2*xi
    x41 = x6*xi_x
    x42 = x40*x41
    x43 = x38 + 9*x39 + x42
    x44 = u*x43
    x45 = Fx*x32
    x46 = x44 + x45
    x47 = sinh(x22)
    x48 = exp(x2 + x20)
    x49 = x47*x48
    x50 = 2*u
    x51 = G*x17
    x52 = Fy*x0
    x53 = Gh + x51*(x52 + xi_y)
    x54 = Fx*x0
    x55 = 4*xi
    x56 = x41*x55
    x57 = x11*x38
    x58 = 18*x5
    x59 = x7*(x16 - 1)
    x60 = x58 + x59 - 3
    x61 = x11*x25
    x62 = 2*x0*x31*x6 + x58 - 6
    x63 = x13 + x32
    x64 = 9*S*u
    x65 = u*xi_y
    x66 = B1h*x1
    x67 = 6*Sh
    x68 = Fyp*x0
    x69 = xi_y^2
    x70 = 2*x6
    x71 = B1*x19*xi_y
    x72 = 9*S
    x73 = B1h*x19*xi_y
    x74 = 24*B2
    x75 = 18*S
    x76 = B2h*x19*xi_y
    x77 = Fyp*xi_y
    x78 = 18*S*x11
    x79 = 36*S*x0
    x80 = u^6
    x81 = 36*B1*S*x80
    x82 = 72*B2*S*x80
    x83 = 6*u*x6*xi
    x84 = x6*xi_y
    x85 = 4*Fyp
    x86 = x85*xi_y
    x87 = x0*x6*xi
    x88 = 8*B1*x19*x6*xi
    x89 = 16*B2*x19*x6*xi
    x90 = 2*x0
    x91 = Fy^2
    x92 = 12*B2
    x93 = x1*x92
    x94 = -x93
    x95 = -x30
    x96 = u^8
    x97 = 36*S*x96
    x98 = B2*x97
    x99 = B2*x6*x80
    x100 = 4*x99
    x101 = x6*xi
    x102 = x101*x11
    x103 = u^7
    x104 = B2*x103*x6*xi
    x105 = 8*x104
    x106 = x32*x35
    x107 = 4*x20
    x108 = x107 + x35 - 1
    x109 = x11*x32
    x110 = x6*(x0 - x11*x40) + x95 + 3
    x111 = B2h*x11
    x112 = 2*xi_y
    x113 = 6*B1
    x114 = x1*x113
    x115 = B1*x97
    x116 = 27*x5
    x117 = 72*S
    x118 = B2*x117*x96
    x119 = -x116 + x118 + x12 + x94 + 3
    x120 = u*xi_x
    x121 = B1t*x1
    x122 = 6*St
    x123 = Fxp*x0
    x124 = xi_x^2
    x125 = St*x19
    x126 = 9*S*xi_x
    x127 = B1t*x19
    x128 = Fxp*xi_x
    x129 = B2t*x1
    x130 = 4*Fxp
    x131 = x130*xi_x
    x132 = Fx^2
    x133 = 3*x11
    x134 = x103*x72
    x135 = x19*x6
    x136 = x27*x80
    x137 = 2*xi_x
    x138 = x110 + x12
    x139 = 3*Fx
    x140 = x139*xi_y
    x141 = Fy*xi_x
    x142 = 3*x141
    x143 = Fyp*u
    x144 = B2h*x1
    x145 = B2t*Fy
    x146 = 2*xi_xy
    x147 = 24*Fx*Fy
    x148 = B2*x80
    x149 = 12*B2*x103
    x150 = Fx*xi_y
    x151 = B2h*Fx
    x152 = B2h*x19
    x153 = B2t*xi_y
    x154 = 9*S*x19
    x155 = Fx*Fy*x80
    x156 = Fx*x19
    x157 = Sp*x11
    x158 = xi_x*xi_y
    x159 = Fx*x11*x6
    x160 = u^10
    x161 = B2*Fx*Fy*x160
    x162 = B2*x80*xi_x
    x163 = x6*xi_x*xi_y
    x164 = B2*Fx*x96
    x165 = 4*Fy*xi_x
    x166 = Fx*Fy
    x167 = Fx*x6*xi
    x168 = (2*Fyp)*xi_x
    x169 = u^9
    x170 = 16*Fy*x6*xi
    x171 = B2*Fx*x103
    x172 = 16*Fy*xi_x
    x173 = 8*u*xi_y
    x174 = Gp*x0
    x175 = 2*Gh
    x176 = 10*Fy*x11
    x177 = 4*G*x11
    x178 = 8*Gh
    x179 = B1*Fy*x103
    x180 = 8*Fy
    x181 = B1h*G*x103
    x182 = 8*G
    x183 = B2*Fy*x103
    x184 = 8*B2*Gh*x19
    x185 = G*x103
    x186 = Fy*xi_y
    x187 = 56*G*x1
    x188 = Gp*x11
    x189 = 36*G*x80
    x190 = 2*Gp*x19
    x191 = 20*G*x0
    x192 = 64*G*xi_y
    x193 = 32*B1*G*x160
    x194 = 32*B1*G*x80
    x195 = 32*B2*G*x160
    x196 = 32*B2*G*x80
    x197 = 8*u*xi_x
    x198 = Gt*x0
    x199 = 2*Gt
    x200 = 10*Fx*x11
    x201 = 8*Gt
    x202 = B1*Fx*x103
    x203 = B1*xi_x
    x204 = B1t*Fx
    x205 = 8*G*x103
    x206 = 8*G*xi_x
    x207 = 8*B2*Gt*x19
    x208 = B2t*Fx
    x209 = B2t*xi_x
    x210 = 8*G*x19
    x211 = Fx*G*x19
    x212 = Fx*xi_x
    x213 = Fx*Gp*x11
    x214 = G*x96
    x215 = 4*B2p*u
    x216 = Fx*u
    x217 = 16*B2*x0
    x218 = 2*x19
    x219 = 20*Fx
    x220 = 2*Fxp*x11
    x221 = 16*x0*xi_x
    x222 = 20*Fy
    x223 = 2*Fyp*x11
    x224 = 16*x0*xi_y
    x225 = 8*B2h
    x226 = 144*B2*x103
    x227 = B2*Fx
    x228 = 8*Fyp*x80
    x229 = 112*x19*xi_y
    x230 = B2*Fy
    x231 = 8*Fxp*x80
    x232 = 8*x20
    x233 = 112*x19
    x234 = 80*B2*x11
    x235 = 4*Fx
    x236 = B2p*x1*xi_y
    x237 = G*xi_x
    x238 = G*xi_y
    x239 = B2^2
    x240 = 32*x239
    x241 = u^11
    x242 = Fx*Fy*x241
    x243 = 32*x169*x239
    x244 = x103*xi_x*xi_y
    x245 = G^2
    x246 = 32*x245
    x247 = 32*x169*x245
    x248 = phi0^4
    x249 = phih*x1*x248
    x250 = phit*x248
    x251 = phi0^6
    x252 = 4*x11*x251
    x253 = 12*phi*x19*x248
    x254 = Fx*phi*x251*x80
    x255 = 8*x1*x6*xi
    x256 = 8*phih*x248*xi
    x257 = phit*x248*xi
    x258 = phi*x1*x251
    x259 = phi*x1*x251*xi_y
    x260 = 48*phi*x248*x96*xi
    x261 = 48*phi*x248*x80*xi
    x262 = 48*phi*x1*x248*xi
    x263 = phi^2
    x264 = 36*x169*x251*x263
    x265 = xi^2
    x266 = 36*x103*x251*x263
    x267 = x19*x265*x6
    x268 = 36*x19*x251*x263
    x269 = -x232
    x270 = B2p*x11
    x271 = 2*x270
    x272 = x269 + x271 + 2
    x273 = 4*B1h*x19
    x274 = 4*B1t*x19
    x275 = G*xi_x*xi_y
    x276 = 64*G
    x277 = 2*B1p*u
    x278 = 8*B1*x0
    x279 = B2h*xi_y
    x280 = 32*x0
    x281 = 16*x11
    x282 = B1*B1h
    x283 = 16*Fy*x96
    x284 = 16*x80*xi_y
    x285 = 16*B1*x80
    x286 = B1*Fy
    x287 = 8*x2
    x288 = 16*B2*Fy*x96
    x289 = 16*B2*x80
    x290 = 16*x20
    x291 = G*Gh
    x292 = 72*B1*x103
    x293 = 40*B1*x11
    x294 = 4*B1p*x80
    x295 = 8*B2p*x80
    x296 = 128*B2*x169
    x297 = 64*B1*B2*x241
    x298 = 64*B1*B2*x103
    x299 = B1^2
    x300 = 64*x169*x299
    x301 = 64*x169*x239
    x302 = 64*x169*x245
    x303 = 32*x241*x299
    x304 = 32*x103*x299
    x305 = 32*x239*x241
    x306 = 32*x103*x239
    x307 = 32*x241*x245
    x308 = 4*x19*x6
    x309 = 32*x103*x245
    x310 = 24*phi*x248
    x311 = 24*phih
    x312 = phih*x248*xi
    x313 = 24*phi*x103*x248
    x314 = 16*x6*x80*xi
    x315 = 96*phi*x248*x80*xi
    x316 = 72*x103*x251*x263
    x317 = 32*x19*x265*x6
    x318 = 16*x103*x265*x6
    x319 = 16*x11*x265*x6
    x320 = B1p*x11
    x321 = B1*B1t
    x322 = 16*Fx*x96
    x323 = 16*x80*xi_x
    x324 = B1*Fx*xi_x
    x325 = 16*B2*x96
    x326 = G*Gt
    x327 = 16*xi

    ABCS[1] = 0

    ABCS[2] = -4*x0*x3*x8^4/27

    ABCS[3] = -2*u*x9*(x14 + x15 + 15*x5 + x7*(x16 - 3))/9

    ABCS[4] = -u*x8^2*(x23*(-2*Fyph*x21 + u*x21*(B1*B2h*x283 + B1*Fy*x296*xi_y + B1h^2*x218 + B1h*B2*x284 + B1h*x223 - B1h*x224 + B1h*x288 + 4*B1p*Fy*x1*xi_y - B1s*x50 - B2h^2*x218 + B2h*x273 - B2h*x288 - B2s*x17 + 16*Fy*Fyp*x148 + Fy*phi*x251*x311*x80 - 40*Fy*x144 + Fy*x19*x310*xi_y - 16*Fy*x19*x312 + 12*Fy*x65 + 2*Fyh*(-4*x2 + x269 + x271 + x320 + 2) + Fyp^2*u - 8*Fyp*x52 + Gh^2*x218 + phih^2*x252 - x1*x170*xi_y - 16*x11*x312*xi_y + x111*x85 + x180*x236 + x180*x249 + x186*x300 - x186*x301 + x186*x302 - x186*x315 + x186*x316 + x186*x317 - 224*x19*x230*xi_y + x215*xi_yy - x217*xi_yy - x222*x66 - x226*x91 + x228*x286 - x229*x286 - x234*x69 + x259*x311 - x260*x91 - x262*x69 + x264*x91 + x268*x69 + x277*xi_yy - x278*xi_yy - x279*x280 + x279*x285 - x279*x289 + x281*x91 + x282*x283 + x282*x284 + x283*x291 + x284*x291 + x287*x77 + x290*x77 - x292*x91 - x293*x69 + x294*x91 + x295*x91 + x297*x91 + x298*x69 + x303*x91 + x304*x69 - x305*x91 - x306*x69 + x307*x91 + x308*x91 + x309*x69 + x313*x91 - x314*x91 + x318*x91 + x319*x69 - x86) + x36*(-2*Fxpt + u*(-B1*B2t*x322 - B1*Fx*x231 + B1b*x50 - 4*B1p*Fx*x1*xi_x + B1t^2*x218 - B1t*B2*x323 - B1t*x220 + B1t*x221 - 224*B2*Fx*x19*xi_x - B2b*x17 + 8*B2p*Fx*x1*xi_x - B2t^2*x218 + B2t*x11*x130 - B2t*x274 + Fx*Fxp*x289 - Fx*phit*x19*x248*x327 + 8*Fx*x1*x250 - Fx*x1*x327*x6*xi_x - 40*Fx*x129 + Fx*x19*x310*xi_x + Fxp^2*u - 8*Fxp*x54 + Fxt*(4*x270 + x287 - x290 - 2*x320 + 4) + Gt^2*x218 + phit^2*x252 + 24*phit*x254 + 24*phit*x258*xi_x - 16*x11*x257*xi_x + x121*x219 - x124*x234 - x124*x262 + x124*x268 + x124*x293 - x124*x298 + x124*x304 - x124*x306 + x124*x309 + x124*x319 - x128*x287 + x128*x290 - x131 - x132*x226 - x132*x260 + x132*x264 + x132*x281 + x132*x292 - x132*x294 + x132*x295 - x132*x297 + x132*x303 - x132*x305 + x132*x307 + x132*x308 + x132*x313 - x132*x314 + x132*x318 - x204*x325 - x208*x325 - x209*x280 - x209*x285 - x209*x289 + x212*x300 - x212*x301 + x212*x302 - x212*x315 + x212*x316 + x212*x317 + x215*xi_xx + 12*x216*xi_x - x217*xi_xx + x233*x324 - x277*xi_xx + x278*xi_xx - x296*x324 + x321*x322 + x321*x323 + x322*x326 + x323*x326)) + x48*x90*(4*B1*Gt*x19*xi_y - 4*B1t*Fy*x185 + B2*Fy*x276*x96*xi_x + 64*B2*x275*x80 + Fxh*x174 - Fxh*x177 + 2*Fy*x188*xi_x + 4*Fyp*x11*x237 + Fyp*x198 + Fyt*x174 - Fyt*x177 - 8*G*u*xi_xy - 72*G*x155 - 2*Gc - Gh*x121 - 4*Gh*x19*x203 - Gh*x197 - Gh*x200 - 4*Gh*x202 + Gp*x146 + 4*Gp*x166*x19 - Gt*x173 - Gt*x176 + 4*Gt*x179 + Gt*x66 - 40*x0*x275 + x112*x213 + x123*x53 + x129*x175 - x141*x187 + x144*x199 + x145*x205 - x150*x187 + x151*x205 + x152*x206 + x153*x210 + x161*x276 + x164*x192 + x171*x178 + x181*x235 + x183*x201 + x184*xi_x + x207*xi_y + x211*x85 + x237*x273 - x238*x274)) + x47*(x48*(2*Fxph + 2*Fypt + x50*(B2*Fy*x233*xi_x + 8*B2*x145*x96 + 8*B2*x153*x80 + B2c*x17 + B2h*B2t*x218 - B2h*x220 + B2h*x221 - 8*B2p*Fx*Fy*x80 - B2p*x1*x165 - B2t*x223 + B2t*x224 - 16*Fx*Fy*x103*x265*x6 - 16*Fx*Fy*x11 - 8*Fx*Gh*x214 - 16*Fx*x267*xi_y - Fxh*x272 - Fxp*x143 - Fxp*x232*xi_y + (2*Fxp)*xi_y - 8*Fy*Gt*x214 - 12*Fy*phi*phit*x251*x80 - 4*Fy*x1*x250 + 16*Fy*x167*x80 - Fy*x19*x235*x6 + 8*Fy*x19*x257 - Fyp*x232*xi_x - Fyt*x272 - 8*Gh*x237*x80 - Gt*x175*x19 - 8*Gt*x238*x80 - phi*x103*x147*x248 - phih*phit*x252 - 12*phih*x254 - 12*phih*x258*xi_x - 12*phit*x259 - 6*u*x141 - 16*x11*x163*x265 + x11*x256*xi_x + 8*x11*x257*xi_y + x129*x222 + x130*x52 + x141*x243 - x141*x247 - x141*x253 + x141*x255 + x141*x261 - x141*x266 + x144*x219 + x150*x243 - x150*x247 - x150*x253 + x150*x255 + x150*x261 - x150*x266 + x156*x256 + x158*x234 + x158*x262 - x158*x268 + x162*x225 + x164*x225 + x166*x226 + x166*x260 - x166*x264 + x168 - x172*x267 - x215*xi_xy - 6*x216*xi_y + x217*xi_xy - x227*x228 + x227*x229 - x230*x231 - x235*x236 - x235*x249 + x240*x242 + x240*x244 - x242*x246 - x244*x246 + x54*x85)) - x90*(x21*(B1*Fy*x192*x96 + B2*Fy*x192*x96 + 8*B2h*Fy*x185 + 2*Fy*x188*xi_y + Fyh*x174 - Fyh*x177 - Gh*x173 - Gh*x176 + Gp*xi_yy - Gs + x144*x175 + x175*x66 + x178*x179 + x178*x183 + x178*x71 + x180*x181 + x182*x73 + x182*x76 + x184*xi_y - x186*x187 - x189*x91 + x190*x91 - x191*x69 + x193*x91 + x194*x69 + x195*x91 + x196*x69 - x51*xi_yy + x53*x68) + x36*(-64*B1*Fx*x214*xi_x + Fxp*x198 + Fxt*x0*(Gp - x51) + G*x11*x131 + 64*G*x164*xi_x - Gb + Gp*xi_xx - 8*Gt*x19*x203 - Gt*x197 - Gt*x200 - x121*x199 - x124*x191 - x124*x194 + x124*x196 - x127*x206 + x129*x199 + x130*x211 - x132*x189 + x132*x190 - x132*x193 + x132*x195 + x137*x213 + x171*x201 - x187*x212 - x201*x202 - x204*x205 + x205*x208 + x207*xi_x + x209*x210 - x51*xi_xx))))/9 - 2*x0*x8*(-x1*x23*x48*(Fx*x51*(x52*x62 + x60*xi_y + x61) + Fy*x51*(x57 + x60*xi_x) + 4*G*x0*(x25*xi_x + xi_y*(x38 + 18*x39 + x56)) + Gh*x44 + Gh*x45 + Gt*x29 + Gt*x33) + x1*x47*(x21*x34*x53 + x36*x46*(Gt + x51*(x54 + xi_x))) + x24*(Fyh*x63 - u*(B2h*Sh*x18 + 4*B2h*x1*x84*xi + Fy*(B1h*x109 + Fyp*x110*x50 + x108*x11*x67 + 2*x108*x59*xi_y + x111*x62 + x112*(-x114 + x115 + x119)) + Sh*x19*x74*xi_y - 18*Sh*x65 + 12*Sh*x71 - 3*Ss + x10*xi_yy + x25*x66 - x27*xi_yy + x28*x66 - x64*xi_yy - x67*x68 - x69*x70 - x69*x79 + x69*x81 + x69*x82 - x69*x83 + x69*x88 + x69*x89 + x72*x73 + x75*x76 - x77*x78 - x86*x87 + x90*x91*(-x100 - x102 + x105 + x106 + x12 + x94 + x95 + x98 - 3))) + x37*(Fxt*x63 + u*(12*B1*x125*xi_x - B2t*St*x18 - B2t*x19*x75*xi_x + Fx*(B1t*x109 + 2*B2t*(x133 - x134 + x135 - x136) + Fxp*x32*x50 + 2*St*x103*(x113 - x92) + x11*x122 - x137*(x114 - x115 + x119) + 2*x59*xi_x*(-x107 + x35 + 1)) + 3*Sb + 18*St*x120 - St*x19*x74*xi_x - x10*xi_xx + x121*x38 + x121*x42 + x122*x123 + x124*x70 + x124*x79 + x124*x81 - x124*x82 + x124*x83 + x124*x88 - x124*x89 + x126*x127 + x128*x78 - x129*x56 + x131*x87 + x132*x90*(x100 + x102 - x105 + x106 + x14 + x30 + x93 - x98) + x27*xi_xx + x64*xi_xx)) + x49*(Fxh*x138 + Fyt*x138 + x50*(B2*Fx*x169*x170 + 12*B2*Sh*x19*xi_x + 12*B2*x125*xi_y - Fx*Fyp*x154 + Fx*Sh*x149 - Fx*x61 - 4*Fx*x99*xi_y - Fxp*u*x34 + 6*Fy*Sp*x156 + Fy*St*x149 - Fy*x11*x55*x6*xi_x - 8*Fy*x164*x6 - 6*Fy*x54 - Fy*x57 - 2*Fyp*x1*x167 - Fyp*x11*x126 + Fyp*x159 + 72*S*x162*xi_y - 3*Sc - 9*Sh*x120 - St*x26 + x10*xi_xy - x101*x146 + x104*x172 - x111*x139 - x112*x41 - x116*x141 - x116*x150 + x117*x161 + x118*x141 + x118*x150 + x126*x152 + x129*x25 + x129*x28 - x133*x145 + x134*x145 + x134*x151 - x135*x145 - x135*x151 + x136*x145 + x136*x151 + x139*x143 + x140*x157 + x140 + x141*x7 - x141*x93 + x142*x157 + x142 + x144*x38 + x144*x42 - x147*x148 - x15*x163 - x150*x93 + x153*x154 - x155*x75 - x158*x79 + x158*x89 - 4*x159*xi*xi_y - x165*x99 - 2*x166*x19*x6*xi - x168*x87 + 16*x171*x6*xi*xi_y - x38*x68 + x54*x84 - x64*xi_xy)))/9 + 4*x19*(x24*x34^2 - x34*x49*(x43*x50 + 2*x45) + x37*x46^2)/9 + x9*(-B2*x17 + B2p)*(Sd*x18 - x7 + 3*(x4 + 1)^2)/9
    
    nothing
end


# this is another coupled equation, for B1d and Gd. the notation used is
#
# ( A11 d_uu B1d + A12 d_uu Gd + B11 d_u B1d + B12 d_u Gd + C11 B1d + C12 Gd ) = -S1
# ( A21 d_uu B1d + A22 d_uu Gd + B21 d_u B1d + B22 d_u Gd + C21 B1d + C22 Gd ) = -S2

function B1dGd_eq_coeff!(AA::Matrix, BB::Matrix, CC::Matrix, SS::Vector, vars::Tuple, ::Inner)
    (
        phi0, u, xi, xi_x, xi_y, xi_xx, xi_yy, xi_xy,
        B1     ,    B2     ,    G      ,    phi    ,    S      ,    Fx     ,    Fy     ,  Sd,
        B1p    ,    B2p    ,    Gp     ,    phip   ,    Sp     ,    Fxp    ,    Fyp    ,
        B1pp   ,    B2pp   ,    Gpp    ,    phipp  ,    Spp    ,    Fxpp   ,    Fypp   ,
        B1_x   ,    B2_x   ,    G_x    ,    phi_x  ,    S_x    ,    Fx_x   ,    Fy_x   ,
        B1_y   ,    B2_y   ,    G_y    ,    phi_y  ,    S_y    ,    Fx_y   ,    Fy_y   ,
        B1p_x  ,    B2p_x  ,    Gp_x   ,    phip_x ,    Sp_x   ,    Fxp_x  ,    Fyp_x  ,
        B1p_y  ,    B2p_y  ,    Gp_y   ,    phip_y ,    Sp_y   ,    Fxp_y  ,    Fyp_y  ,
        B1_xx  ,    B2_xx  ,    G_xx   ,    phi_xx ,    S_xx   ,
        B1_yy  ,    B2_yy  ,    G_yy   ,    phi_yy ,    S_yy   ,
                    B2_xy  ,    G_xy   ,                S_xy
    ) = vars

    @tilde_inner("B1")
    @tilde_inner("B2")
    @tilde_inner("G")
    @tilde_inner("phi")
    @tilde_inner("S")
    @tilde_inner("Fx")
    @tilde_inner("Fy")

    @hat_inner("B1")
    @hat_inner("B2")
    @hat_inner("G")
    @hat_inner("phi")
    @hat_inner("S")
    @hat_inner("Fx")
    @hat_inner("Fy")

    @bar_inner("B1")
    @bar_inner("B2")
    @bar_inner("G")
    @bar_inner("phi")
    @bar_inner("S")

    @star_inner("B1")
    @star_inner("B2")
    @star_inner("G")
    @star_inner("phi")
    @star_inner("S")

    @tilde_inner("Fxp")
    @tilde_inner("Fyp")

    @hat_inner("Fxp")
    @hat_inner("Fyp")

    @cross_inner("B2")
    @cross_inner("G")
    @cross_inner("S")

    x0 = u^2
    x1 = u^4
    x2 = B1*x1
    x3 = exp(x2)
    x4 = u*xi
    x5 = S*x1
    x6 = x4 - 1
    x7 = phi0^2
    x8 = x0*x7
    x9 = 3*x4 + 3*x5 + x6*x8 + 3
    x10 = -4*x0*x3*x9^4/27
    x11 = x9^3
    x12 = u*x11*x3
    x13 = 4*u*xi/3
    x14 = 10*S*x1/3
    x15 = u^3
    x16 = -2*Sp*x15/3
    x17 = G*x1
    x18 = tanh(x17)
    x19 = u^8
    x20 = u^5
    x21 = x4 + 1
    x22 = x21 + x5
    x23 = B1*u
    x24 = x15*x7
    x25 = x3*(3*u*x22 + x24*x6)^4/27
    x26 = x11*x3*(6*Sd*x1 + 3*x21^2 - x8)/9
    x27 = 4*x20/3
    x28 = sech(x17)
    x29 = B2*x1
    x30 = exp(x29)
    x31 = 3*Sh
    x32 = 9*S*u
    x33 = 2*x7*xi
    x34 = u*(x31 + x32*xi_y + x33*xi_y)
    x35 = 9*x5
    x36 = x35 + x8*(2*x4 - 1) - 3
    x37 = Fy*x36 + x34
    x38 = x30*x37^2
    x39 = exp(2*x2 + x29)
    x40 = 3*St
    x41 = x32*xi_x + x33*xi_x + x40
    x42 = u*x41
    x43 = Fx*x36
    x44 = x39*(x42 + x43)^2
    x45 = 2*x0*x9/3
    x46 = exp(x2 + x29)
    x47 = St*xi_y
    x48 = Sh*xi_x
    x49 = 4*G*u
    x50 = x8 + 3
    x51 = Fx*Gh
    x52 = Fy*Gt
    x53 = 3*Sp
    x54 = -x15*x53
    x55 = x36 + x54
    x56 = Fyh*x55
    x57 = 18*u
    x58 = Sh*xi_y
    x59 = 3*Sh*x1
    x60 = 6*Sh
    x61 = xi_y^2
    x62 = 2*x7
    x63 = 12*B2*x20
    x64 = 9*S*xi_y
    x65 = 36*S*x0
    x66 = u^6
    x67 = 36*B2*S*x66
    x68 = 6*u*x7*xi
    x69 = 2*x1*x7*xi*xi_y
    x70 = 4*Fyp
    x71 = x70*xi_y
    x72 = x0*x7*xi
    x73 = 8*B2*x20*x7*xi
    x74 = Fy^2
    x75 = x7*xi
    x76 = x15*x75
    x77 = 2*x29
    x78 = 2*x0*(x35 + x36*x77 + x54 + x76 + 3)
    x79 = u^7
    x80 = B2*x79
    x81 = x15 + 2*x80
    x82 = 4*x4
    x83 = 2*x0*x7*(x77 + 1)*(x82 - 1)
    x84 = 2*u*x36
    x85 = B2h*x15
    x86 = 12*x29
    x87 = 54*x5
    x88 = 6*Sp
    x89 = 72*B2*S*x19
    x90 = x15*x88 + x86 - x87 - x89 + 6
    x91 = u*(B2h*x20*x64 + B2h*x59 + B2h*x69 + Fy*(Fyp*x84 + x36*x85 + x60*x81 + x83*xi_y - x90*xi_y) + 18*Fyp*S*x15*xi_y + Fyp*x0*x60 + 3*Ss + x32*xi_yy + x33*xi_yy - x53*xi_yy + x57*x58 + x58*x63 + x61*x62 + x61*x65 + x61*x67 + x61*x68 + x61*x73 + x71*x72 + x74*x78)
    x92 = Fxt*x55
    x93 = St*xi_x
    x94 = 3*St*x1
    x95 = 6*St*x0
    x96 = xi_x^2
    x97 = 18*S*xi_x
    x98 = Fxp*x15
    x99 = 2*x1*x7*xi*xi_x
    x100 = 4*xi_x
    x101 = Fxp*x100
    x102 = Fx^2
    x103 = 6*St
    x104 = B2t*x15
    x105 = u*(9*B2t*S*x20*xi_x + B2t*x94 + B2t*x99 + Fx*(Fxp*x84 + x103*x81 + x104*x36 + x83*xi_x - x90*xi_x) + Fxp*x95 + 3*Sb + x101*x72 + x102*x78 + x32*xi_xx + x33*xi_xx - x53*xi_xx + x57*x93 + x62*x96 + x63*x93 + x65*x96 + x67*x96 + x68*x96 + x73*x96 + x97*x98)
    x106 = u*x9^2/3
    x107 = 2*x1
    x108 = 2*Fx*u
    x109 = 2*Fy*u
    x110 = B2t*x0
    x111 = Fyp*u
    x112 = 4*G*xi_x
    x113 = B2h*x0
    x114 = 4*B2*x20
    x115 = B2*x15
    x116 = 4*B2h*Fx
    x117 = B2t*Fy*x20
    x118 = 4*xi_y
    x119 = 4*Fx*Fyp
    x120 = Fy*x0
    x121 = 2*Fyph
    x122 = 2*u
    x123 = 2*B2p*u
    x124 = 12*u
    x125 = Fy*xi_y
    x126 = 8*B2*x0
    x127 = B2h*Fy
    x128 = 20*x1
    x129 = 2*Fyp
    x130 = 16*xi_y
    x131 = 4*x20
    x132 = 16*x15
    x133 = 32*B2*x19
    x134 = B2*B2h*x66
    x135 = B2*x66
    x136 = 112*B2*x20
    x137 = 8*x29
    x138 = 4*B2p
    x139 = Fy*x1*xi_y
    x140 = 72*B2*x79
    x141 = 40*B2*x15
    x142 = 4*B2p*x66
    x143 = B2^2
    x144 = u^9
    x145 = 128*x143*x144
    x146 = phi0^4
    x147 = Fy*x1*x146
    x148 = 64*u^11*x143
    x149 = 64*x143*x79
    x150 = 4*x20*x7
    x151 = phi0^6
    x152 = 4*x15*x151
    x153 = 24*phi*x146*x20
    x154 = 24*Fy
    x155 = phi*phih*x151*x66
    x156 = Fy*x146*x20*xi
    x157 = phi*phih*x1*x151
    x158 = phih*x146*x15*xi
    x159 = 24*phi*x146*x79
    x160 = 16*x66*x7*xi
    x161 = 96*phi*x146*x66*xi
    x162 = phi^2
    x163 = 72*x151*x162*x79
    x164 = xi^2
    x165 = 32*x164*x20*x7
    x166 = 48*phi*x146*x19*xi
    x167 = 48*phi*x1*x146*xi
    x168 = 36*x144*x151*x162
    x169 = 16*x164*x7*x79
    x170 = 36*x151*x162*x20
    x171 = 16*x15*x164*x7
    x172 = B2p*x15
    x173 = x137 - 2*x172 + 4
    x174 = u*(B2h^2*x131 + B2s*x122 - 8*Fy*Fyp*x135 + Fyh*x173 + Fyp^2*u - 8*Fyp*x120 - Fyp*x137*xi_y + phih^2*x152 + 8*phih*x147 - 16*phih*x156 + x113*x130 - x123*xi_yy + x124*x125 + x125*x136 + x125*x145 + x125*x153 - x125*x161 + x125*x163 + x125*x165 + x126*xi_yy + x127*x128 + x127*x133 - x129*x85 - x130*x158 + x132*x74 + 32*x134*xi_y - x138*x139 - 16*x139*x7*xi + x140*x74 + x141*x61 - x142*x74 + x148*x74 + x149*x61 + x150*x74 + x154*x155 + 24*x157*xi_y + x159*x74 - x160*x74 - x166*x74 - x167*x61 + x168*x74 + x169*x74 + x170*x61 + x171*x61 - x71)
    x175 = 2*Fxpt
    x176 = Fx*xi_x
    x177 = B2t*Fx
    x178 = 2*Fxp
    x179 = 16*xi_x
    x180 = 8*Fxp
    x181 = Fx*x0
    x182 = B2*B2t*x66
    x183 = B2*Fx*x66
    x184 = Fx*x1
    x185 = 8*phit
    x186 = Fx*x1*x146
    x187 = phi*phit*x151*x66
    x188 = x1*x7*xi
    x189 = phi*phit*x1*x151
    x190 = phit*x146*x15*xi
    x191 = u*(B2b*x122 - 4*B2p*x184*xi_x + B2t^2*x131 - 16*Fx*phit*x146*x20*xi + 24*Fx*x187 - 16*Fx*x188*xi_x + Fxp^2*u - Fxp*x137*xi_x + Fxt*x173 + phit^2*x152 - x101 + x102*x132 + x102*x140 - x102*x142 + x102*x148 + x102*x150 + x102*x159 - x102*x160 - x102*x166 + x102*x168 + x102*x169 - x104*x178 + x110*x179 - x123*xi_xx + x124*x176 + x126*xi_xx + x128*x177 + x133*x177 + x136*x176 + x141*x96 + x145*x176 + x149*x96 + x153*x176 - x161*x176 + x163*x176 + x165*x176 - x167*x96 + x170*x96 + x171*x96 - x179*x190 - x180*x181 - x180*x183 + 32*x182*xi_x + x185*x186 + 24*x189*xi_x)
    x192 = sinh(x17)
    x193 = x46*cosh(x17)
    x194 = 6*Fx
    x195 = x194*xi_y
    x196 = 6*Fy*xi_x
    x197 = B1h*Fx
    x198 = 3*x15
    x199 = B1t*Fy
    x200 = B2h*Fx
    x201 = B2t*Fy
    x202 = 12*Fx
    x203 = Fy*x15
    x204 = x7*xi_y
    x205 = 12*Fx*Sh
    x206 = B1*x79
    x207 = 12*Fx*xi_y
    x208 = 12*Fy*St
    x209 = 12*Fy*xi_x
    x210 = 12*B1*x20
    x211 = 9*S*x79
    x212 = 9*S*x20*xi_x
    x213 = B1t*x20
    x214 = Fx*xi_y
    x215 = Fy*xi_x
    x216 = B2t*x20
    x217 = Fx*Fy*x66
    x218 = Sp*x15
    x219 = Fyp*x15
    x220 = 72*S*xi_x*xi_y
    x221 = Fy*x7
    x222 = Fx*Fy
    x223 = x20*x7
    x224 = 4*x7*xi_y
    x225 = 4*B1*x66*xi_x
    x226 = 2*x66*x7*xi
    x227 = 4*B2*x66*xi_x
    x228 = 4*Fyp*xi_x
    x229 = B2*Fy
    x230 = 16*x7*xi*xi_y
    x231 = Fy*x7*xi
    x232 = 10*x1
    x233 = 8*xi_x
    x234 = 8*xi_y
    x235 = B1*x19
    x236 = B1*x66
    x237 = 4*Fxp*Fy
    x238 = Fxp*xi_y
    x239 = 4*x2
    x240 = 4*B2*x19
    x241 = 2*B1p*x1
    x242 = 16*B2*x19
    x243 = 4*x29
    x244 = xi_x*xi_y
    x245 = 2*B2p*x1
    x246 = 64*x143*x144
    x247 = phi*x146*x20
    x248 = 48*phi*x146*x66*xi
    x249 = 36*x151*x162*x79
    x250 = x164*x7
    x251 = B1p*x15
    x252 = -x172 + x243 + 2

    AA[1,1] = 0

    AA[1,2] = 0

    BB[1,1] = x10

    BB[1,2] = 0

    CC[1,1] = -x12*(16*G*S*x18*x19/9 + 16*G*x1*x18/9 + 16*G*x18*x20*xi/9 - 4*Gp*x15*x18*x22/9 + 2*x0*x7*(8*G*x1*x18*x6 - 2*Gp*x15*x18*x6 + 12*x4 - 9)/27 + x13 + x14 + x16 + 2/3)

    CC[1,2] = x18*x25*(4*B1p - 16*x23)

    SS[1] = -x106*x28*(x107*x46*(B2t*G*x118*x15 + Fxh*Gp - Fxh*x49 - Fxp*(Gh + x49*(x120 + xi_y)) + Fyp*Gt - Fyt*Gp + Fyt*x49 - G*x116*x20 + 4*G*x117 + G*x119*x15 + Gh*x108 + Gh*x110 + 4*Gh*x115*xi_x + Gp*x108*xi_y - Gp*x109*xi_x - Gt*x109 - Gt*x113 - 4*Gt*x115*xi_y + x111*x112 - x112*x85 + x114*x51 - x114*x52) + x30*(x121 - x174) + x39*(-x175 + x191)) + x26*(B1p - 4*x23) + x27*x28*(-x38 + x44) - x28*x45*(x1*x46*(-Fx*x49*(x15*x31 + x50*xi_y) + Fy*x49*(x15*x40 + x50*xi_x) + 12*G*x0*(x47 - x48) + Gh*x42 - Gt*x34 + x36*x51 + x52*(-x35 + x7*(x0 - 2*x15*xi) + 3)) - x30*(x56 + x91) + x39*(x105 + x92))

    AA[2,1] = 0

    AA[2,2] = 0

    BB[2,1] = 0

    BB[2,2] = x10

    CC[2,1] = x25*(-2*B1p + 8*x23)*sinh(2*x17)

    CC[2,2] = -x12*(2*x0*x7*(x82 - 3)/9 + x13 + x14 + x16 + 2/3)

    SS[2] = -x106*(x192*(-x121*x30 + x174*x30 - x175*x39 + x191*x39) + x193*(2*Fxph + 2*Fypt - x122*(B1*B2t*x118*x66 + B1h*x216 + B1h*x227 - B1h*x98 - 4*B1t*x135*xi_y + B1t*x219 + 56*B2*Fx*x20*xi_y + B2c*x122 - B2h*x213 + 4*B2h*x216 - B2h*x225 - B2h*x98 + 4*B2t*Fy*x235 - B2t*x219 + 24*Fx*Fy*phi*x146*x79 + 16*Fx*Fy*x164*x7*x79 - 8*Fx*phih*x146*x20*xi + 16*Fx*x20*x250*xi_y + 16*Fx*x203 - 16*Fx*x231*x66 + Fxh*(-x239 + x251 + x252) + Fxp*x111 - 4*Fxp*x120 + 12*Fy*x187 + 16*Fy*x20*x250*xi_x + Fyt*(x239 - x251 + x252) + 4*phih*phit*x15*x151 + 4*phih*x186 + 4*phit*x147 + u*x195 + u*x196 - 8*x1*x231*xi_x + x107*x197 - x107*x199 + x110*x234 + x113*x233 - x116*x235 - x119*x135 + x119*x236 - x123*xi_xy + x126*xi_xy - x129*xi_x + x130*x182 + x134*x179 - x135*x237 - x138*x217 + x140*x222 + x141*x244 + x148*x222 + x149*x244 + x150*x222 + 36*x151*x162*x20*xi_x*xi_y + x155*x202 - x156*x185 + 12*x157*xi_x - x158*x233 + 16*x164*x24*xi_x*xi_y - x166*x222 - x167*x244 + x168*x222 - x178*xi_y - x181*x70 - 8*x184*x7*xi*xi_y + 12*x189*xi_y - x190*x234 + x197*x240 - x199*x240 + x2*x228 + 56*x20*x229*xi_x + x200*x232 + x200*x242 + x201*x232 + x201*x242 + x207*x247 + x209*x247 + x214*x241 - x214*x245 + x214*x246 - x214*x248 + x214*x249 - x215*x241 - x215*x245 + x215*x246 - x215*x248 + x215*x249 - x228*x29 - x236*x237 - x238*x239 - x238*x243))) + x26*(Gp - x49) + x27*(x192*x38 + x192*x44 - x193*x37*(x122*x41 + 2*x43)) + x45*(-x192*(x105*x39 + x30*x56 + x30*x91 + x39*x92) + x193*(Fxh*x55 + Fyt*x55 + u*(-B1*Fx*x224*x66 + B1h*x212 + B1h*x94 + B1h*x99 - B1t*x59 - B1t*x69 + B2*Fx*x230*x79 + 72*B2*S*u^10*x222 - 8*B2*x19*x222*x7 + B2*x20*x230*xi_x + 16*B2*x231*x79*xi_x + B2h*x212 + B2h*x94 + B2h*x99 + B2t*x59 + B2t*x69 - 12*Fx*Fy*Sp*x20 + 18*Fx*Fyp*S*x20 - 2*Fx*Fyp*x24 + 16*Fx*x144*x229*x7*xi + 8*Fx*x76*xi_y + Fxp*x122*x37 - Fy*x0*x62*xi_x + 8*Fy*x15*x75*xi_x + Fyp*x95 + 36*S*x217 + S*x57*xi_xy + 6*Sc + Sh*x15*x194 + 12*u*x204*xi*xi_x - 2*x0*x214*x7 + x0*x220 + x100*x204 + x103*x203 - x111*x194 - x117*x7 + x119*x188 + x120*x202 + x135*x220 - x154*x183 - x183*x224 - x195*x218 - x195 - x196*x218 - x196 - x197*x198 + x197*x211 - x197*x223 + x197*x226 + x198*x199 - x198*x200 - x198*x201 - x199*x211 - x199*x226 - x2*x207 + x2*x209 + 4*x20*x222*x7*xi + x200*x211 - x200*x223 + x200*x226 + x201*x211 + x201*x226 - x205*x206 + x205*x80 + x206*x208 + x208*x80 + x210*x47 - x210*x48 + x213*x221 - x213*x64 - x214*x86 + x214*x87 + x214*x89 - x215*x86 + x215*x87 + x215*x89 + x216*x64 + x219*x97 + x221*x225 - x221*x227 + x228*x72 + x47*x57 + x47*x63 + x48*x57 + x48*x63 + 4*x7*xi*xi_xy - x88*xi_xy)))
    
    nothing
end

function phid_eq_coeff!(ABCS::Vector, vars::Tuple, ::Inner)
    (
        potential, phi0, u, xi, xi_x, xi_y, xi_xx, xi_yy, xi_xy,
        B1     ,    B2     ,    G      ,    phi    ,    S      ,    Fx     ,    Fy     ,  Sd,
        B1p    ,    B2p    ,    Gp     ,    phip   ,    Sp     ,    Fxp    ,    Fyp    ,
        B1pp   ,    B2pp   ,    Gpp    ,    phipp  ,    Spp    ,    Fxpp   ,    Fypp   ,
        B1_x   ,    B2_x   ,    G_x    ,    phi_x  ,    S_x    ,    Fx_x   ,    Fy_x   ,
        B1_y   ,    B2_y   ,    G_y    ,    phi_y  ,    S_y    ,    Fx_y   ,    Fy_y   ,
        B1p_x  ,    B2p_x  ,    Gp_x   ,    phip_x ,    Sp_x   ,    Fxp_x  ,    Fyp_x  ,
        B1p_y  ,    B2p_y  ,    Gp_y   ,    phip_y ,    Sp_y   ,    Fxp_y  ,    Fyp_y  ,
        B1_xx  ,    B2_xx  ,    G_xx   ,    phi_xx ,    S_xx   ,
        B1_yy  ,    B2_yy  ,    G_yy   ,    phi_yy ,    S_yy   ,
                    B2_xy  ,    G_xy   ,    phi_xy,     S_xy
    ) = vars

    @tilde_inner("B1")
    @tilde_inner("B2")
    @tilde_inner("G")
    @tilde_inner("phi")
    @tilde_inner("S")
    @tilde_inner("Fx")
    @tilde_inner("Fy")

    @hat_inner("B1")
    @hat_inner("B2")
    @hat_inner("G")
    @hat_inner("phi")
    @hat_inner("S")
    @hat_inner("Fx")
    @hat_inner("Fy")

    # @bar_inner("B1")
    # @bar_inner("B2")
    # @bar_inner("G")
    @bar_inner("phi")
    # @bar_inner("S")

    # @star_inner("B1")
    # @star_inner("B2")
    # @star_inner("G")
    @star_inner("phi")
    # @star_inner("S")

    # @tilde_inner("Fxp")
    # @tilde_inner("Fyp")

    # @hat_inner("Fxp")
    # @hat_inner("Fyp")

    # @cross_inner("B2")
    # @cross_inner("G")
    # @cross_inner("S")
    @cross_inner("phi")

    phiouter = u * phi0 - u^2 * phi0 * xi + u^3 * phi0^3 * phi

    x0 = u^2
    x1 = u^4
    x2 = B1*x1
    x3 = exp(x2)
    x4 = phi0^3
    x5 = u*xi
    x6 = x5 - 1
    x7 = x0*x6
    x8 = S*x1
    x9 = x5 + 1
    x10 = x8 + x9
    x11 = 3*x10
    x12 = phi0^2
    x13 = 3*x5
    x14 = x13 + 3
    x15 = x12*x7 + x14 + 3*x8
    x16 = 12*x5
    x17 = u^3
    x18 = Sp*x17
    x19 = x0*x12
    x20 = u^11
    x21 = S^3*x20
    x22 = 81*xi
    x23 = S^2
    x24 = 6*x5
    x25 = u^7
    x26 = 81*x25
    x27 = x10^2
    x28 = 4*x5
    x29 = x0*xi^2
    x30 = 3*x29
    x31 = S*x17
    x32 = x6^2
    x33 = 2*x5
    x34 = x33 - 1
    x35 = u^5
    x36 = phi0^6
    x37 = 3*x36
    x38 = x1*xi
    x39 = phi0^4
    x40 = 9*x39
    x41 = x17*x6
    x42 = 2*x6
    x43 = 2*x29
    x44 = x43 + 1
    x45 = u*x12
    x46 = 27*x10
    x47 = 9*S
    x48 = u*xi_y
    x49 = 2*xi
    x50 = x49*xi_y
    x51 = x19*x34 + 9*x8 - 3
    x52 = Fy*x51 + u*(3*Sh + x12*x50 + x47*x48)
    x53 = G*x1
    x54 = cosh(x53)
    x55 = x0*x49
    x56 = u - x55
    x57 = Fy*x0
    x58 = x57 + xi_y
    x59 = 3*phi
    x60 = u*x59
    x61 = phih + x58*x60
    x62 = Fy*x56 + x12*x61 - x50
    x63 = B2*x1
    x64 = exp(x63)
    x65 = x62*x64
    x66 = u*xi_x
    x67 = x49*xi_x
    x68 = Fx*x51 + u*(3*St + x12*x67 + x47*x66)
    x69 = Fx*x0
    x70 = x69 + xi_x
    x71 = Fx*x56 + x12*(phit + x60*x70) - x67
    x72 = exp(2*x2 + x63)
    x73 = x71*x72
    x74 = exp(x2 + x63)
    x75 = sinh(x53)
    x76 = x74*x75
    x77 = u^6
    x78 = 18*x77
    x79 = 162*x34
    x80 = phi0^8
    x81 = x17*xi^3
    x82 = 2*x1
    x83 = 6*x1
    x84 = Sd*x6
    x85 = Sd*x82
    x86 = -x13
    x87 = x86 - 3
    x88 = x29 + x6*x85 - 2*x8 + x81 + x87
    x89 = 9*phi
    x90 = x44 + x86
    x91 = x9^2
    x92 = x85 + x91
    x93 = phip*x11
    x94 = x10*x89
    x95 = 2*x81
    x96 = 4*x1*x84
    x97 = x14 - x43 + x8 - x95 - x96
    x98 = phi^3
    x99 = x6^3
    x100 = 2*UUp(phiouter,potential)
    x101 = x6^4
    x102 = x10^3
    x103 = phi*x11
    x104 = phi^2
    x105 = UU(phiouter,potential)*x104
    x106 = x32*x94
    x107 = -x32
    x108 = 4*x10
    x109 = 54*x27
    x110 = UU(phiouter,potential)*x6
    x111 = x108*x110
    x112 = u^9
    x113 = x104*x27
    x114 = 9*x113
    x115 = 24*phi
    x116 = 12*x10
    x117 = -phi*x116*x32
    x118 = 6*x6
    x119 = 54*x10
    x120 = u^10
    x121 = x6^6
    x122 = 4*x121
    x123 = x102*x98
    x124 = phi*x101
    x125 = x113*x32
    x126 = u^8
    x127 = 108*x123
    x128 = 36*x10
    x129 = x110*x116 - 1
    x130 = u*x58
    x131 = 4*G
    x132 = u*x131
    x133 = 4*xi
    x134 = xi_x*xi_y
    x135 = phih*xi_x
    x136 = phit*xi_y
    x137 = 6*phi
    x138 = 6*x66
    x139 = 2*x57
    x140 = xi_y^2
    x141 = Fy*x35
    x142 = Fy*x17
    x143 = Fy^2
    x144 = 4*B1
    x145 = Fy*xi_y
    x146 = x145*x77
    x147 = Fy*x77
    x148 = x147*x49
    x149 = 2*x38
    x150 = x149*xi_y
    x151 = 4*B2
    x152 = Fy*x149
    x153 = 20*xi
    x154 = x142*xi_y
    x155 = x55*xi_y
    x156 = 4*x126
    x157 = B1*x143
    x158 = B2*x156
    x159 = 14*xi
    x160 = x143*x35
    x161 = Fy*x25
    x162 = 16*xi
    x163 = x162*xi_y
    x164 = B2*x161
    x165 = 8*xi
    x166 = x112*x165
    x167 = x165*x35
    x168 = x140*x167
    x169 = B2*x143
    x170 = phih*x1
    x171 = 8*x142
    x172 = x17*x59
    x173 = phip*x0
    x174 = phih*x144
    x175 = x35*xi_y
    x176 = B1h*x59
    x177 = phih*x151
    x178 = B2h*x59
    x179 = 36*phi
    x180 = x1*x179
    x181 = 2*phip
    x182 = x115*x77
    x183 = 12*phi
    x184 = x0*x183
    x185 = x115*x126
    x186 = x145*x185
    x187 = x120*x183
    x188 = x183*x77
    x189 = x140*x188
    x190 = Fyp*x0
    x191 = xi_x^2
    x192 = Fx*x35
    x193 = Fx*x17
    x194 = 2*x69
    x195 = Fx^2
    x196 = Fx*x77
    x197 = x196*xi_x
    x198 = x196*x49
    x199 = x149*xi_x
    x200 = Fx*x149
    x201 = x193*xi_x
    x202 = Fxp*xi_x
    x203 = B1*x195
    x204 = x195*x35
    x205 = Fx*x25
    x206 = B1*x205
    x207 = x162*xi_x
    x208 = B2*x205
    x209 = B1*x191
    x210 = B2*x195
    x211 = B2*x191
    x212 = phit*x1
    x213 = 8*x193
    x214 = Fxp*x0
    x215 = 4*phit
    x216 = x35*xi_x
    x217 = x215*x216
    x218 = B1t*x59
    x219 = B2t*x59
    x220 = x192*x59
    x221 = Fx*x1
    x222 = x179*x221
    x223 = Fx*x126
    x224 = x115*xi_x
    x225 = B2*x223
    x226 = Fy*x192
    x227 = x193*xi_y
    x228 = x142*xi_x
    x229 = Fyp*xi_x
    x230 = B2*Fx*Fy
    x231 = B2*x134
    x232 = 6*u
    x233 = x151*x35
    x234 = Fy*xi_x

    ABCS[1] = 0

    ABCS[2] = -8*x0*x3*(phi0*x11 + x4*x7)^3/27

    ABCS[3] = -4*u*x15^2*x3*x4*(x16 - 9*x18 + x19*(10*x5 - 7) + 39*x8 + 3)/27

    ABCS[4] = 2*phi0*(18*x15*x17*(-x1*x54*x74*(Gh*x71 + Gt*x62 + x132*(Fx*u*(x139*(1 - x33) - x28*xi_y + x45*(phih + x130*x137) + xi_y) + Fy*u*(-x28*xi_x + x45*(phi*x138 + phit) + xi_x) + x12*x135 + x12*x136 - x133*x134 + x134*x137*x45)) + x1*x75*(x65*(Gh + x130*x131) + x73*(Gt + x132*x70)) + x54*x64*(B1*x161*x163 + B1*x168 - B1h*x141 + B1h*x148 + B1h*x150 - B2*x168 + B2h*x141 - B2h*x148 - B2h*x150 + Fyh*u - Fyh*x55 - Fyp*x142 + Fyp*x152 + Fyp*x155 + x12*(-B1*x186 - B1*x189 - B1h*x170 + B2*x186 + B2*x189 + B2h*x170 + Fyh*x172 - Fyh*x173 + phih*x171 + 6*phih*x48 - phip*xi_yy + phis + x140*x184 + x143*x182 + x145*x180 - x154*x181 - x157*x187 - x160*x181 - x161*x174 - x161*x176 + x161*x177 + x161*x178 + x169*x187 - x174*x175 - x175*x176 + x175*x177 + x175*x178 - x190*x61 + x60*xi_yy) + x139*xi_y - x140*x24 - 2*x140 + x143*x158 + x143*x83 - x144*x146 + x146*x151 - x153*x154 - x156*x157 + x157*x166 - x159*x160 - x163*x164 - x166*x169 - x49*xi_yy) + x54*x72*(B1t*x192 - B1t*x198 - B1t*x199 + B2t*x192 - B2t*x198 - B2t*x199 - Fxp*x193 + Fxp*x200 + Fxt*u - Fxt*x55 + x12*(B1*x217 + B1*x223*x224 + B1t*x212 + B2*x217 + B2t*x212 - Fxp*x220 + Fxt*x0*(-phip + x60) + phib - phip*xi_xx + phit*x138 + phit*x213 - phit*x214 - x172*x202 - x181*x201 - x181*x204 + x182*x195 + x184*x191 + x187*x203 + x187*x210 + x188*x209 + x188*x211 + x205*x218 + x205*x219 + x206*x215 + x208*x215 + x216*x218 + x216*x219 + x222*xi_x + x224*x225 + x60*xi_xx) + x144*x197 + x151*x197 - x153*x201 + x156*x203 + x158*x195 - x159*x204 - x166*x203 - x166*x210 - x167*x209 - x167*x211 - x191*x24 - 2*x191 + x194*xi_x + x195*x83 + x202*x55 - x206*x207 - x207*x208 - x49*xi_xx) - x76*(B2h*x192 - B2h*x198 - B2h*x199 + B2t*x141 - B2t*x148 - B2t*x150 + Fxh*u - Fxh*x55 - Fxp*x142 + Fxp*x152 + Fxp*x155 + 12*Fy*x221 + 8*Fy*x225 - Fyp*x193 + Fyp*x200 + Fyt*u - Fyt*x55 - x112*x162*x230 + x12*(B2*x185*x234 + B2h*x212 + B2t*x170 + Fxh*x172 - Fxh*x173 + 48*Fy*phi*x196 - Fyp*x220 + Fyt*x172 - Fyt*x173 + phi*x232*xi_xy + 2*phic + 4*phih*x208 + phih*x213 - 4*phip*x226 + phit*x171 - phit*x190 + x0*x115*x134 + x115*x120*x230 + x115*x225*xi_y + x135*x232 + x135*x233 + x136*x232 + x136*x233 + x161*x219 + x164*x215 - x172*x229 + x175*x219 + x178*x205 + x178*x216 + x180*x234 - x181*x227 - x181*x228 - x181*xi_xy + x182*x231 - x214*x61 + x222*xi_y) - x133*xi_xy - x134*x16 - 4*x134 + x139*xi_x + x147*x151*xi_x + x151*x196*xi_y - x153*x227 - x153*x228 - x162*x231*x35 - x163*x208 - x164*x207 + x194*xi_y - 28*x226*xi + x229*x55)) + x3*(-81*Sp*x0*x27 + 243*x21 - x22*(x5 + 2) + x23*x26*(x24 + 5) + 81*x31*(x28 + x30 + 1) + x32*x34*x35*x37 + x40*x41*(Sp*(x17 - x38) + 4*x29 + x6 + x8*(7*x5 - 5)) + x45*x46*(-x18*x42 + x44 - x5 + x8*(8*x5 - 7))) - x3*(-Sd*x17*x27*x79 + x0*x40*(u*(-x88*x90 - x94*x97) + x93*x97) + x1*x37*x6*(phip*(-x30 + 9*x5 + 6*x8 - 3*x81 - x83*x84 + 9) + u*(x88*x89 + x90)) - x12*x46*(u*x34*(x43 - x8 + x87 + x95 + x96) - u*x92*x94 + x92*x93) - x22*(8*x29 + x33 + 7*x81 + x82*xi^4 - 2) - x26*x34*(S*x5 + S)^2 - x31*x79*x9^3 + 3*x32*x77*x80*(phip - x60)) - x3*(8*UU(phiouter,potential)*phi0^14*u^13*x98*x99 + UUp(phiouter,potential)*phi0^13*u^12*x104*x118*(2*x101 + x114 + x117) - UUp(phiouter,potential)*phi0^9*x126*x42*(-x121 + x124*x128 - 162*x125 + x127) + UUp(phiouter,potential)*phi0^7*x10*x32*x78*(x101 + 18*x113 + x117) - UUp(phiouter,potential)*phi0^5*x1*x109*x99*(phi*x108 + x107) + 54*UUp(phiouter,potential)*x0*x101*x102*x4 + phi^4*phi0^17*u^16*x100*x99 + phi*phi0^11*x100*x120*(x119*x124 - x122 + 27*x123 - 108*x125) + phi0^15*u^14*x100*x98*(-4*x101 + x106) + 24*phi0^12*x105*x20*(-x101 + x103*x32) + phi0^10*x110*x112*x115*(x101 - x106 + x114) - x109*x45*(x103 + x32*(x111 - 3)) - x118*x35*x36*(108*x102*x105 - x129*x6*x94 + x129*x99) + x119*x39*x41*(x103 + x107)*(x111 - 1) + 162*x21*x6 + 486*x23*x25*(x29 - 1) + 2*x25*x80*(-UU(phiouter,potential)*x122 - 324*UU(phiouter,potential)*x125 + UU(phiouter,potential)*x127 + x59*x99*(x110*x128 - 1)) + 486*x31*x6*x91 + 162*xi*(x43 + x81 - 2)) + x78*(x52*x54*x65 + x54*x68*x73 - x76*(x52*x71 + x62*x68)))/27


    nothing
end


function A_eq_coeff!(ABCS::Vector, vars::Tuple, ::Inner)
    (
        potential, phi0, u, xi, xi_x, xi_y, xi_xx, xi_yy, xi_xy,
        B1   , B2   , G   , phi   , S    , Fx    , Fy    , Sd, B1d, B2d, Gd, phid,
        B1p  , B2p  , Gp  , phip  , Sp   , Fxp   , Fyp   ,
        B1pp , B2pp , Gpp , phipp , Spp  , Fxpp  , Fypp  ,
        B1_x , B2_x , G_x , phi_x , S_x  , Fx_x  , Fy_x  ,
        B1_y , B2_y , G_y , phi_y , S_y  , Fx_y  , Fy_y  ,
        B1p_x, B2p_x, Gp_x, phip_x, Sp_x , Fxp_x , Fyp_x ,
        B1p_y, B2p_y, Gp_y, phip_y, Sp_y , Fxp_y , Fyp_y ,
        B1_xx, B2_xx, G_xx, phi_xx, S_xx ,
        B1_yy, B2_yy, G_yy, phi_yy, S_yy ,
               B2_xy, G_xy, phi_xy, S_xy
    ) = vars

    @tilde_inner("B1")
    @tilde_inner("B2")
    @tilde_inner("G")
    @tilde_inner("phi")
    @tilde_inner("S")
    @tilde_inner("Fx")
    @tilde_inner("Fy")

    @hat_inner("B1")
    @hat_inner("B2")
    @hat_inner("G")
    @hat_inner("phi")
    @hat_inner("S")
    @hat_inner("Fx")
    @hat_inner("Fy")

    @bar_inner("B1")
    @bar_inner("B2")
    @bar_inner("G")
    @bar_inner("phi")
    @bar_inner("S")

    @star_inner("B1")
    @star_inner("B2")
    @star_inner("G")
    @star_inner("phi")
    @star_inner("S")

    @tilde_inner("Fxp")
    @tilde_inner("Fyp")

    @hat_inner("Fxp")
    @hat_inner("Fyp")

    @cross_inner("B2")
    @cross_inner("G")
    @cross_inner("S")
    @cross_inner("phi")

    phiouter = u * phi0 - u^2 * phi0 * xi + u^3 * phi0^3 * phi

    x0 = u^4
    x1 = B1*x0
    x2 = exp(x1)
    x3 = u^3
    x4 = u*xi
    x5 = x4 - 1
    x6 = phi0^2
    x7 = x5*x6
    x8 = S*x0
    x9 = x4 + 1
    x10 = x8 + x9
    x11 = 3*x10
    x12 = 3*x4
    x13 = 3*x8
    x14 = u^2
    x15 = x12 + x13 + x14*x7 + 3
    x16 = 4*x15^4*x2/9
    x17 = u^14
    x18 = S^4*x17
    x19 = u^6
    x20 = x19*(S*x4 + S)^2
    x21 = S*x14
    x22 = x21*x9^3
    x23 = S^3
    x24 = u^10
    x25 = x23*x24*x9
    x26 = x5^4
    x27 = phi0^8
    x28 = 4*x19*x27
    x29 = xi^2
    x30 = 4*x4
    x31 = x14*x29
    x32 = x10^3
    x33 = x5^3
    x34 = phi0^6
    x35 = x0*x34
    x36 = x5^2
    x37 = phi0^4
    x38 = x10^2
    x39 = x14*x38
    x40 = 3*Sh
    x41 = 9*S
    x42 = u*x41
    x43 = 2*x6
    x44 = x43*xi
    x45 = u*(x40 + x42*xi_y + x44*xi_y)
    x46 = 9*x8
    x47 = 2*x4
    x48 = x47 - 1
    x49 = x14*x6
    x50 = x46 + x48*x49 - 3
    x51 = Fy*x50
    x52 = x45 + x51
    x53 = G*x0
    x54 = cosh(x53)
    x55 = B2*x0
    x56 = exp(x55)
    x57 = x54*x56
    x58 = 3*St
    x59 = 2*xi_x
    x60 = x6*xi
    x61 = u*(x42*xi_x + x58 + x59*x60)
    x62 = Fx*x50
    x63 = x61 + x62
    x64 = 2*x1
    x65 = exp(x55 + x64)
    x66 = x54*x65
    x67 = sinh(x53)
    x68 = exp(x1 + x55)
    x69 = x67*x68
    x70 = 324*x14
    x71 = xi^3
    x72 = u*x71
    x73 = 1944*xi
    x74 = S*x3
    x75 = S^2
    x76 = x19*x75
    x77 = 972*x23
    x78 = xi^4
    x79 = x29*x8
    x80 = u^5
    x81 = S*x71
    x82 = x80*x81
    x83 = S*x19
    x84 = x78*x83
    x85 = u^7
    x86 = x85*xi
    x87 = x75*x86
    x88 = u^11
    x89 = u^8
    x90 = x75*x89
    x91 = u^9
    x92 = 1944*x71
    x93 = 2*x0
    x94 = Sd*x93
    x95 = u*x38
    x96 = x3*x71
    x97 = 7*x29
    x98 = x78*x93
    x99 = Sp*x3
    x100 = x80*xi
    x101 = -x12
    x102 = 2*x31
    x103 = x102 + 1
    x104 = -x47
    x105 = x104 + 2
    x106 = 12*x5
    x107 = 3*x96
    x108 = Sd*x5
    x109 = 4*x0
    x110 = x101 - 3
    x111 = x19*x29
    x112 = 8*x85
    x113 = -x4
    x114 = 8*x4
    x115 = -x31 - 2
    x116 = 2*xi
    x117 = x14*x78
    x118 = xi^5
    x119 = 4*x3
    x120 = 4*x31
    x121 = 2*x14
    x122 = 36*x0
    x123 = Fy*x14
    x124 = 4*u
    x125 = G*x124
    x126 = x0*x67
    x127 = 18*u
    x128 = S*x127
    x129 = 4*x60
    x130 = 4*x14
    x131 = x30 - 1
    x132 = x131*x49 + 18*x8 - 3
    x133 = x54*x68
    x134 = 3*x3
    x135 = Sp*x134
    x136 = -x135
    x137 = x136 + x50
    x138 = 3*Sp
    x139 = Sh*x127
    x140 = x0*x40
    x141 = xi_y^2
    x142 = 12*x80
    x143 = Sh*x142
    x144 = B1*xi_y
    x145 = x41*x80
    x146 = x145*xi_y
    x147 = B2*x143
    x148 = 36*x21
    x149 = 36*B1
    x150 = x141*x83
    x151 = 36*B2
    x152 = 6*x4
    x153 = x152*x6
    x154 = x60*x93
    x155 = x154*xi_y
    x156 = 8*B1
    x157 = x6*x80
    x158 = x157*xi
    x159 = x141*x158
    x160 = 8*B2
    x161 = Fy^2
    x162 = 27*x8
    x163 = 5*x4
    x164 = x48*x64
    x165 = 2*x55
    x166 = S*x89
    x167 = B2*x166
    x168 = 18*x167
    x169 = -6*x1
    x170 = 18*B1*x166 + x169
    x171 = 6*x55
    x172 = x135 + x171 + 3
    x173 = 9*x3
    x174 = 6*B1
    x175 = 6*B2
    x176 = 2*Sh
    x177 = x3*x50
    x178 = B2h*x177
    x179 = x41*x85
    x180 = x19*x44
    x181 = 45*x8
    x182 = x131*x64
    x183 = -x30
    x184 = 36*x167
    x185 = x149*x166 + x169
    x186 = 2*xi_y
    x187 = St*xi_x
    x188 = x0*x58
    x189 = xi_x^2
    x190 = x142*x187
    x191 = B1t*xi_x
    x192 = x145*xi_x
    x193 = x189*x83
    x194 = x154*xi_x
    x195 = x158*x189
    x196 = Fx^2
    x197 = x163 + x165*x48 - 2
    x198 = x136 - x171 - 3
    x199 = x114 + x131*x165 - 1
    x200 = x181 + x184 + x198
    x201 = Fy*xi_x
    x202 = 6*x201
    x203 = St*xi_y
    x204 = B2t*Fy
    x205 = Fy*St
    x206 = xi_x*xi_y
    x207 = x206*x6
    x208 = B2*x85
    x209 = x201*x55
    x210 = 72*x206
    x211 = 72*x201
    x212 = 12*x4
    x213 = 4*Fy
    x214 = B2*x19
    x215 = 16*x3
    x216 = x215*x60
    x217 = x201*x216
    x218 = 16*B2
    x219 = x201*x6
    x220 = x0*x6
    x221 = 11664*x29
    x222 = 5832*x19
    x223 = B2*B2d
    x224 = B2d*B2p
    x225 = 1458*x224
    x226 = G*Gd
    x227 = 1944*x19
    x228 = Gd*Gp
    x229 = 486*x80
    x230 = 23328*x223
    x231 = 7776*x226
    x232 = x19*x73
    x233 = x71*x91
    x234 = x24*x78
    x235 = B2d*x85
    x236 = x78*x91
    x237 = x221*x89
    x238 = 1944*x234
    x239 = x228*x85
    x240 = 2916*x29
    x241 = x228*x89
    x242 = 486*x236
    x243 = B1d*x54^2
    x244 = B1*x243
    x245 = B1p*x243
    x246 = 7776*x244
    x247 = x243*x85
    x248 = B1p*x247
    x249 = x245*x89
    x250 = phi^4
    x251 = 8*UU(phiouter,potential)
    x252 = phi^3
    x253 = 4*x89
    x254 = x226*x253
    x255 = x235*(-B2*x124 + B2p)
    x256 = 3*x255
    x257 = B1*x124
    x258 = x239 - x254
    x259 = x247*(B1p - x257) + x256 + x258 + 4
    x260 = 1944*x259
    x261 = phi^2
    x262 = x5^6
    x263 = 8*x10
    x264 = x261*x38
    x265 = x264*x36
    x266 = phi*u
    x267 = 6*phid
    x268 = 216*x252
    x269 = UU(phiouter,potential)*x32
    x270 = UU(phiouter,potential)*x10
    x271 = x270*x5
    x272 = 3*phi
    x273 = x272*x33
    x274 = 8*x226
    x275 = x91*xi
    x276 = 6*x5
    x277 = 8*x244
    x278 = x47 - 2
    x279 = phid*x10
    x280 = 24*x279
    x281 = x113 + 1
    x282 = 3*phip
    x283 = x10^4
    x284 = UU(phiouter,potential)*x283
    x285 = -x152
    x286 = 2*UU(phiouter,potential)
    x287 = 36*x10
    x288 = x270*x36
    x289 = x104 + 5
    x290 = 4*x91
    x291 = 9*x279
    x292 = x10*x5
    x293 = 12*x292
    x294 = 12*UU(phiouter,potential)
    x295 = x10*x36
    x296 = UU(phiouter,potential)*x33
    x297 = x115 + x30
    x298 = x38*x48
    x299 = 4*x96
    x300 = x0*x78
    x301 = x118*x80
    x302 = x296*x9
    x303 = x33*x90
    x304 = x290*xi
    x305 = x226*x304 - x241*xi - x244*x253 + x244*x304 + x248 - x249*xi - x256*x5 + x258
    x306 = 18*x38
    x307 = x261*x32
    x308 = phid*x298
    x309 = 24*x96
    x310 = UU(phiouter,potential)*x31
    x311 = x3*x34
    x312 = 72*UU(phiouter,potential)
    x313 = 2*x57
    x314 = 2*x67
    x315 = x314*x56
    x316 = x313*xi_yy
    x317 = x156*x57
    x318 = u*xi_yy
    x319 = 16*u
    x320 = B1h*x57
    x321 = x320*xi_y
    x322 = x160*x57
    x323 = x319*xi_y
    x324 = B2h*x57
    x325 = Fy*x124
    x326 = G*x67
    x327 = 8*x326
    x328 = x327*x56
    x329 = Gh*x56
    x330 = x329*x67
    x331 = Fyh*x3
    x332 = x57*x93
    x333 = Fy*x215
    x334 = 4*Gh
    x335 = B1h*x56
    x336 = Fyh*x121
    x337 = B1p*x57
    x338 = x67*x93
    x339 = B2p*x57
    x340 = Gp*x67
    x341 = x340*x56
    x342 = x161*x57
    x343 = B1*Fy
    x344 = x343*x85
    x345 = 16*x320
    x346 = x144*x80
    x347 = x112*x343
    x348 = 8*B2h
    x349 = 16*x330
    x350 = 96*xi_y
    x351 = Fy*x57
    x352 = x350*x351
    x353 = B2*Fy*x112
    x354 = x160*x80
    x355 = 16*x85
    x356 = Fy*x326
    x357 = 16*x80
    x358 = x357*xi_y
    x359 = Fy*x119*xi_y
    x360 = 32*x324
    x361 = x80*xi_y
    x362 = x354*xi_y
    x363 = x356*x56
    x364 = G*x355
    x365 = Gh*x364
    x366 = x53*x67
    x367 = G*Gh
    x368 = 4*x69
    x369 = 4*x133
    x370 = 56*x342
    x371 = B1*x19
    x372 = x141*x57
    x373 = 40*x14
    x374 = B1*x373
    x375 = 4*x80
    x376 = x342*x375
    x377 = B2*x373
    x378 = 56*x19
    x379 = x161*x56
    x380 = x326*x379
    x381 = x340*x375
    x382 = x326*x373
    x383 = x141*x56
    x384 = 64*x89
    x385 = B2*x384
    x386 = 128*x89
    x387 = 2*x66
    x388 = x314*x65
    x389 = Fyp*x69
    x390 = x24*x342
    x391 = 32*B1
    x392 = B2*x391
    x393 = x214*x391
    x394 = x24*x380
    x395 = 64*B1
    x396 = x326*x383
    x397 = 64*x371
    x398 = B1^2
    x399 = x384*x398
    x400 = x351*xi_y
    x401 = 32*B2
    x402 = 32*x214
    x403 = B2^2
    x404 = 128*x403
    x405 = x404*x89
    x406 = G^2
    x407 = 64*x406
    x408 = x407*x89
    x409 = x3*x37
    x410 = phih*x409
    x411 = x387*xi_xx
    x412 = 32*x398
    x413 = x19*x412
    x414 = x319*x69
    x415 = 64*x403
    x416 = x19*x415
    x417 = B2h*xi_x
    x418 = B2t*xi_y
    x419 = Fx*x124
    x420 = 32*x406
    x421 = 4*x220
    x422 = G*x133
    x423 = x19*x420
    x424 = Gh*x133
    x425 = x319*xi_x
    x426 = Gt*x133
    x427 = x130*x34
    x428 = phi*x37
    x429 = 24*x428
    x430 = x0*x429
    x431 = phih*x34
    x432 = 24*phi
    x433 = x432*x80
    x434 = x431*x433
    x435 = 16*x37*xi
    x436 = x0*x435
    x437 = x311*x432
    x438 = x437*xi_y
    x439 = phih*x57
    x440 = x14*x435
    x441 = x156*x66
    x442 = u*xi_xx
    x443 = x191*x66
    x444 = x160*x66
    x445 = B2t*x66
    x446 = x327*x65
    x447 = Gt*x65
    x448 = x447*x67
    x449 = x19*x428
    x450 = 24*x449
    x451 = 16*x158
    x452 = x100*x428
    x453 = x426*x93
    x454 = x424*x93
    x455 = Fxh*x3
    x456 = x160*x69
    x457 = Fyt*x3
    x458 = x126*x68
    x459 = x215*x69
    x460 = B2h*Fx
    x461 = Fxh*x121
    x462 = B2p*x69
    x463 = Fyt*x121
    x464 = Fx*x69
    x465 = 8*x464
    x466 = Fx*x215
    x467 = 8*x422
    x468 = Gp*x133
    x469 = x261*x34
    x470 = x19*x469
    x471 = 72*x470
    x472 = 32*x29
    x473 = x220*x472
    x474 = x428*x86
    x475 = 48*x474
    x476 = phi*x409*xi
    x477 = 48*x476
    x478 = Fxt*x3
    x479 = Fxt*x121
    x480 = B1p*x66
    x481 = x66*x93
    x482 = x466*x66
    x483 = B2p*x66
    x484 = x340*x65
    x485 = Fx*x112
    x486 = x424*x485
    x487 = x80*xi_x
    x488 = x424*x487
    x489 = Fy*x422
    x490 = 32*x69
    x491 = x208*x490
    x492 = x417*x80
    x493 = B2*x490
    x494 = x418*x80
    x495 = Fy*x464
    x496 = x464*xi_y
    x497 = 96*x55
    x498 = x206*x69
    x499 = 80*x14
    x500 = x112*x422
    x501 = 8*Fy
    x502 = Fx*x80
    x503 = x501*x502
    x504 = Fx*x119
    x505 = x504*xi_y
    x506 = x119*x201
    x507 = Fx*x19
    x508 = x133*x53
    x509 = Gt*x69
    x510 = x469*x89
    x511 = 36*x510
    x512 = 16*x6
    x513 = x111*x512
    x514 = x357*xi_x
    x515 = x206*x422
    x516 = x122*x469
    x517 = x31*x512
    x518 = x196*x66
    x519 = B1t*x66
    x520 = Fx*x355
    x521 = B1*x520
    x522 = B2t*x487
    x523 = Fx*x66
    x524 = x523*xi_x
    x525 = 96*x524
    x526 = x504*xi_x
    x527 = B2*x485
    x528 = x326*x65
    x529 = 32*x445
    x530 = Fx*xi_x
    x531 = x385*x422
    x532 = 56*x518
    x533 = x189*x66
    x534 = x375*x518
    x535 = x196*x65
    x536 = x326*x535
    x537 = x189*x65
    x538 = x24*x495
    x539 = x201*x69
    x540 = x19*x498
    x541 = Fx*x458
    x542 = x541*x6
    x543 = phit*x69
    x544 = x14*x543
    x545 = x528*x530
    x546 = x24*x518
    x547 = x24*x536
    x548 = x326*x537
    x549 = phit*x523
    x550 = 32*x495
    x551 = phih*x69*xi_x
    x552 = phit*x66
    x553 = x552*xi_x
    x554 = 96*x452

    ABCS[1] = 2*x2*(u*x11 + x3*x7)^4/27

    ABCS[2] = x16*x3

    ABCS[3] = x14*x16

    ABCS[4] = -x0*x15^2*(B1*x357*x443 + B1*x386*x545 + B1*x445*x485 + B1*x448*x514 - B1*x486 + B1b*x387 + B1h^2*x332 - B1h*B2h*x332 + B1h*x422*x485 + B1h*x453 + B1h*x467*x487 + B1p*x316 + B1p*x376 - B1p*x411 - B1p*x534 - B1s*x313 + B1t^2*x481 + B1t*B2t*x481 - B1t*x112*x489 + 4*B1t*x126*x447 - B1t*x361*x467 - B1t*x454 + B1t*x482 + B1t*x520*x528 - 64*B2*Fx*x24*x489 + B2*x360*x361 + B2*x395*x524*x89 - B2*x414*xi_xy - B2*x486 + B2*x487*x529 - B2*x498*x499 + B2b*x387 - B2c*x368 + B2h^2*x109*x57 + B2h*x112*x363 + B2h*x329*x338 - B2h*x453 - B2p*x316 + B2p*x368*xi_xy - B2p*x376 - B2p*x411 - B2p*x534 + B2s*x313 + B2t^2*x109*x66 + B2t*x338*x447 - B2t*x348*x458 - B2t*x454 + B2t*x482 + B2t*x485*x528 + Fx*x208*x529 - Fx*x350*x508 - Fx*x531*xi_y - Fxp^2*x66 - Fxp*x325*x69 + 2*Fxp*x389 + Fxp*x419*x66 + Fy*phit*x435*x458 + Fy*x208*x360 - Fy*x34*x433*x543 + Fy*x350*x366*x56 - Fy*x364*x509 - Fyp^2*x57 + Fyp*x325*x57 + G*Gt*x514*x66 - G*x358*x509 + Gb*x388 - Gc*x369 + Gh^2*x332 - Gp*x315*xi_yy + Gp*x369*xi_xy - Gp*x388*xi_xx + Gs*x315 + Gt^2*x481 - Gt*x334*x458 + Gt*x364*x523 + phih^2*x427*x57 - phih*x351*x436 + phih*x435*x541 + phit^2*x427*x66 - x1*x352 + x1*x525 - x111*x550*x6 + x123*x465 - x126*x334*x335 - x130*x342 - x130*x518 - x144*x351*x385 - x144*x363*x386 - x156*x488 + x158*x550 + x160*x448*x487 - x160*x488 + x191*x357*x528 + x191*x444*x80 - x201*x429*x458 - 96*x201*x508 - x201*x531 - x204*x459 - x204*x491 - x204*x500 - x207*x31*x490 - 96*x209*x69 - x210*x458*x469 - x211*x470*x69 + x214*x370 - 112*x214*x495 - 64*x214*x515 + x214*x532 - x216*x400 + x216*x496 - x216*x524 + x217*x69 - x219*x458*x472 - x317*x318 - x317*x331 + x318*x322 + x318*x328 - x319*x321 - x319*x422*xi_xy + x319*x443 - x320*x333 - x320*x353 - x321*x354 + x322*x331 + x323*x324 + x323*x330 - x323*x426 + x324*x333 - x324*x347 - x326*x335*x358 + x326*x348*x361*x56 + x328*x331 + x330*x333 + x330*x353 + x330*x362 - x333*x426 - x335*x355*x356 + x336*x337 - x336*x339 - x336*x341 + x337*x359 - x339*x359 + x34*x432*x502*x552 - x341*x359 + x342*x421 + x342*x450 - x342*x451 - x342*x475 + x342*x511 + x342*x513 + x344*x345 - x344*x349 + x345*x346 - x346*x348*x57 - x346*x349 + 8*x346*x426 + x347*x426 + x351*x365 + 8*x351*x410 + x351*x434 - x352*x452 + x352*x55 - x353*x426 + x358*x367*x57 - x362*x426 + x363*x385*xi_y - x365*x464 + 96*x366*x530*x65 - x367*x514*x69 - x370*x371 + x371*x532 - x372*x374 + x372*x377 - x372*x393 + x372*x413 + x372*x416 + x372*x423 - x372*x477 + x372*x516 + x372*x517 + x374*x533 + x377*x533 + x378*x380 + x378*x536 - x379*x381 - x381*x535 + x382*x383 + x382*x537 + x385*x545 - x389*x419 - x390*x392 + x390*x412 + x390*x415 + x390*x420 + x392*x546 + x393*x533 - x394*x395 + x394*x401 + x395*x547 - x396*x397 + x396*x402 + x397*x548 + x399*x400 + x399*x524 + x400*x405 + x400*x408 + x400*x430 + x400*x471 + x400*x473 + x401*x547 + x402*x548 - x404*x538 - x404*x540 - x405*x496 + x405*x524 - x405*x539 - x407*x538 - x407*x540 - x408*x496 + x408*x524 - x408*x539 - x409*x501*x543 + 8*x409*x549 - x410*x465 + x412*x546 + x413*x533 - x414*x417 - x414*x418 + x415*x546 + x416*x533 + x420*x546 + x421*x518 + x423*x533 - x424*x425 - x424*x466 + x425*x445 + x425*x448 - x429*x541*xi_y + x430*x524 - 8*x431*x544 - x434*x464 + x435*x544*xi_y - x436*x549 - x437*x551 + x437*x553 + x438*x439 - x438*x543 - x439*x440*xi_y + x440*x551 - x440*x553 + x441*x442 + x441*x478 + x441*x522 + x442*x444 + x442*x446 + x444*x478 + x446*x478 + x446*x522 + x448*x466 + x448*x521 + x448*x527 - 48*x449*x495 + x450*x518 - x451*x518 - x452*x525 - x455*x456 - x455*x467 - x456*x457 - x457*x467 - x459*x460 - x460*x491 - x460*x500 + x461*x462 + x461*x468 + x462*x463 + x462*x503 + x462*x505 + x462*x506 + x463*x468 - x467*x492 - x467*x494 + x468*x503 + x468*x505 + x468*x506 + 72*x469*x507*x66*xi_x - x471*x496 - x472*x542*xi_y + x473*x524 + 96*x474*x495 - x475*x518 + 96*x476*x498 - x477*x533 - x479*x480 - x479*x483 - x479*x484 - x480*x526 - x483*x526 - x484*x526 - 112*x489*x507 - x492*x493 - x493*x494 - 72*x495*x510 - x496*x497 + x496*x554 + x497*x524 - x499*x515 - x501*x542 + x511*x518 + x513*x518 + x516*x533 + x517*x533 + x519*x521 + x519*x527 + x539*x554)/3 - 8*x15*x3*(-x0*x133*(Fx*x125*(x121*x51 + x132*xi_y + x3*x40) + Fy*x125*(x132*xi_x + x3*x58) + G*x130*(x40*xi_x + xi_y*(x128*xi_x + x129*xi_x + x58)) + Gh*x61 + Gh*x62 + Gt*x45 + Gt*x51) + x126*(x52*x56*(Gh + x125*(x123 + xi_y)) + x63*x65*(Gt + x125*(Fx*x14 + xi_x))) + x57*(Fyh*x137 + u*(-B1h*x140 - B1h*x146 - B1h*x155 + B2h*x140 + B2h*x146 + B2h*x155 + Fy*(B1h*(x134 + x157 - x179 - x180) - x176*(-x173 + x85*(x174 - x175)) + x178 - x186*(x172 - x181 - x184 + x185 + x49*(-x114 + x165*(x183 + 1) + x182 + 1))) + 3*Ss - x121*x161*(-x162 - x168 + x170 + x172 + x49*(-x163 + x164 + x165*(x104 + 1) + 2)) - x138*xi_yy + x139*xi_y + x141*x148 + x141*x153 + x141*x43 - x143*x144 + x147*xi_y - x149*x150 + x150*x151 - x156*x159 + x159*x160 + x42*xi_yy + x44*xi_yy)) + x66*(Fxt*x137 + u*(B1*x190 + B1t*x188 + B2*x190 + B2t*x188 + B2t*x192 + B2t*x194 + Fx*(B1t*x177 + B2t*x177 + 2*St*(x173 + x85*(x174 + x175)) + x59*(x185 + x200 + x49*(x182 + x199))) + 3*Sb + x121*x196*(x162 + x168 + x170 + x198 + x49*(x164 + x197)) + x127*x187 - x138*xi_xx + x145*x191 + x148*x189 + x149*x193 + x151*x193 + x153*x189 + x154*x191 + x156*x195 + x160*x195 + x189*x43 + x42*xi_xx + x44*xi_xx)) - x69*(Fxh*x137 + Fyt*x137 + u*(B2*x142*x203 + B2*x210*x83 + B2h*x188 + B2h*x192 + B2h*x194 + B2t*x140 + B2t*x146 + B2t*x155 + Fx*(x176*(x173 + x175*x85) + x178 + x186*(x199*x49 + x200) + x213*(-3*x14*(x165 - x166*x175 - x46 + x99 + 1) + x197*x220)) + 6*Sc - 6*Sp*xi_xy - x123*x59*x6 + x127*x203 + x128*xi_xy + x129*xi_xy - x134*x204 + x139*xi_x + x147*xi_x - x157*x204 + x158*x206*x218 + x167*x211 + x179*x204 + x180*x204 + 90*x201*x8 - x202*x99 - x202 + 12*x205*x208 + 18*x205*x3 + x207*x212 + 4*x207 - 12*x209 + x21*x210 - x213*x214*x6*xi_x + x217 + x218*x219*x86)))/3 + 4*x19*(x52^2*x57 - x52*x69*(2*x61 + 2*x62) + x63^2*x66)/3 + x2*(48*x10*x33*x35 + 324*x18 + 1944*x20 + 1296*x22 + 1296*x25 + x26*x28 + 324*x29*(x30 + x31 + 6) + 432*x32*x7 + 216*x36*x37*x39)/27 - x2*(-S*x70 - 648*Sd*x39*(x13 - 1) + 324*Sp*x95*(x9^2 + x94) - u^12*x29*x77 + 108*x10*x6*(x100*x41 + x104 - x107 - x111*x41 - x112*x81 + x115 + x46 + 3*x90 - x94*(x103 + x113 + x8*(x114 - 7)) - x98 + x99*(x102 + x108*x109 + x110 - x8 + 2*x96)) - x106*x35*(-7*S*x100 + x105 - x14*x97 + x5*x99 + 5*x8 + x94*(x101 + x103) + x96 + x98) + x122*x37*(Sp*u*x5*(x110 + x31 + x5*x94 - 2*x8 + x96) - x108*x121*(x120 + x5 + x8*(7*x4 - 5)) + x116*x74 - 5*x117 - x118*x119 - 11*x21 + 6*x72 - 7*x76 + 22*x79 - 2*x82 - 7*x84 + 8*x87 + x97) - x23*x73*x88 - x24*x77 + x28*x36*x48 - 5508*x29*x90 + 1944*x29 + x70*x78 + 1296*x72 - x73*x74 - x75*x91*x92 - 1620*x76 - 3888*x79 - 3240*x82 - 972*x84 - 5184*x87)/9 - x2*(-8748*B2p*x235*x29 + 32*UU(phiouter,potential)*phi0^18*u^16*x252*(phi*x11*x33 - x5^5) + 48*UU(phiouter,potential)*phi0^16*x17*x261*(-phi*x26*x263 + x262 + 9*x265) + phi0^20*u^18*x250*x251*x26 + 4*phi0^14*x5*x88*(x266*(-432*UU(phiouter,potential)*x265 - x251*x262 + x268*x269 + x273*(48*x271 - 1)) - x267*x33*(phip - 3*x266)) + phi0^12*x290*(u*(162*x250*x284 - 864*x252*x269*x36 + x26*(2*phid*(x285 + 3) + x26*x286) + x261*x287*x33*(18*x271 - 1) + x273*(72*x279 - x5*(32*x288 + x289))) - x282*x33*(x280 + x281)) - phi0^10*x106*x85*(phip*x293*(x281 + x291) + u*(-phi*x293*(27*x279 - x5*(x289 + x294*x295)) - 54*x264*x5*(x251*x292 - 1) + x268*x284 + x36*(x280*x48 - x5*(x263*x296 + x297)))) + 324*u*x37*(x282*x283 + x95*(-x267*x298 + x272*x38*(x47 - 5) + x5*(UU(phiouter,potential)*x120 - UU(phiouter,potential)*x299 - x286*x300 + x286*x301 + x286*x303 + x286*x4 - x286 - x299 + x30 + x305 + 12*x31 + 4*x8*(x297 + x302) - 4))) - 1944*x117 + 486*x18*(-x239 + x247*(-B1p + x257) + x254 - x256 - 4) - 2916*x20*x259 - x22*x260 - x221 + x222*x223 - x222*x224*xi + 5832*x223*x234 + 34992*x223*x29*x89 - 5832*x224*x71*x89 - x225*x236 - x225*x80 + x226*x227 + x226*x237 + x226*x238 + x227*x244 - x228*x229 - x228*x232 - x228*x242 - x229*x245 + x230*x233 + x230*x86 + x231*x233 + x231*x86 - x232*x245 + x233*x246 + x237*x244 + x238*x244 - x239*x240 - x240*x248 - x241*x92 - x242*x245 + x246*x86 - x249*x92 - x25*x260 - x27*x276*x80*(108*phip*x38*(4*x279 + x281) + u*(-36*phi*x38*(36*x279 - x5*(x285 + 16*x288 + 15)) - 216*x307*(UU(phiouter,potential)*x11*x5 - 1) + x5*(216*x308 - x5*(-144*UU(phiouter,potential)*x96 - x300*x312 + x301*x312 + x303*x312 + x305 - x309 + 72*x31 + 144*x310 + x312*x4 - x312 + 44*x4 + 24*x8*(x297 + 6*x302) - 44)))) - x287*x311*(phip*x306*(phid*x11 + x105) + u*(-phi*x306*(-x278*(x286*x295 + x289) + x291) + x278*(36*x308 - x5*(-UU(phiouter,potential)*x309 + x13*(x212 + 8*x302 - 3*x31 - 6) - x294*x300 + x294*x301 + x294*x303 + x294*x4 - x294 + x305 + 27*x31 + 24*x310 + 14*x4 - 9*x96 - 14)) + 27*x307)) - 324*x32*x6*(x107 + x116*x241 + x116*x249 + x13*(x183 + x31 + 2) - 2*x239 - 2*x248 + x255*x276 - x274*x275 + x274*x89 - x275*x277 + x277*x89 + x278 - 9*x31) - 7776*x72)/81

    nothing
end

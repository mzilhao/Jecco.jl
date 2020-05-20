
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

function S_eq_coeff!(ABCS::Vector, vars, ::Inner)
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

function Fxy_eq_coeff!(AA::Matrix, BB::Matrix, CC::Matrix, SS::Vector, vars, ::Inner)
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


function Sd_eq_coeff!(ABCS::Vector, vars, ::Inner)
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

    x0 = u^2
    x1 = u^4
    x2 = B1*x1
    x3 = exp(x2)
    x4 = S*x1
    x5 = 3*u
    x6 = x5*xi
    x7 = u*xi
    x8 = x7 - 1
    x9 = phi0^2
    x10 = x0*x9
    x11 = x10*x8 + x6 + 3
    x12 = x11 + 3*x4
    x13 = x12^3*x3
    x14 = 8*u/9
    x15 = x12^2
    x16 = u^3
    x17 = Sp*x16
    x18 = 2*x7
    x19 = x18 - 1
    x20 = x10*x19
    x21 = x20 + 9*x4 - 3
    x22 = xi^3
    x23 = S*u
    x24 = x7 + 1
    x25 = x24^3
    x26 = u^5
    x27 = u^9
    x28 = S^3*x27
    x29 = phi0^6
    x30 = 4*x16*x29
    x31 = x8^3
    x32 = x24*x8^2
    x33 = phi0^4
    x34 = 108*x9
    x35 = S^2
    x36 = u^7
    x37 = xi^2
    x38 = x0*x37
    x39 = x24^2
    x40 = x39*x8
    x41 = x16*x22
    x42 = B2*x1
    x43 = exp(x42)
    x44 = G*x1
    x45 = cosh(x44)
    x46 = x43*x45
    x47 = 3*Sh
    x48 = x9*xi
    x49 = u*(9*x23*xi_y + x47 + 2*x48*xi_y)
    x50 = Fy*x21
    x51 = x49 + x50
    x52 = exp(2*x2 + x42)
    x53 = x45*x52
    x54 = 3*St
    x55 = S*u*xi_x
    x56 = x9*xi*xi_x
    x57 = x54 + 9*x55 + 2*x56
    x58 = u*x57
    x59 = Fx*x21
    x60 = x58 + x59
    x61 = sinh(x44)
    x62 = exp(x2 + x42)
    x63 = x61*x62
    x64 = 2*u
    x65 = UU(phi, potential)*x0
    x66 = x16*x37
    x67 = UU(phi, potential)*x1
    x68 = xi^4
    x69 = x26*x68
    x70 = 3*x20 + 2*x67 - 6
    x71 = x16*xi
    x72 = 8*x26
    x73 = phi0^8
    x74 = 15*x7 - 15
    x75 = UU(phi, potential)*x26*xi
    x76 = x26*x35
    x77 = u*x37
    x78 = u^6
    x79 = 16*x78
    x80 = 8*UU(phi, potential)
    x81 = x36*x68
    x82 = xi^5
    x83 = 4*x26
    x84 = x16*x9
    x85 = x1*x9
    x86 = x26*x37
    x87 = u^8
    x88 = 12*Sp
    x89 = -12*x7
    x90 = x1*x68
    x91 = -x18
    x92 = 4*x1
    x93 = 4*G*u
    x94 = Fy*x0
    x95 = Gh + x93*(x94 + xi_y)
    x96 = Fx*x0
    x97 = 18*x4
    x98 = x10*(4*x7 - 1) + x97 - 3
    x99 = 36*S*x1
    x100 = x16*x88
    x101 = -x100
    x102 = 4*x0*x19*x9
    x103 = x101 + x102 + x99 - 12
    x104 = 36*S*u
    x105 = 72*Sh*u
    x106 = 12*Sh
    x107 = B1h*x1
    x108 = B2h*x1
    x109 = Fyp*x0
    x110 = 8*x9*xi
    x111 = xi_y^2
    x112 = 8*x9
    x113 = B1*x26*xi_y
    x114 = 36*S*x26*xi_y
    x115 = 48*B2*Sh*x26
    x116 = 144*S*x0
    x117 = 144*B1*S*x78
    x118 = 144*B2*S*x78
    x119 = 24*u*x9*xi
    x120 = 8*x9*xi*xi_y
    x121 = x0*x9*xi
    x122 = 32*B1*x26*x9*xi
    x123 = 32*B2*x26*x9*xi
    x124 = 2*x0
    x125 = Fy^2
    x126 = 24*x42
    x127 = 117*x4
    x128 = 72*B2*S*x87
    x129 = x2*(72*S*x1 - 24)
    x130 = 22*x7
    x131 = 8*x2
    x132 = x131*x19
    x133 = 8*x42
    x134 = 24*xi_y
    x135 = 78*x16
    x136 = Sh*x135
    x137 = 48*B1*x36
    x138 = 48*B1*x1
    x139 = 48*B2*x36
    x140 = Sh*x139
    x141 = 48*B2*x1
    x142 = x141*xi_y
    x143 = 378*S*x1
    x144 = x143*xi_y
    x145 = B1*xi_y
    x146 = 288*S*x87
    x147 = 288*B2*S*x87
    x148 = 8*x0*x9
    x149 = 16*x78*x9*xi_y
    x150 = x16*xi_y
    x151 = 68*x150*x9*xi
    x152 = 64*x36*x9*xi
    x153 = B1h*x16
    x154 = x102 + x99 - 12
    x155 = 9*S*x26
    x156 = x1*x9*xi
    x157 = -x155 - 2*x156 + x5 + x84
    x158 = B2h*x16
    x159 = 72*St*u
    x160 = 12*St
    x161 = B1t*x1
    x162 = B2t*x1
    x163 = 3*St*x0
    x164 = xi_x^2
    x165 = B1*x26*xi_x
    x166 = 36*S*x26*xi_x
    x167 = 48*B2*St*x26
    x168 = 9*S
    x169 = 8*x9*xi*xi_x
    x170 = Fx^2
    x171 = 24*xi_x
    x172 = St*x135
    x173 = St*x139
    x174 = x141*xi_x
    x175 = x143*xi_x
    x176 = B1*xi_x
    x177 = x9*xi_x
    x178 = 16*B2*x78
    x179 = x177*x178
    x180 = 68*x16*x56
    x181 = 64*x36*x9*xi*xi_x
    x182 = B1t*x16
    x183 = B2t*x16
    x184 = Fx*x134
    x185 = Fy*x171
    x186 = Fx*Fyp
    x187 = 12*Fx
    x188 = 16*x9*xi
    x189 = 16*xi_y
    x190 = B2*x78
    x191 = B2h*xi_x
    x192 = 36*S*x26
    x193 = B2t*Fy
    x194 = B2t*xi_y
    x195 = Fx*Fy*x78
    x196 = Fx*Fy*x26
    x197 = 288*S*xi_y
    x198 = x0*xi_x
    x199 = u^10
    x200 = B2*Fx*Fy*x199
    x201 = B2*Fx*xi_y
    x202 = B2*Fy*xi_x
    x203 = B2*x78*xi_x
    x204 = 4*B2h*x26
    x205 = 4*Fy
    x206 = 36*Fx*Fy
    x207 = B2*Fx
    x208 = 8*x78*x9*xi
    x209 = 2*Fyp
    x210 = x209*xi_x
    x211 = 64*B2*x27
    x212 = B2*Fx*x36
    x213 = B2*Fy*x36
    x214 = 64*x9*xi*xi_x
    x215 = Fxp*u
    x216 = 8*Gh*u
    x217 = Gh*x1
    x218 = 2*B1h
    x219 = 10*Fy*x16
    x220 = 4*G*x16
    x221 = Gp*x0
    x222 = 8*Gh
    x223 = B1*Fy*x36
    x224 = 8*Fy
    x225 = B1h*G*x36
    x226 = 4*Gh
    x227 = 4*B2*Gh*x26
    x228 = B2h*G*x36
    x229 = 4*B2h*G*x26
    x230 = Fy*xi_y
    x231 = 56*G*x1
    x232 = 2*Fy*Gp
    x233 = 36*G*x78
    x234 = 2*Gp*x26
    x235 = 20*G*x0
    x236 = 64*G
    x237 = 32*B2*Fy*G*x87
    x238 = 32*B1*G*x199
    x239 = 32*B1*G*x78
    x240 = 16*B2*G*x199
    x241 = 16*B2*G*x78
    x242 = 8*Gt*u
    x243 = Gt*x1
    x244 = 2*B1t
    x245 = 10*Fx*x16
    x246 = Gt*x0
    x247 = 8*Gt
    x248 = B1*Fx*x36
    x249 = 8*Fx
    x250 = B1t*G*x36
    x251 = B1t*G
    x252 = 4*Gt
    x253 = 4*B2*Gt*x26
    x254 = 4*Fx
    x255 = B2t*G*x36
    x256 = 4*B2t*G*x26
    x257 = 4*Fxp
    x258 = Fx*G*x26
    x259 = 56*Fx*G*x1
    x260 = 2*Fx*Gp
    x261 = x16*xi_x
    x262 = x257*xi_x
    x263 = G*x16
    x264 = 32*B2*Fx*G*x87
    x265 = 2*Fxp
    x266 = 2*B2p*u
    x267 = Fx*u
    x268 = Fy*u
    x269 = 8*B2*x0
    x270 = 8*x0
    x271 = Fx*Fy
    x272 = 16*x16
    x273 = 4*Fyp
    x274 = 2*x26
    x275 = 16*Fx*x87
    x276 = 16*Fy*x87
    x277 = 72*B2*x36
    x278 = 56*x26
    x279 = 4*Fyp*xi_x
    x280 = B2p*x1
    x281 = Fy*xi_x
    x282 = 8*G*Gh
    x283 = Fx*x87
    x284 = 8*G*Gt
    x285 = Fy*x87
    x286 = x78*xi_x
    x287 = x78*xi_y
    x288 = B2^2
    x289 = 64*x288
    x290 = u^11
    x291 = Fx*Fy*x290
    x292 = Fx*xi_y
    x293 = 64*x27*x288
    x294 = x36*xi_x*xi_y
    x295 = G^2
    x296 = 32*x295
    x297 = 4*Fx*Fy*x26
    x298 = 32*x27*x295
    x299 = phih*x1*x33
    x300 = phit*x1*x33
    x301 = 24*phi*x33*x36
    x302 = 12*xi_y
    x303 = Fx*phi*x26*x33
    x304 = phi*phih*x29*x78
    x305 = 12*xi_x
    x306 = Fy*phi*x26*x33
    x307 = x1*x9*xi*xi_x
    x308 = phi*phih*x1*x29
    x309 = phi*phit*x1*x29
    x310 = phih*x16*x33*xi
    x311 = phit*x16*x33*xi
    x312 = 48*phi*x33*xi
    x313 = xi_x*xi_y
    x314 = 48*phi*x1*x33*xi
    x315 = phi^2
    x316 = 16*x36*x37*x9
    x317 = 36*x29*x315*x36
    x318 = x26*x37*x9*xi_x
    x319 = 36*x26*x29*x315
    x320 = B2p*x16
    x321 = -x320 + 4*x42 + 2
    x322 = 2*B1p*u
    x323 = 8*B1*x0
    x324 = 20*Fy
    x325 = 16*x0*xi_y
    x326 = B1*B1h
    x327 = 16*x78*xi_y
    x328 = 8*B1*B2h
    x329 = B1*Fy
    x330 = 8*Fyp*x78
    x331 = 112*x26*xi_y
    x332 = Fyp*xi_y
    x333 = 8*B1h*B2
    x334 = 4*Fy*xi_y
    x335 = B1p*x1
    x336 = 32*B2*B2h
    x337 = B2*Fy
    x338 = G*Gh
    x339 = 72*B1*x36
    x340 = 40*B1*x16
    x341 = 4*B1p*x78
    x342 = 40*B2*x16
    x343 = 4*B2p*x78
    x344 = 32*B1*B2*x290
    x345 = 32*B1*B2*x36
    x346 = B1^2
    x347 = 64*x27*x346
    x348 = 128*x27*x288
    x349 = 64*x27*x295
    x350 = 32*x290*x346
    x351 = 32*x346*x36
    x352 = 64*x288*x290
    x353 = 64*x288*x36
    x354 = 32*x290*x295
    x355 = 4*x26*x9
    x356 = 32*x295*x36
    x357 = 16*x26
    x358 = 16*x78*x9*xi
    x359 = 96*phi*x33*xi
    x360 = 72*x29*x315*x36
    x361 = 48*phi*x33*x87*xi
    x362 = 36*x27*x29*x315
    x363 = 16*x16*x37*x9
    x364 = B1p*x16
    x365 = 20*Fx
    x366 = 16*x0*xi_x
    x367 = B1*B1t
    x368 = 16*x78*xi_x
    x369 = 8*B1*B2t
    x370 = B1*Fx
    x371 = 8*Fxp*x78
    x372 = 112*x26*xi_x
    x373 = Fxp*xi_x
    x374 = 4*Fx*xi_x
    x375 = 8*B1t*B2
    x376 = 32*B2*B2t
    x377 = G*Gt
    x378 = Fx*xi_x

    ABCS[1] = 0

    ABCS[2] = -4*x0*x13/9

    ABCS[3] = -x13*x14 - x14*x15*x3*(-3*x17 + x21)

    ABCS[4] = u*x15*(x45*(-2*Fyph*x43 + u*x43*(-B1*Fy*x211*xi_y + B1h^2*x274 - B1h*x325 - B1s*x64 + B2h^2*x83 - B2h*x218*x26 + B2h*x325 + B2s*x64 - Fy*phih*x33*x357*xi + 24*Fy*x304 - Fy*x359*x78*xi_y + 32*Fy*x86*x9*xi_y + 2*Fyh*(-4*x2 + x321 + x364) + Fyp^2*u - 8*Fyp*x94 + Gh^2*x274 + phih^2*x30 - 16*x1*x230*x9*xi - x107*x324 + x108*x324 - x111*x314 + x111*x319 - x111*x340 + x111*x342 - x111*x345 + x111*x351 + x111*x353 + x111*x356 + x111*x363 + x125*x272 + x125*x277 + x125*x301 + x125*x316 - x125*x339 + x125*x341 - x125*x343 - x125*x344 + x125*x350 + x125*x352 + x125*x354 + x125*x355 - x125*x358 - x125*x361 + x125*x362 + x131*x332 - x133*x332 + x134*x306 + x134*x308 + x153*x209 - x158*x209 - x189*x310 + x224*x299 + x230*x347 + x230*x348 + x230*x349 + x230*x360 - x266*xi_yy + x268*x302 + x269*xi_yy - x273*xi_y + x276*x326 + x276*x338 - x280*x334 - x285*x328 - x285*x333 + x285*x336 - x287*x328 - x287*x333 + x287*x336 + x322*xi_yy - x323*xi_yy + x326*x327 + x327*x338 + x329*x330 - x329*x331 - x330*x337 + x331*x337 + x334*x335) - x124*x62*(-4*B1h*G*x26*xi_x - B1h*x243 + B1t*x217 + B2h*x243 + B2t*x217 + Fxh*x220 - Fxh*x221 - Fxp*x0*x95 - Fyp*x246 + Fyt*x220 - Fyt*x221 + 72*G*x195 + 40*G*x198*xi_y + 32*G*x200 + 32*G*x203*xi_y + 2*Gc + 4*Gh*x176*x26 + Gh*x245 - Gp*x297 - 2*Gp*xi_xy - 4*Gt*x145*x26 + Gt*x219 + (8*u)*(G*xi_xy) - x150*x260 + x205*x250 + x205*x255 + x212*x226 + x213*x252 + x216*xi_x - x223*x252 - x225*x254 + x226*x248 + x227*xi_x + x228*x254 + x229*xi_x + x231*x281 - x232*x261 + x237*xi_x + x242*xi_y + 4*x251*x26*xi_y + x253*xi_y + x256*xi_y - x258*x273 + x259*xi_y - x263*x279 + x264*xi_y) + x52*(-2*Fxpt + u*(B1*Fx*x211*xi_x + B1b*x64 + B1t^2*x274 + B1t*x366 + B2b*x64 + B2t^2*x83 + B2t*x244*x26 + B2t*x366 + 24*Fx*phi*phit*x29*x78 - Fx*phit*x33*x357*xi - 16*Fx*x307 + 32*Fx*x318 - Fx*x359*x78*xi_x + Fxp^2*u - 8*Fxp*x96 + Fxt*(x131 + x133 - 2*x320 - 2*x364 + 4) + Gt^2*x274 + phit^2*x30 - x131*x373 - x133*x373 + x161*x365 + x162*x365 - x164*x314 + x164*x319 + x164*x340 + x164*x342 + x164*x345 + x164*x351 + x164*x353 + x164*x356 + x164*x363 + x170*x272 + x170*x277 + x170*x301 + x170*x316 + x170*x339 - x170*x341 - x170*x343 + x170*x344 + x170*x350 + x170*x352 + x170*x354 + x170*x355 - x170*x358 - x170*x361 + x170*x362 + x171*x303 + x171*x309 - x182*x265 - x183*x265 - x207*x371 + x207*x372 + x249*x300 - x262 - x266*xi_xx + x267*x305 + x269*xi_xx + x275*x367 + x275*x377 - x280*x374 + x283*x369 + x283*x375 + x283*x376 + x286*x369 + x286*x375 + x286*x376 - 16*x311*xi_x - x322*xi_xx + x323*xi_xx - x335*x374 + x347*x378 + x348*x378 + x349*x378 + x360*x378 + x367*x368 + x368*x377 - x370*x371 + x370*x372))) + x61*(x124*(x43*(-B1*Fy*x236*x87*xi_y - B1h*G*x72*xi_y + B2h*x217 + Fyh*x220 - Fyh*x221 + Gh*x219 - Gp*xi_yy + Gs - x109*x95 + x111*x235 - x111*x239 + x111*x241 - x113*x222 + x125*x233 - x125*x234 - x125*x238 + x125*x240 - x150*x232 + x205*x228 + x213*x226 + x216*xi_y - x217*x218 - x222*x223 - x224*x225 + x227*xi_y + x229*xi_y + x230*x231 + x237*xi_y + x93*xi_yy) + x52*(B1*Fx*x236*x87*xi_x + B2t*x243 - Fxp*x246 + Fxt*x0*(-Gp + x93) + Gb - Gp*xi_xx + Gt*x245 + x164*x235 + x164*x239 + x164*x241 + x165*x247 + x170*x233 - x170*x234 + x170*x238 + x170*x240 + x212*x252 + x242*xi_x + x243*x244 + x247*x248 + x249*x250 + 8*x251*x26*xi_x + x253*xi_x + x254*x255 + x256*xi_x - x257*x258 + x259*xi_x - x260*x261 - x262*x263 + x264*xi_x + x93*xi_xx)) + x62*(2*Fxph + 2*Fypt - x64*(B2*B2h*x275 + B2*B2t*x276 + 40*B2*x16*xi_x*xi_y + B2c*x64 - 2*B2p*x1*x281 - 4*B2p*x195 + B2t*x204 - Fx*Fy*x312*x87 - 8*Fx*phih*x26*x33*xi + 10*Fx*x108 - 2*Fx*x280*xi_y - Fx*x312*x78*xi_y + 16*Fx*x86*x9*xi_y + Fxh*x321 - 4*Fxp*Fy*x190 - Fxp*x158 - 4*Fxp*x42*xi_y + 12*Fy*phi*phit*x29*x78 - 8*Fy*phit*x26*x33*xi + 10*Fy*x162 - Fy*x312*x78*xi_x + 16*Fy*x318 - Fyp*x183 - 4*Fyp*x207*x78 + Fyp*x215 + Fyt*x321 + Gh*Gt*x274 + phih*phit*x30 - 8*x1*x292*x9*xi + 16*x16*x37*x9*xi_x*xi_y + x178*x191 + x178*x194 + x187*x304 - x188*x195 + x191*x270 + x194*x270 + x201*x278 + x202*x278 + x205*x300 + x206*x27*x29*x315 - x210 - x224*x307 + x254*x299 - x257*x94 - x265*xi_y - x266*xi_xy + 6*x267*xi_y + 6*x268*xi_x + x269*xi_xy + x271*x272 + x271*x277 + x271*x301 + x271*x316 - x273*x96 - x279*x42 + x281*x293 + x281*x298 + x281*x317 + x282*x283 + x282*x286 + x284*x285 + x284*x287 + x289*x291 + x289*x294 + x291*x296 + x292*x293 + x292*x298 + x292*x317 + x294*x296 + x297*x9 + x302*x303 + x302*x309 + x305*x306 + x305*x308 - 8*x310*xi_x - 8*x311*xi_y - x313*x314 + x313*x319))))/9 + 2*x0*x12*(-x45*x62*x92*(Fx*x93*(x16*x47 + x94*(2*x0*x19*x9 + x97 - 6) + x98*xi_y) + Fy*x93*(x16*x54 + x98*xi_x) + 4*G*x0*(x47*xi_x + xi_y*(x54 + 18*x55 + 4*x56)) + Gh*x58 + Gh*x59 + Gt*x49 + Gt*x50) + x46*(Fyh*x103 + u*(-B1h*x114 + B2h*x114 + Fy*(B1*x149 - B2*x149 + B2*x152*xi_y + Fyp*x157 - Sh*x137 - x134*x17 - x134 + x136 + x138*xi_y + x140 - x142 + x144 - x145*x146 - x145*x152 + x147*xi_y - x148*xi_y + x151 - x153*x154 - x158*(-36*S*x1 + 4*x9*(x0 - 2*x71) + 12)) - 9*Fyp*S*x16*xi_y - 2*Fyp*x121*xi_y - 48*Sh*x113 + 12*Ss + x104*xi_yy + x105*xi_y - x106*x107 + x106*x108 - x107*x120 + x108*x120 - x109*x47 + x110*xi_yy + x111*x112 + x111*x116 - x111*x117 + x111*x118 + x111*x119 - x111*x122 + x111*x123 + x115*xi_y - x124*x125*(x10*(-x130 + x132 + x133*(x91 + 1) + 9) + x100 + x126 - x127 - x128 + x129 + 15) - x88*xi_yy)) + x53*(Fxt*x103 + u*(B1t*x166 + B2t*x166 + Fx*(B1*x181 - B1*x79*x9*xi_x + B2*x181 + Fxp*x157 + St*x137 - x138*xi_x + x146*x176 + x147*xi_x - x148*xi_x + x154*x182 + x154*x183 - x17*x171 - x171 + x172 + x173 - x174 + x175 - x179 + x180) - 2*Fxp*x0*x56 - Fxp*x16*x168*xi_x - Fxp*x163 + 12*Sb + 48*St*x165 + x104*xi_xx + x110*xi_xx + x112*x164 + x116*x164 + x117*x164 + x118*x164 + x119*x164 + x122*x164 + x123*x164 + x124*x170*(x10*(x130 + x132 + x133*x19 - 9) + x101 - x126 + x127 + x128 + x129 - 15) + x159*xi_x + x160*x161 + x160*x162 + x161*x169 + x162*x169 + x167*xi_x - x88*xi_xx)) + x61*x92*(x43*x51*x95 + x52*x60*(Gt + x93*(x96 + xi_x))) - x63*(Fxh*x103 + Fyt*x103 + u*(-32*B2*Fx*Fy*x87*x9 + B2*x214*x26*xi_y + 36*B2h*Fx*S*x36 + B2h*Fx*x208 - B2t*x205*x26*x9 - 96*Fx*Fy*x190 + Fx*Fy*x211*x9*xi - 2*Fx*Fyp*x156 - Fx*x0*x112*xi_y + Fx*x136 + Fx*x140 - Fx*x142 + Fx*x144 + Fx*x151 - Fx*x204*x9 - 8*Fy*x0*x177 + Fy*x172 + Fy*x173 - Fy*x174 + Fy*x175 - Fy*x179 + Fy*x180 - 12*Fy*x183 - 60*Fy*x96 - Fyp*x16*x168*xi_x - Fyp*x163 + 36*S*x193*x36 + 468*S*x195 + 288*S*x200 + 24*Sc - 48*Sp*x196 - 24*Sp*xi_xy + 48*u*x177*xi*xi_y + x105*xi_x + x106*x162 + x108*x160 + x108*x169 + x115*xi_x + x120*x162 - x121*x210 + x146*x201 + x146*x202 - x149*x207 - x155*x186 - x158*x187 + x159*xi_y + x167*xi_y - x17*x184 - x17*x185 + x177*x189 - x184 - x185 + x186*x5 + x186*x84 + x188*xi_xy + x191*x192 + x192*x194 + x193*x208 + 88*x196*x9*xi + x197*x198 + x197*x203 - x206*x85 + 64*x212*x9*xi*xi_y + x213*x214 - x215*x51 + 72*x23*xi_xy)))/9 + 4*x26*(-x46*x51^2 + x51*x63*(x57*x64 + 2*x59) - x53*x60^2)/9 + x3*(36*u*x32*x33*(x24 + x4) + 108*x22*(x7 + 4) + 324*x23*x25 + x24*x30*x31 + 324*x26*(S*x7 + S)^2 + x28*(108*u*xi + 108) + x34*(2*S*x16*x40 + x35*x36*(x38 - 1) + xi*(2*x38 + x41 - 2)))/9 + x3*(324*S^4*u^13*x70 + 3888*UU(phi, potential)*x66 + 648*UU(phi, potential)*x69 + 12*phi0^10*x19*x36*x8^4 + 648*u*(UU(phi, potential) - 3*x68) + 48*x11^3*x23*x70 + 216*x11^2*x70*x76 + 432*x11*x28*x70 + 24*x16*x29*x32*(54*x38 - 4*x67 + x74 + 4*x75) + 432*x16*x33*x40*(UU(phi, potential)*x71 + 6*x37 - x65) + 2592*x22*x67 - 7776*x22 + x31*x72*x73*(36*x38 - x67 + x74 + x75) + x34*(UU(phi, potential)*x22*x79 + 24*x0*x22 - 2*x1*(x80*xi - 9*x82) + x16*(39*x68 - x80) + 18*x77 + x80*x81 + 30*xi) + 2592*x65*xi)/81 - x3*(-324*u*x68 - 216*x1*x22*x33 + 48*x1*x29*xi + 432*x10*x22 + x15*x88*(x10 - 3*x39) - 24*x16*x29 - 96*x22*x29*x78 - 8*x22*x73*x87 - 1296*x22 - x23*(-108*x0*x9*(x38 + 20*x41 - 16*x7 + 10*x90 - 7) - 36*x1*x33*(-22*x38 + 2*x41 + 7*x90 + x91 + 11) - 324*x25*(x6 + 1) + 12*x29*x78*(7*x38 + x89 + 5)) - x28*(324*x0*x9 - 972*x39) - 12*x29*x81 + 24*x29*x82*x87 + 60*x29*x86 - 252*x33*x66 + 180*x33*x69 + 144*x33*x78*x82 + 20*x36*x37*x73 + 432*x48 + 540*x68*x84 - x73*x79*xi + x73*x83 - x76*(-108*x0*x9*(9*x38 + 8*x41 + x89 - 12) + 36*x1*x33*(8*x7 - 7) - 324*x39*(6*x7 + 5)) + 324*x77*x9 + 216*x82*x85)/27
    nothing
end


function B2d_eq_coeff!(ABCS::Vector, vars, ::Inner)
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

function B1dGd_eq_coeff!(AA::Matrix, BB::Matrix, CC::Matrix, SS::Vector, vars, ::Inner)
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

function phid_eq_coeff!(ABCS::Vector, vars, ::Inner)
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

    x0 = u.^2
    x1 = u.^4
    x2 = B1.*x1
    x3 = exp(x2)
    x4 = phi0.^3
    x5 = u.*xi
    x6 = x5 - 1
    x7 = x0.*x6
    x8 = S.*x1
    x9 = x5 + 1
    x10 = x8 + x9
    x11 = 3*x10
    x12 = 3*x5
    x13 = phi0.^2
    x14 = x12 + x13 .* x7 + 3*x8 + 3
    x15 = 12*x5
    x16 = u.^3
    x17 = Sp.*x16
    x18 = x0.*x13
    x19 = u.*x4
    x20 = S.^3 .*u.^11
    x21 = 27*xi
    x22 = 6*x5
    x23 = u.^7
    x24 = S.^2 .*x23
    x25 = x10.^2
    x26 = 4*x5
    x27 = x0.*xi.^2
    x28 = 3*x27
    x29 = S.*x16
    x30 = 2*x5
    x31 = x30 - 1
    x32 = phi0.^6
    x33 = u.^5
    x34 = x6.^2
    x35 = x33 .* x34
    x36 = x1.*xi
    x37 = phi0.^4
    x38 = 3*x37
    x39 = x16.*x6
    x40 = 2*x27
    x41 = x40 + 1
    x42 = u.*x13
    x43 = 9*x10
    x44 = 0.222222222222222*phi0.*x3
    x45 = 4*u
    x46 = x23 .* x6.^3
    x47 = 4*x39
    x48 = phi.*x10
    x49 = 9*x48
    x50 = -x34
    x51 = x9.^2
    x52 = x16.*xi.^3
    x53 = u.^6
    x54 = 9*S
    x55 = u.*xi_y
    x56 = 2*xi
    x57 = x56.*xi_y
    x58 = x18 .* x31 + 9*x8 - 3
    x59 = Fy.*x58 + u.*(3*Sh + x13 .* x57 + x54.*x55)
    x60 = G.*x1
    x61 = cosh(x60)
    x62 = x0.*x56
    x63 = u - x62
    x64 = Fy.*x0
    x65 = x64 + xi_y
    x66 = 3*phi
    x67 = u.*x66
    x68 = phih + x65 .* x67
    x69 = Fy.*x63 + x13 .* x68 - x57
    x70 = B2 .*x1
    x71 = exp(x70)
    x72 = x69 .* x71
    x73 = u.*xi_x
    x74 = x56.*xi_x
    x75 = u.*(3*St + x13 .* x74 + x54.*x73)
    x76 = Fx.*x58
    x77 = Fx.*x0
    x78 = x77 + xi_x
    x79 = Fx.*x63 + x13 .* (phit + x67 .* x78) - x74
    x80 = exp(2*x2 + x70)
    x81 = x79 .* x80
    x82 = exp(x2 + x70)
    x83 = sinh(x60)
    x84 = x82 .*x83
    x85 = 1.33333333333333*phi0
    x86 = 54*x31
    x87 = 2*x1
    x88 = 6*x1
    x89 = Sd.*x6
    x90 = -x12
    x91 = Sd.*x87
    x92 = x90 - 3
    x93 = x27 + x52 + x6.*x91 - 2*x8 + x92
    x94 = x51 + x91
    x95 = phip.*x11
    x96 = 2*x52
    x97 = 4*x1.*x89
    x98 = x12 - x40
    x99 = x8 - x96 - x97 + x98 + 3
    x100 = G.*x45
    x101 = 4*xi
    x102 = xi_x.*xi_y
    x103 = phih.*xi_x
    x104 = phit.*xi_y
    x105 = 6*x73
    x106 = 2*x64
    x107 = 6*u
    x108 = phi.*x107
    x109 = xi_y.^2
    x110 = Fy.*x33
    x111 = Fy.*x16
    x112 = Fy.^2
    x113 = 4*B1
    x114 = Fy.*xi_y
    x115 = x114.*x53
    x116 = Fy.*x53
    x117 = x116.*x56
    x118 = 2*x36
    x119 = x118 .* xi_y
    x120 = 4*B2
    x121 = Fy.*x118
    x122 = 20*xi
    x123 = x111.*xi_y
    x124 = x62 .*xi_y
    x125 = u.^8
    x126 = 4*x125
    x127 = B1.*x112
    x128 = B2 .*x126
    x129 = 14*xi
    x130 = x112 .*x33
    x131 = Fy.*x23
    x132 = 16*xi
    x133 = x132 .*xi_y
    x134 = B2 .*x131
    x135 = u.^9
    x136 = 8*xi
    x137 = x135 .* x136
    x138 = x136.*x33
    x139 = x109 .* x138
    x140 = B2 .*x112
    x141 = phih.*x1
    x142 = 8*x111
    x143 = x16.*x66
    x144 = phip.*x0
    x145 = phih.*x113
    x146 = x33 .* xi_y
    x147 = B1h.*x66
    x148 = phih.*x120
    x149 = B2h.*x66
    x150 = 36*phi
    x151 = x1.*x150
    x152 = 2*phip
    x153 = 24*phi
    x154 = x153 .* x53
    x155 = 12*phi
    x156 = x0.*x155
    x157 = x125 .* x153
    x158 = x114.*x157
    x159 = u.^10
    x160 = x155 .* x159
    x161 = x155 .* x53
    x162 = x109 .* x161
    x163 = Fyp.*x0
    x164 = xi_x.^2
    x165 = Fx.*x33
    x166 = Fx.*x16
    x167 = 2*x77
    x168 = Fx.^2
    x169 = Fx.*x53
    x170 = x169 .* xi_x
    x171 = x169 .* x56
    x172 = x118 .* xi_x
    x173 = Fx.*x118
    x174 = x166.*xi_x
    x175 = Fxp.*xi_x
    x176 = B1.*x168
    x177 = x168 .* x33
    x178 = Fx.*x23
    x179 = B1.*x178
    x180 = x132 .*xi_x
    x181 = B2 .*x178
    x182 = B1.*x164
    x183 = B2 .*x168
    x184 = B2 .*x164
    x185 = phit.*x1
    x186 = 8*x166
    x187 = Fxp.*x0
    x188 = 4*phit
    x189 = x33 .* xi_x
    x190 = x188 .* x189
    x191 = B1t.*x66
    x192 = B2t.*x66
    x193 = x165 .* x66
    x194 = Fx.*x1
    x195 = x150.*x194
    x196 = Fx.*x125
    x197 = x153 .* xi_x
    x198 = B2 .*x196
    x199 = Fy.*x165
    x200 = x166.*xi_y
    x201 = x111.*xi_x
    x202 = Fyp.*xi_x
    x203 = B2 .*Fx.*Fy
    x204 = B2 .*x102
    x205 = x120.*x33
    x206 = Fy.*xi_x

    ABCS[1] = 0

    ABCS[2] = -0.296296296296296*x0.*x3 .* (phi0.*x11 + x4.*x7).^3

    ABCS[3] = -0.148148148148148*x14.^2 .*x19 .* x3 .* (x15 - 9*x17 + x18 .* (10*x5 - 7) + 39*x8 + 3)

    ABCS[4] = x14.*x16.*x85 .* (-x1.*x61.*x82 .*(Gh.*x79 + Gt.*x69 + x100.*(Fx.*u.*(x106.*(1 - x30) - x26.*xi_y + x42 .*(phih + x108 .* x65) + xi_y) + Fy.*u.*(-x26.*xi_x + x42 .*(phi.*x105 + phit) + xi_x) + 6*phi.*x102 .*x42 - x101.*x102 + x103 .* x13 + x104.*x13)) + x1.*x83 .* (x72 .*(Gh + x100.*x65) + x81.*(Gt + x100.*x78)) + x61.*x71.*(B1.*x131.*x133 + B1.*x139 - B1h.*x110 + B1h.*x117 + B1h.*x119 - B2 .*x139 + B2h.*x110 - B2h.*x117 - B2h.*x119 + Fyh.*u - Fyh.*x62 - Fyp.*x111 + Fyp.*x121 + Fyp.*x124 + x106.*xi_y - x109 .* x22 - 2*x109 + x112 .*x128 + x112 .*x88 - x113 .* x115 + x115 .* x120 - x122 .*x123 - x126.*x127 + x127 .* x137 - x129 .* x130 + x13 .* (-B1.*x158 - B1.*x162 - B1h.*x141 + B2 .*x158 + B2 .*x162 + B2h.*x141 + Fyh.*x143 - Fyh.*x144 + phih.*x142 + 6*phih.*x55 - phip.*xi_yy + phis + x109 .* x156 + x112 .*x154 + x114.*x151 - x123 .* x152 - x127 .* x160 - x130.*x152 - x131.*x145 - x131.*x147 + x131.*x148 + x131.*x149 + x140.*x160 - x145 .* x146 - x146.*x147 + x146.*x148 + x146.*x149 - x163 .* x68 + x67 .* xi_yy) - x133 .* x134 - x137 .* x140 - x56.*xi_yy) + x61.*x80.*(B1t.*x165 - B1t.*x171 - B1t.*x172 + B2t.*x165 - B2t.*x171 - B2t.*x172 - Fxp.*x166 + Fxp.*x173 + Fxt.*u - Fxt.*x62 + x113 .* x170 + x120.*x170 - x122 .*x174 + x126.*x176 + x128 .* x168 - x129 .* x177 + x13 .* (B1.*x190 + B1.*x196.*x197 + B1t.*x185 + B2 .*x190 + B2t.*x185 - Fxp.*x193 + Fxt.*x0.*(-phip + x67) + phib - phip.*xi_xx + phit.*x105 + phit.*x186 - phit.*x187 - x143 .* x175 - x152 .*x174 - x152 .*x177 + x154.*x168 + x156.*x164 + x160.*x176 + x160.*x183 + x161.*x182 + x161.*x184 + x178 .* x191 + x178 .* x192 + x179 .* x188 + x181.*x188 + x189 .* x191 + x189 .* x192 + x195 .* xi_x + x197 .* x198 + x67 .* xi_xx) - x137 .* x176 - x137 .* x183 - x138 .* x182 - x138 .* x184 - x164.*x22 - 2*x164 + x167 .* xi_x + x168 .* x88 + x175 .* x62 - x179 .* x180 - x180.*x181 - x56.*xi_xx) + x84.*(-B2h.*x165 + B2h.*x171 + B2h.*x172 - B2t.*x110 + B2t.*x117 + B2t.*x119 - Fxh.*u + Fxh.*x62 + Fxp.*x111 - Fxp.*x121 - Fxp.*x124 - 12*Fy.*x194 - 8*Fy.*x198 + Fyp.*x166 - Fyp.*x173 - Fyt.*u + Fyt.*x62 + x101.*xi_xy + x102 .*x15 + 4*x102 - x106.*xi_x - x116.*x120.*xi_x - x120.*x169 .* xi_y + x122 .*x200 + x122 .*x201 - x13 .* (B2 .*x157 .* x206 + B2h.*x185 + B2t.*x141 + Fxh.*x143 - Fxh.*x144 + 48*Fy.*phi.*x169 - Fyp.*x193 + Fyt.*x143 - Fyt.*x144 + 2*phic + 4*phih.*x181 + phih.*x186 - 4*phip.*x199 + phit.*x142 - phit.*x163 + x0.*x102 .*x153 + x103 .* x107 + x103 .* x205 + x104.*x107 + x104.*x205 + x108 .* xi_xy + x131.*x192 + x134.*x188 - x143 .* x202 + x146.*x192 + x149 .* x178 + x149 .* x189 + x151.*x206 - x152 .*x200 - x152 .*x201 - x152 .*xi_xy + x153 .* x159 .* x203 + x153 .* x198 .* xi_y + x154.*x204 - x187 .* x68 + x195 .* xi_y) + x132 .*x135 .* x203 + x132 .*x204.*x33 + x133 .* x181 + x134.*x180 - x167 .* xi_y + 28*x199 .* xi - x202 .*x62)) + x3 .* (-UUp(phi, potential).*x10.^3 .* x45 - 1.33333333333333*UUp(phi, potential).*x10.*x35 .* x37 - UUp(phi, potential).*x13 .* x25 .* x47 - 0.148148148148148*UUp(phi, potential).*x32 .*x46 + 0.444444444444444*phi.*phi0.^9 .* x46 + 0.444444444444444*phi0.^7 .* x33 .* (x34.*x49 - x6.^4) + phi0.^5 .* x10.*x47 .* (3*x48 + x50) - 12*phi0.*(x20.*x6 + 3*x24.*(x27 - 1) + 3*x29 .* x51.*x6 + xi.*(x40 + x52 - 2)) + 12*x19 .* x25 .* (x48 + x50)) + x44.*(-27*Sp.*x0.*x25 + 81*x20 - x21.*(x5 + 2) + 27*x24.*(x22 + 5) + 27*x29 .* (x26 + x28 + 1) + x31.*x32 .*x35 + x38 .* x39 .* (Sp.*(x16 - x36) + 4*x27 + x6 + x8 .* (7*x5 - 5)) + x42 .*x43 .* (-2*x17 .* x6 + x41 - x5 + x8 .* (8*x5 - 7))) - x44.*(-Sd.*x16.*x25 .* x86 + phi0.^8 .* x34.*x53 .* (phip - x67) + x0.*x38 .* (u.*(-x49 .* x99 + x93 .* (x98 - 1)) + x95 .* x99) + x1.*x32 .*x6.*(phip.*(-x28 + 9*x5 - 3*x52 + 6*x8 - x88 .* x89 + 9) + u.*(9*phi.*x93 + x41 + x90)) - x13 .* x43 .* (u.*x31.*(x40 - x8 + x92 + x96 + x97) - u.*x49 .* x94 + x94.*x95) - x21.*(8*x27 + x30 + 7*x52 + x87 .* xi.^4 - 2) - 27*x23 .* x31.*(S.*x5 + S).^2 - x29 .* x86.*x9.^3) + x53 .* x85 .* (x59 .* x61.*x72 + x61.*x81.*(x75 + x76) + x84.*(-x59 .* x79 + x69 .* (-x75 - x76)))

    nothing
end


function A_eq_coeff!(ABCS::Vector, vars, ::Inner)
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

    x0 = u^4
    x1 = B1*x0
    x2 = exp(x1)
    x3 = u^3
    x4 = phi0^2
    x5 = u*xi
    x6 = x5 - 1
    x7 = x4*x6
    x8 = S*x0
    x9 = x5 + 1
    x10 = x8 + x9
    x11 = u*x10
    x12 = 3*x5
    x13 = 3*x8
    x14 = u^2
    x15 = x14*x7
    x16 = x12 + x13 + x15 + 3
    x17 = 4*x16^4*x2/9
    x18 = S^4*u^14
    x19 = u^6
    x20 = 1944*x19
    x21 = (S*x5 + S)^2
    x22 = S*x14*x9^3
    x23 = S^3
    x24 = u^10
    x25 = x23*x24*x9
    x26 = x6^4
    x27 = phi0^8
    x28 = 4*x19*x27
    x29 = xi^2
    x30 = 4*x5
    x31 = x14*x29
    x32 = x10^3
    x33 = phi0^6
    x34 = x0*x33
    x35 = x6^3
    x36 = phi0^4
    x37 = x6^2
    x38 = x10^2
    x39 = x14*x38
    x40 = B2*x0
    x41 = exp(x40)
    x42 = G*x0
    x43 = cosh(x42)
    x44 = x41*x43
    x45 = 3*Sh
    x46 = 9*S*u
    x47 = 2*xi
    x48 = x47*xi_y
    x49 = u*(x4*x48 + x45 + x46*xi_y)
    x50 = 9*x8
    x51 = 2*x5
    x52 = x51 - 1
    x53 = x14*x4
    x54 = x50 + x52*x53 - 3
    x55 = Fy*x54
    x56 = x49 + x55
    x57 = 2*x1
    x58 = exp(x40 + x57)
    x59 = x43*x58
    x60 = 3*St
    x61 = S*u*xi_x
    x62 = x47*xi_x
    x63 = u*(x4*x62 + x60 + 9*x61)
    x64 = Fx*x54
    x65 = x63 + x64
    x66 = sinh(x42)
    x67 = exp(x1 + x40)
    x68 = x66*x67
    x69 = 2*Fx
    x70 = S*x14
    x71 = xi^3
    x72 = u*x71
    x73 = 1944*xi
    x74 = S*x3
    x75 = S^2
    x76 = x19*x75
    x77 = 972*x23
    x78 = xi^4
    x79 = x14*x78
    x80 = S*x0*x29
    x81 = u^5
    x82 = S*x71
    x83 = x19*x78
    x84 = u^7
    x85 = x75*x84*xi
    x86 = u^11
    x87 = u^8
    x88 = x75*x87
    x89 = u^9
    x90 = x71*x89
    x91 = Sp*u
    x92 = 2*x0
    x93 = Sd*x92
    x94 = -x51
    x95 = x94 + 2
    x96 = x3*x71
    x97 = 7*S
    x98 = x81*xi
    x99 = 7*x29
    x100 = x78*x92
    x101 = Sp*x3
    x102 = -x12
    x103 = 2*x31
    x104 = x103 + 1
    x105 = x102 + x104
    x106 = 9*S
    x107 = 3*x3
    x108 = -x5
    x109 = 8*x5
    x110 = x102 - 3
    x111 = 2*x81
    x112 = 2*x14
    x113 = 648*UU(phi, potential)
    x114 = B2*B2d
    x115 = 1458*B2d*B2p
    x116 = G*Gd
    x117 = 486*Gd*Gp
    x118 = 2592*UU(phi, potential)
    x119 = 23328*B2*B2d
    x120 = x84*xi
    x121 = 7776*G*Gd
    x122 = Gd*Gp
    x123 = 1944*x19*xi
    x124 = x78*x89
    x125 = 11664*x29*x87
    x126 = 1944*x24*x78
    x127 = 2916*x29*x84
    x128 = 1944*x71*x87
    x129 = x43^2
    x130 = B1*B1d*x129
    x131 = 486*B1d*B1p*x129
    x132 = 7776*B1*B1d*x129
    x133 = B1d*B1p*x129
    x134 = 3*phi*u
    x135 = B2*u
    x136 = 4*x135
    x137 = B2p - x136
    x138 = B1*u
    x139 = 4*x138
    x140 = B1p - x139
    x141 = 4*UU(phi, potential)
    x142 = x0*x141
    x143 = 12*G*Gd
    x144 = x143*x87
    x145 = 3*Gd*Gp
    x146 = x145*x84
    x147 = x142 + x144 - x146
    x148 = 9*B2d*x137*x84
    x149 = -1944*B1d*x129*x140*x84 - 5832*B2d*x137*x84 + 7776*G*Gd*x87 - 1944*Gd*Gp*x84 + 2592*UU(phi, potential)*x0 - 7776
    x150 = phid*x10
    x151 = u*x105
    x152 = x108 + 1
    x153 = 9*x150 + x152
    x154 = x141*x98
    x155 = x89*xi
    x156 = x143*x155
    x157 = x87*xi
    x158 = x145*x157
    x159 = x94 + 1
    x160 = 12*B1*B1d*x129
    x161 = x160*x87
    x162 = 3*B1d*B1p*x129
    x163 = x162*x84
    x164 = x155*x160
    x165 = x157*x162
    x166 = x148*x6
    x167 = phi*x38
    x168 = phid*x38*x52
    x169 = 12*x5
    x170 = S*x0*x52
    x171 = -x142 - x144 + x146 + x154 + x156 - x158 - x161 + x163 + x164 - x165 - x166
    x172 = Fy*x14
    x173 = x172 + xi_y
    x174 = G*u*x173
    x175 = Gh + 4*x174
    x176 = Fx*x14
    x177 = x176 + xi_x
    x178 = G*u*x177
    x179 = Gt + 4*x178
    x180 = x43*x67
    x181 = 4*x4*xi
    x182 = 4*Fy
    x183 = x30 - 1
    x184 = x183*x53 + 18*x8 - 3
    x185 = Fx*u
    x186 = 4*x185
    x187 = 2*Fy*x14
    x188 = 3*Sp
    x189 = x188*x3
    x190 = -x189
    x191 = x190 + x54
    x192 = 18*Sh*u
    x193 = 3*Sh*x0
    x194 = 2*x4*xi
    x195 = xi_y^2
    x196 = 2*x4
    x197 = Sh*xi_y
    x198 = 12*B1*x81
    x199 = 9*S*x81*xi_y
    x200 = 12*B2*x81
    x201 = 36*S*x14
    x202 = 36*B1*S*x19
    x203 = 36*B2*S*x19
    x204 = 6*u*x4*xi
    x205 = 2*B1h
    x206 = x0*x4*xi*xi_y
    x207 = 2*B2h
    x208 = 8*B1*x4*x81*xi
    x209 = 8*B2*x4*x81*xi
    x210 = Fy^2
    x211 = 6*B2
    x212 = x0*x211
    x213 = 6*B1
    x214 = -x0*x213
    x215 = 27*x8
    x216 = B1*S*x87
    x217 = 18*x216
    x218 = B2*S*x87
    x219 = 18*x218
    x220 = 5*x5
    x221 = x52*x57
    x222 = 2*x40
    x223 = 2*Sh
    x224 = 9*x3
    x225 = x3*x54
    x226 = B2h*x225
    x227 = x106*x84
    x228 = x4*x81
    x229 = 2*xi_y
    x230 = x214 + 36*x216
    x231 = 45*x8
    x232 = 36*x218
    x233 = x183*x57
    x234 = 18*St*u
    x235 = 3*St*x0
    x236 = xi_x^2
    x237 = St*xi_x
    x238 = 9*S*x81*xi_x
    x239 = 2*B1t
    x240 = x0*x4*xi*xi_x
    x241 = 2*B2t
    x242 = Fx^2
    x243 = -x212
    x244 = x220 + x222*x52 - 2
    x245 = x109 + x183*x222 - 1
    x246 = 6*(Fy*xi_x)
    x247 = B2t*Fy
    x248 = x4*xi_x*xi_y
    x249 = 72*xi_x*xi_y
    x250 = 16*x4
    x251 = 2*Fy
    x252 = u*x251
    x253 = B1h + x139*x173
    x254 = 2*x0*x41*x43
    x255 = B2h + x136*x173
    x256 = B1t + x139*x177
    x257 = 2*x0*x43*x58
    x258 = B2t + x136*x177
    x259 = 8*u*x173
    x260 = 5*B1*u
    x261 = 5*B2*u
    x262 = x41*x66
    x263 = 5*x174
    x264 = 2*B2p - 8*x135
    x265 = Fyh*x14 + x111*x210 + x251*x3*xi_y + xi_yy
    x266 = x265*x41*x43
    x267 = G*u
    x268 = 2*Gp - 8*x267
    x269 = x175*x41*x66
    x270 = 2*B1p - 8*x138
    x271 = Fxh*x14 + x173*x3*x69 + xi_xy
    x272 = x264*x66*x67
    x273 = x268*x43*x67
    x274 = x3*xi_x
    x275 = Fx*x251*x81 + Fyt*x14 + x251*x274 + xi_xy
    x276 = 2*x0*x175*x43*x67
    x277 = 8*x66*x67
    x278 = 4*x0*x179
    x279 = 8*u*x177
    x280 = x177*x261
    x281 = x58*x66
    x282 = Fxt*x14 + x111*x242 + x274*x69 + xi_xx
    x283 = x282*x43*x58
    x284 = 4*x14*x4
    x285 = u - x14*x47
    x286 = Fy*x285 + x4*(phih + x134*x173) - x48
    x287 = 16*u
    x288 = Fx*x285 + x4*(phit + x134*x177) - x62
    x289 = 5*x267*xi_x
    ABCS[1] = 2*x2*(3*x11 + x3*x7)^4/27
    ABCS[2] = x17*x3
    ABCS[3] = x14*x17
    ABCS[4] = 3*x0*(x10 + x15/3)^2*(2*x0*x179*x180*x255 - 2*x0*x179*x253*x43*x67 - 2*x0*x179*x258*x281 + 4*x0*x253*x269 - 4*x0*x255^2*x44 + x0*x255*x258*x277 - 2*x0*x255*x269 - 4*x0*x258^2*x59 + x14*x277*x286*x288*x4 - x175^2*x254 + x175*x278*x66*x67 - x179^2*x257 + x180*(4*Gc + x287*(Gh*xi_x + Gt*xi_y + x172*(Gt + x289) + x176*(Gh + x263) + x289*xi_y)) - x253^2*x254 + x253*x254*x255 - x256^2*x257 - x256*x257*x258 + x256*x276 - x256*x278*x58*x66 + x258*x276 + x262*x265*x268 - x262*(2*Gs + x259*(2*Gh + x263)) + x264*x266 + x264*x283 - x266*x270 + x268*x281*x282 + x270*x283 - x271*x272 - x271*x273 - x272*x275 - x273*x275 - x281*(2*Gb + x279*(2*Gt + 5*x178)) - x284*x286^2*x41*x43 - x284*x288^2*x43*x58 + x44*(2*B1s + x259*(x173*x260 + x205)) - x44*(2*B2s + x259*(x173*x261 + x207)) + x44*(Fyp - x252)^2 - x59*(2*B1b + x279*(x177*x260 + x239)) - x59*(2*B2b + x279*(x241 + x280)) + x59*(Fxp - 2*x185)^2 + x68*(4*B2c + x287*(B2h*x177 + x173*(B2t + x280))) + x68*(2*Fxp - x186)*(-Fyp + x252)) - 8*x16*x3*(-x0*x180*(G*u*x182*(x184*xi_x + x3*x60) + 4*G*x14*(x45*xi_x + xi_y*(x181*xi_x + x60 + 18*x61)) + G*x186*(x184*xi_y + x187*x54 + x3*x45) + Gh*x63 + Gh*x64 + Gt*x49 + Gt*x55) + x0*x66*(x175*x41*x56 + x179*x58*x65) + x44*(Fyh*x191 + u*(-B1h*x193 - B1h*x199 + B2h*x193 + B2h*x199 + Fy*(B1h*(x107 - x19*x194 - x227 + x228) - x223*(-x224 + x84*(-x211 + x213)) + x226 - x229*(x189 + x212 + x230 - x231 - x232 + x53*(-x109 + x222*(-x30 + 1) + x233 + 1) + 3)) + 3*Ss - x112*x210*(x189 + x212 + x214 - x215 + x217 - x219 + x53*(x159*x222 - x220 + x221 + 2) + 3) - x188*xi_yy + x192*xi_y + x194*xi_yy + x195*x196 + x195*x201 - x195*x202 + x195*x203 + x195*x204 - x195*x208 + x195*x209 - x197*x198 + x197*x200 - x205*x206 + x206*x207 + x46*xi_yy)) + x59*(Fxt*x191 + u*(B1t*x235 + B1t*x238 + B2t*x235 + B2t*x238 + Fx*(B1t*x225 + B2t*x225 + 2*St*(x224 + x84*(x211 + x213)) + 2*xi_x*(x190 + x230 + x231 + x232 + x243 + x53*(x233 + x245) - 3)) + 3*Sb + x112*x242*(x190 + x214 + x215 + x217 + x219 + x243 + x53*(x221 + x244) - 3) - x188*xi_xx + x194*xi_xx + x196*x236 + x198*x237 + x200*x237 + x201*x236 + x202*x236 + x203*x236 + x204*x236 + x208*x236 + x209*x236 + x234*xi_x + x239*x240 + x240*x241 + x46*xi_xx)) - x68*(Fxh*x191 + Fyt*x191 + u*(12*B2*Fy*St*x84 - 4*B2*Fy*x19*x4*xi_x + B2*Fy*x250*x84*xi*xi_x + B2*S*x19*x249 + 16*B2*x4*x81*xi*xi_x*xi_y + B2h*x235 + B2h*x238 + 2*B2t*Fy*x19*x4*xi + B2t*x193 + B2t*x199 + Fx*(x182*(x0*x244*x4 - 3*x14*(-S*x211*x87 + x101 + x222 - x50 + 1)) + x223*(x211*x84 + x224) + x226 + x229*(x190 + x231 + x232 + x243 + x245*x53 - 3)) + 18*Fy*St*x3 + 72*Fy*x218*xi_x + Fy*x250*x3*xi*xi_x - 12*Fy*x40*xi_x + 90*Fy*x8*xi_x + 6*Sc + Sh*x200*xi_x - 6*Sp*xi_xy + St*x200*xi_y - x101*x246 - x107*x247 + x169*x248 + x181*xi_xy - x187*x4*xi_x + x192*xi_x + x206*x241 + x207*x240 + x227*x247 - x228*x247 + x234*xi_y - x246 + 4*x248 + x249*x70 + (18*xi_xy)*(S*u))))/3 + 4*x19*(x44*x56^2 - x56*x68*(x54*x69 + 2*x63) + x59*x65^2)/3 + x2*(48*x10*x34*x35 + 324*x18 + x20*x21 + 1296*x22 + 1296*x25 + x26*x28 + 324*x29*(x30 + x31 + 6) + 432*x32*x7 + 216*x36*x37*x39)/27 - x2*(-972*S*x83 - 648*Sd*x39*(x13 - 1) - u^12*x29*x77 + 36*x0*x36*(-Sd*x112*x6*(4*x31 + x6 + x8*(7*x5 - 5)) - x111*x82 - 4*x3*xi^5 + x47*x74 + x6*x91*(x110 + x31 + x6*x93 - 2*x8 + x96) - 11*x70 + 6*x72 - 7*x76 - 5*x79 + 22*x80 - x83*x97 + 8*x85 + x99) + 108*x10*x4*(-x100 + x101*(4*Sd*x0*x6 + x103 + x110 - x8 + 2*x96) - x106*x19*x29 + x106*x98 - x107*x71 - x31 + x50 - 8*x82*x84 + 3*x88 - x93*(x104 + x108 + x8*(x109 - 7)) + x94 - 2) - x23*x73*x86 - x24*x77 + x28*x37*x52 - 5508*x29*x88 + 1944*x29 - 12*x34*x6*(x100 + x101*x6 + x105*x93 - x14*x99 + 5*x8 + x95 + x96 - x97*x98) + 324*x38*x91*(x9^2 + x93) - 324*x70 + 1296*x72 - x73*x74 - 1944*x75*x90 - 1620*x76 + 324*x79 - 3888*x80 - 3240*x81*x82 - 5184*x85)/9 + x2*(-5832*B2*B2d*x24*x78 + 5832*B2d*B2p*x19*xi + 8748*B2d*B2p*x29*x84 + 5832*B2d*B2p*x71*x87 - 3888*UU(phi, potential)*x0*x29 + 24*phi0^14*phid*x26*x86*(phip - x134) - 12*phi0^12*x35*x89*(phi*u*(x102 + 72*x150 + 3) - 2*phid*x151 - 24*phip*x150 + phip*x5 - phip) + 24*phi0^10*x37*x84*(-18*phi*x11*x153 + 6*phip*x10*x153 + x151*(12*x150 + x152)) - 108*u*x36*(9*phip*x10^4 + u*x38*(-27*x167 - 18*x168 + x6*(x169 + 24*x170 + x171 + 48*x31 - 12))) + 24*x10*x3*x33*(27*phip*x38*(3*x150 + x95) - u*(x167*(243*phid*x10 - 162*u*xi + 162) + x6*(-108*x168 + x6*(54*x170 + x171 + 108*x31 + 42*x5 - 42)))) + x111*x27*x6*(324*phip*x38*(4*x150 + x152) - u*(x167*(3888*phid*x10 - 972*u*xi + 972) + x6*(-648*x168 + x6*(144*x170 + x171 + 288*x31 + 132*x5 - 132)))) - x113*x14 - x113*x83 - 5832*x114*x19 - 34992*x114*x29*x87 + x115*x124 + x115*x81 - x116*x125 - x116*x126 - x116*x20 + x117*x124 + x117*x81 - x118*x3*xi - x118*x71*x81 - x119*x120 - x119*x90 - x120*x121 - x120*x132 - x121*x90 + x122*x123 + x122*x127 + x122*x128 + x123*x133 + x124*x131 - x125*x130 - x126*x130 + x127*x133 + x128*x133 - x130*x20 + x131*x81 - x132*x90 - x149*x22 - x149*x25 - x18*(-486*B1d*x129*x140*x84 - 1458*B2d*x137*x84 + 1944*G*Gd*x87 - 486*Gd*Gp*x84 + 648*UU(phi, potential)*x0 - 1944) - 972*x19*x21*(-3*B1d*x129*x140*x84 + x147 - x148 - 12) + 11664*x29 + 216*x32*x4*(x12 + x147 - x154 - x156 + x158 + x159*x50 + x161 - x163 - x164 + x165 + x166 - 18*x31 - 3) + 7776*x72 + 1944*x79)/81

    nothing
end

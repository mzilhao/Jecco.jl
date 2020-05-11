
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

function S_eq_coeff!(ABCS::Vector, vars::SVars, ::Inner)
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

function Fxy_eq_coeff!(AA::Matrix, BB::Matrix, CC::Matrix, SS::Vector, vars::FVars, ::Inner)
    @unpack (
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


function Sd_eq_coeff!(ABCS::Vector, vars::SdVars, ::Inner)
    @unpack (
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

    ABCS[4] = *(u, *(-4, (xi_x + *(u3, St + *(Sp, xi_x + *(Fx, u2))) + *(-1, xi_x + *(Fx, u2), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3)) + *(1/3, xi_x, phi0 ^ 2, u2)) ^ 2, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-4, (xi_y + *(u3, Sh + *(Sp, xi_y + *(Fy, u2))) + *(-1, xi_y + *(Fy, u2), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3)) + *(1/3, xi_y, phi0 ^ 2, u2)) ^ 2, coshGu4, expB2u4) + *(8, xi_x + *(u3, St + *(Sp, xi_x + *(Fx, u2))) + *(-1, xi_x + *(Fx, u2), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3)) + *(1/3, xi_x, phi0 ^ 2, u2), xi_y + *(u3, Sh + *(Sp, xi_y + *(Fy, u2))) + *(-1, xi_y + *(Fy, u2), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3)) + *(1/3, xi_y, phi0 ^ 2, u2), exp(*(B1, u4) + *(B2, u4)), sinhGu4)) + *(1 + *(S, u4) + *(u, xi) + *(-1/3, phi0 ^ 2, u2) + *(1/3, xi, phi0 ^ 2, u3), *(-16, xi_xy + *(u3, Sc + *(Spt, xi_y + *(Fy, u2)) + *(Sph + *(Spp, xi_y + *(Fy, u2)), xi_x + *(Fx, u2))) + *(-1, xi_x + *(Fx, u2), *(u3, Sph + *(Spp, xi_y + *(Fy, u2))) + *(-3, u4, Sh + *(Sp, xi_y + *(Fy, u2))) + *(-2/3, xi_y, phi0 ^ 2, u3)) + *(-1, xi_y + *(Fy, u2), *(u3, Spt + *(Spp, xi_x + *(Fx, u2))) + *(-1, xi_x + *(Fx, u2), *(Spp, u3) + *(-6, Sp, u4) + *(12, S, u5) + *(-2/3, phi0 ^ 2, u3) + *(2, xi, phi0 ^ 2, u4)) + *(-3, u4, St + *(Sp, xi_x + *(Fx, u2))) + *(-2/3, xi_x, phi0 ^ 2, u3)) + *(1/3, xi_xy, phi0 ^ 2, u2), exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(8, xi_xx + *(u3, Sb + *(Spt, xi_x + *(Fx, u2)) + *(Spt + *(Spp, xi_x + *(Fx, u2)), xi_x + *(Fx, u2))) + *(-1, xi_x + *(Fx, u2), *(u3, Spt + *(Spp, xi_x + *(Fx, u2))) + *(-3, u4, St + *(Sp, xi_x + *(Fx, u2))) + *(-2/3, xi_x, phi0 ^ 2, u3)) + *(-1, xi_x + *(Fx, u2), *(u3, Spt + *(Spp, xi_x + *(Fx, u2))) + *(-1, xi_x + *(Fx, u2), *(Spp, u3) + *(-6, Sp, u4) + *(12, S, u5) + *(-2/3, phi0 ^ 2, u3) + *(2, xi, phi0 ^ 2, u4)) + *(-3, u4, St + *(Sp, xi_x + *(Fx, u2))) + *(-2/3, xi_x, phi0 ^ 2, u3)) + *(1/3, xi_xx, phi0 ^ 2, u2), coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(8, xi_yy + *(u3, Ss + *(Sph, xi_y + *(Fy, u2)) + *(Sph + *(Spp, xi_y + *(Fy, u2)), xi_y + *(Fy, u2))) + *(-1, xi_y + *(Fy, u2), *(u3, Sph + *(Spp, xi_y + *(Fy, u2))) + *(-3, u4, Sh + *(Sp, xi_y + *(Fy, u2))) + *(-2/3, xi_y, phi0 ^ 2, u3)) + *(-1, xi_y + *(Fy, u2), *(u3, Sph + *(Spp, xi_y + *(Fy, u2))) + *(-1, xi_y + *(Fy, u2), *(Spp, u3) + *(-6, Sp, u4) + *(12, S, u5) + *(-2/3, phi0 ^ 2, u3) + *(2, xi, phi0 ^ 2, u4)) + *(-3, u4, Sh + *(Sp, xi_y + *(Fy, u2))) + *(-2/3, xi_y, phi0 ^ 2, u3)) + *(1/3, xi_yy, phi0 ^ 2, u2), coshGu4, expB2u4) + *(-8, *(u4, B1h + *(B1p, xi_y + *(Fy, u2))) + *(-1, xi_y + *(Fy, u2), *(B1p, u4) + *(-4, B1, u5)), xi_y + *(u3, Sh + *(Sp, xi_y + *(Fy, u2))) + *(-1, xi_y + *(Fy, u2), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3)) + *(1/3, xi_y, phi0 ^ 2, u2), coshGu4, expB2u4) + *(-8, *(u4, B2h + *(B2p, xi_y + *(Fy, u2))) + *(-1, xi_y + *(Fy, u2), *(B2p, u4) + *(-4, B2, u5)), xi_x + *(u3, St + *(Sp, xi_x + *(Fx, u2))) + *(-1, xi_x + *(Fx, u2), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3)) + *(1/3, xi_x, phi0 ^ 2, u2), exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-8, *(u4, B2t + *(B2p, xi_x + *(Fx, u2))) + *(-1, xi_x + *(Fx, u2), *(B2p, u4) + *(-4, B2, u5)), xi_y + *(u3, Sh + *(Sp, xi_y + *(Fy, u2))) + *(-1, xi_y + *(Fy, u2), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3)) + *(1/3, xi_y, phi0 ^ 2, u2), exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-8, *(u4, Gh + *(Gp, xi_y + *(Fy, u2))) + *(-1, xi_y + *(Fy, u2), *(Gp, u4) + *(-4, G, u5)), xi_x + *(u3, St + *(Sp, xi_x + *(Fx, u2))) + *(-1, xi_x + *(Fx, u2), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3)) + *(1/3, xi_x, phi0 ^ 2, u2), coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-8, *(u4, Gt + *(Gp, xi_x + *(Fx, u2))) + *(-1, xi_x + *(Fx, u2), *(Gp, u4) + *(-4, G, u5)), xi_y + *(u3, Sh + *(Sp, xi_y + *(Fy, u2))) + *(-1, xi_y + *(Fy, u2), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3)) + *(1/3, xi_y, phi0 ^ 2, u2), coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-8, xi_xx + *(u2, Fxt + *(Fxp, xi_x + *(Fx, u2))) + *(-1, xi_x + *(Fx, u2), *(Fxp, u2) + *(-2, Fx, u3)), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3), coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-8, xi_yy + *(u2, Fyh + *(Fyp, xi_y + *(Fy, u2))) + *(-1, xi_y + *(Fy, u2), *(Fyp, u2) + *(-2, Fy, u3)), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3), coshGu4, expB2u4) + *(-2, *(Fxp, u2) + *(-2, Fx, u3), xi_x + *(u3, St + *(Sp, xi_x + *(Fx, u2))) + *(-1, xi_x + *(Fx, u2), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3)) + *(1/3, xi_x, phi0 ^ 2, u2), coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-2, *(Fyp, u2) + *(-2, Fy, u3), xi_y + *(u3, Sh + *(Sp, xi_y + *(Fy, u2))) + *(-1, xi_y + *(Fy, u2), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3)) + *(1/3, xi_y, phi0 ^ 2, u2), coshGu4, expB2u4) + *(2, *(Fxp, u2) + *(-2, Fx, u3), xi_y + *(u3, Sh + *(Sp, xi_y + *(Fy, u2))) + *(-1, xi_y + *(Fy, u2), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3)) + *(1/3, xi_y, phi0 ^ 2, u2), exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(2, *(Fyp, u2) + *(-2, Fy, u3), xi_x + *(u3, St + *(Sp, xi_x + *(Fx, u2))) + *(-1, xi_x + *(Fx, u2), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3)) + *(1/3, xi_x, phi0 ^ 2, u2), exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(8, *(u4, B1t + *(B1p, xi_x + *(Fx, u2))) + *(-1, xi_x + *(Fx, u2), *(B1p, u4) + *(-4, B1, u5)), xi_x + *(u3, St + *(Sp, xi_x + *(Fx, u2))) + *(-1, xi_x + *(Fx, u2), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3)) + *(1/3, xi_x, phi0 ^ 2, u2), coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(8, *(u4, B2h + *(B2p, xi_y + *(Fy, u2))) + *(-1, xi_y + *(Fy, u2), *(B2p, u4) + *(-4, B2, u5)), xi_y + *(u3, Sh + *(Sp, xi_y + *(Fy, u2))) + *(-1, xi_y + *(Fy, u2), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3)) + *(1/3, xi_y, phi0 ^ 2, u2), coshGu4, expB2u4) + *(8, *(u4, B2t + *(B2p, xi_x + *(Fx, u2))) + *(-1, xi_x + *(Fx, u2), *(B2p, u4) + *(-4, B2, u5)), xi_x + *(u3, St + *(Sp, xi_x + *(Fx, u2))) + *(-1, xi_x + *(Fx, u2), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3)) + *(1/3, xi_x, phi0 ^ 2, u2), coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(8, *(u4, Gh + *(Gp, xi_y + *(Fy, u2))) + *(-1, xi_y + *(Fy, u2), *(Gp, u4) + *(-4, G, u5)), xi_y + *(u3, Sh + *(Sp, xi_y + *(Fy, u2))) + *(-1, xi_y + *(Fy, u2), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3)) + *(1/3, xi_y, phi0 ^ 2, u2), expB2u4, sinhGu4) + *(8, *(u4, Gt + *(Gp, xi_x + *(Fx, u2))) + *(-1, xi_x + *(Fx, u2), *(Gp, u4) + *(-4, G, u5)), xi_x + *(u3, St + *(Sp, xi_x + *(Fx, u2))) + *(-1, xi_x + *(Fx, u2), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3)) + *(1/3, xi_x, phi0 ^ 2, u2), exp(*(B2, u4) + *(2, B1, u4)), sinhGu4) + *(8, xi_xy + *(u2, Fxh + *(Fxp, xi_y + *(Fy, u2))) + *(-1, xi_y + *(Fy, u2), *(Fxp, u2) + *(-2, Fx, u3)), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3), exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(8, xi_xy + *(u2, Fyt + *(Fyp, xi_x + *(Fx, u2))) + *(-1, xi_x + *(Fx, u2), *(Fyp, u2) + *(-2, Fy, u3)), 1 + *(Sp, u3) + *(-3, S, u4) + *(1/3, phi0 ^ 2, u2) + *(-2/3, xi, phi0 ^ 2, u3), exp(*(B1, u4) + *(B2, u4)), sinhGu4)) + *(-4/27, *(-324, xi ^ 3) + *(phi0 ^ 8, u5) + *(-81, u, xi ^ 4) + *(-6, phi0 ^ 6, u3) + *(108, xi, phi0 ^ 2) + *(-81, S ^ 3, u ^ 9, *(-3, (1 + *(u, xi)) ^ 2) + *(phi0 ^ 2, u2)) + *(-63, phi0 ^ 4, u3, xi ^ 2) + *(-54, phi0 ^ 4, u4, xi ^ 3) + *(-24, phi0 ^ 6, u6, xi ^ 3) + *(-9, S ^ 2, u5, *(-9, (1 + *(u, xi)) ^ 2, 5 + *(6, u, xi)) + *(phi0 ^ 4, u4, -7 + *(8, u, xi)) + *(-3, phi0 ^ 2, u2, -12 + *(-12, u, xi) + *(8, u3, xi ^ 3) + *(9, u2, xi ^ 2))) + *(-4, xi, phi0 ^ 8, u6) + *(-3, S, u, *(-27, (1 + *(u, xi)) ^ 3, 1 + *(3, u, xi)) + *(phi0 ^ 6, u6, 5 + *(-12, u, xi) + *(7, u2, xi ^ 2)) + *(-9, phi0 ^ 2, u2, -7 + *(u2, xi ^ 2) + *(-16, u, xi) + *(10, u4, xi ^ 4) + *(20, u3, xi ^ 3)) + *(-3, phi0 ^ 4, u4, 11 + *(-22, u2, xi ^ 2) + *(-2, u, xi) + *(2, u3, xi ^ 3) + *(7, u4, xi ^ 4))) + *(-3, phi0 ^ 6, u ^ 7, xi ^ 4) + *(-2, phi0 ^ 8, u8, xi ^ 3) + *(3, Sp, (3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi))) ^ 2, *(-3, (1 + *(u, xi)) ^ 2) + *(phi0 ^ 2, u2)) + *(5, phi0 ^ 8, u ^ 7, xi ^ 2) + *(6, phi0 ^ 6, u8, xi ^ 5) + *(12, xi, phi0 ^ 6, u4) + *(15, phi0 ^ 6, u5, xi ^ 2) + *(36, phi0 ^ 4, u6, xi ^ 5) + *(45, phi0 ^ 4, u5, xi ^ 4) + *(54, phi0 ^ 2, u4, xi ^ 5) + *(81, u, phi0 ^ 2, xi ^ 2) + *(108, phi0 ^ 2, u2, xi ^ 3) + *(135, phi0 ^ 2, u3, xi ^ 4), expB1u4) + *(4/9, *(27, phi0 ^ 2, *(xi, -2 + *(u3, xi ^ 3) + *(2, u2, xi ^ 2)) + *(S ^ 2, u ^ 7, -1 + *(u2, xi ^ 2)) + *(2, S, u3, (1 + *(u, xi)) ^ 2, -1 + *(u, xi))) + *(27, xi ^ 3, 4 + *(u, xi)) + *(81, u5, (S + *(S, u, xi)) ^ 2) + *(27, S ^ 3, u ^ 9, 1 + *(u, xi)) + *(81, S, u, (1 + *(u, xi)) ^ 3) + *(phi0 ^ 6, u3, (-1 + *(u, xi)) ^ 3, 1 + *(u, xi)) + *(9, u, phi0 ^ 4, (-1 + *(u, xi)) ^ 2, 1 + *(u, xi), 1 + *(S, u4) + *(u, xi)), expB1u4) + *(4/81, *(-1944, xi ^ 3) + *(27, phi0 ^ 2, *(30, xi) + *(u3, *(-8, UU(phi, potential)) + *(39, xi ^ 4)) + *(-2, u4, *(-9, xi ^ 5) + *(8, UU(phi, potential), xi)) + *(18, u, xi ^ 2) + *(24, u2, xi ^ 3) + *(8, UU(phi, potential), u ^ 7, xi ^ 4) + *(16, UU(phi, potential), u6, xi ^ 3)) + *(162, u, UU(phi, potential) + *(-3, xi ^ 4)) + *(81, S ^ 4, u ^ 13, -6 + *(2, UU(phi, potential), u4) + *(3, phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(162, UU(phi, potential), u5, xi ^ 4) + *(648, UU(phi, potential), xi, u2) + *(648, UU(phi, potential), u4, xi ^ 3) + *(972, UU(phi, potential), u3, xi ^ 2) + *(2, phi0 ^ 8, u5, (-1 + *(u, xi)) ^ 3, -15 + *(-1, UU(phi, potential), u4) + *(15, u, xi) + *(36, u2, xi ^ 2) + *(UU(phi, potential), xi, u5)) + *(3, phi0 ^ 10, u ^ 7, (-1 + *(u, xi)) ^ 4, -1 + *(2, u, xi)) + *(12, S, u, (3 + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi))) ^ 3, -6 + *(2, UU(phi, potential), u4) + *(3, phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(54, S ^ 2, u5, (3 + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi))) ^ 2, -6 + *(2, UU(phi, potential), u4) + *(3, phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(108, S ^ 3, u ^ 9, -6 + *(2, UU(phi, potential), u4) + *(3, phi0 ^ 2, u2, -1 + *(2, u, xi)), 3 + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi))) + *(6, phi0 ^ 6, u3, (-1 + *(u, xi)) ^ 2, 1 + *(u, xi), -15 + *(-4, UU(phi, potential), u4) + *(15, u, xi) + *(54, u2, xi ^ 2) + *(4, UU(phi, potential), xi, u5)) + *(108, phi0 ^ 4, u3, (1 + *(u, xi)) ^ 2, -1 + *(u, xi), *(6, xi ^ 2) + *(-1, UU(phi, potential), u2) + *(UU(phi, potential), xi, u3)), expB1u4) + *(u, (1 + *(S, u4) + *(u, xi) + *(-1/3, phi0 ^ 2, u2) + *(1/3, xi, phi0 ^ 2, u3)) ^ 2, *(u2, *(Fxp ^ 2, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(Fyp ^ 2, coshGu4, expB2u4) + *(-4, B2c, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-4, Gc, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-2, B1s, coshGu4, expB2u4) + *(2, B1b, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(2, B2s, coshGu4, expB2u4) + *(2, B2b, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(2, Gs, expB2u4, sinhGu4) + *(2, Gb, exp(*(B2, u4) + *(2, B1, u4)), sinhGu4) + *(-12, Fx, xi_y, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-12, Fy, xi_x, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-2, B1p, xi_xx, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-2, B2p, xi_xx, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-2, B2p, xi_yy, coshGu4, expB2u4) + *(-2, Fxp, Fyp, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-2, Gp, xi_xx, exp(*(B2, u4) + *(2, B1, u4)), sinhGu4) + *(-2, Gp, xi_yy, expB2u4, sinhGu4) + *(2, B1p, xi_yy, coshGu4, expB2u4) + *(4, B2p, xi_xy, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(4, Gp, xi_xy, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(12, Fx, xi_x, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(12, Fy, xi_y, coshGu4, expB2u4)) + *(-8, u3, *(*(*(B1, xi_yy) + *(Fy, Fyp) + *(-1, B2, xi_yy) + *(-2, B2h, xi_y) + *(2, B1h, xi_y), expB2u4) + *(*(Fx, Fxp) + *(-1, B1, xi_xx) + *(-1, B2, xi_xx) + *(-2, B1t, xi_x) + *(-2, B2t, xi_x), exp(*(B2, u4) + *(2, B1, u4))) + *(2, *(G, xi_xy) + *(Gh, xi_x) + *(Gt, xi_y), exp(*(B1, u4) + *(B2, u4))), coshGu4) + *(-1, *(*(Fx, Fyp) + *(Fxp, Fy) + *(-2, B2, xi_xy) + *(-2, B2h, xi_x) + *(-2, B2t, xi_y), exp(*(B1, u4) + *(B2, u4))) + *(G, xi_xx, exp(*(B2, u4) + *(2, B1, u4))) + *(G, xi_yy, expB2u4) + *(2, Gh, xi_y, expB2u4) + *(2, Gt, xi_x, exp(*(B2, u4) + *(2, B1, u4))), sinhGu4)) + *(-4, u ^ 7, *(B2p, *(Fx ^ 2, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(Fy ^ 2, coshGu4, expB2u4) + *(-2, Fx, Fy, exp(*(B1, u4) + *(B2, u4)), sinhGu4)) + *(B1p, *(Fx ^ 2, exp(*(B2, u4) + *(2, B1, u4))) + *(-1, Fy ^ 2, expB2u4), coshGu4) + *(Gp, Fx ^ 2, exp(*(B2, u4) + *(2, B1, u4)), sinhGu4) + *(Gp, Fy ^ 2, expB2u4, sinhGu4) + *(-8, B2, B2h, xi_y, coshGu4, expB2u4) + *(-8, B2, B2t, xi_x, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-4, B1, B1h, xi_y, coshGu4, expB2u4) + *(-4, B1, B1t, xi_x, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-4, B1, Gt, xi_x, exp(*(B2, u4) + *(2, B1, u4)), sinhGu4) + *(-4, B1t, G, xi_x, exp(*(B2, u4) + *(2, B1, u4)), sinhGu4) + *(-4, G, Gh, xi_y, coshGu4, expB2u4) + *(-4, G, Gt, xi_x, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-2, B1, B2t, xi_x, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-2, B1, Fy, Fyp, coshGu4, expB2u4) + *(-2, B1, Gt, xi_y, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-2, B1h, G, xi_x, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-2, B1t, B2, xi_x, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-2, B2, Fx, Fyp, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-2, B2, Fxp, Fy, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-2, B2, Gh, xi_y, expB2u4, sinhGu4) + *(-2, B2, Gt, xi_x, exp(*(B2, u4) + *(2, B1, u4)), sinhGu4) + *(-2, B2h, G, xi_y, expB2u4, sinhGu4) + *(-2, B2t, G, xi_x, exp(*(B2, u4) + *(2, B1, u4)), sinhGu4) + *(-2, Fx, Fy, Gp, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-2, Fx, Fyp, G, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-2, Fxp, Fy, G, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(2, B1, B2h, xi_y, coshGu4, expB2u4) + *(2, B1, Fx, Fxp, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(2, B1, Gh, xi_x, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(2, B1h, B2, xi_y, coshGu4, expB2u4) + *(2, B1t, G, xi_y, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(2, B2, Fx, Fxp, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(2, B2, Fy, Fyp, coshGu4, expB2u4) + *(2, B2, Gh, xi_x, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(2, B2, Gt, xi_y, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(2, B2h, G, xi_x, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(2, B2t, G, xi_y, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(2, Fx, Fxp, G, exp(*(B2, u4) + *(2, B1, u4)), sinhGu4) + *(2, Fy, Fyp, G, expB2u4, sinhGu4) + *(4, B1, Gh, xi_y, expB2u4, sinhGu4) + *(4, B1h, G, xi_y, expB2u4, sinhGu4) + *(4, G, Gh, xi_x, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(4, G, Gt, xi_y, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(4, xi, Fx ^ 2, phi0 ^ 2, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(4, xi, Fy ^ 2, phi0 ^ 2, coshGu4, expB2u4) + *(8, B2, B2h, xi_x, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(8, B2, B2t, xi_y, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-8, Fx, Fy, xi, phi0 ^ 2, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-6, Fx, phi, phit, phi0 ^ 6, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-6, Fy, phi, phih, phi0 ^ 6, coshGu4, expB2u4) + *(6, Fx, phi, phih, phi0 ^ 6, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(6, Fy, phi, phit, phi0 ^ 6, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-24, Fx, phi, xi, xi_y, phi0 ^ 4, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-24, Fy, phi, xi, xi_x, phi0 ^ 4, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(24, Fx, phi, xi, xi_x, phi0 ^ 4, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(24, Fy, phi, xi, xi_y, phi0 ^ 4, coshGu4, expB2u4)) + *(2, u4, *(*(*(Fxh, Gp) + *(Fxp, Gh) + *(Fyp, Gt) + *(Fyt, Gp) + *(-40, G, xi_x, xi_y), exp(*(B1, u4) + *(B2, u4))) + *(*(8, Fx ^ 2) + *(-1, B1p, Fxt) + *(-1, B1t, Fxp) + *(-1, B2p, Fxt) + *(-1, B2t, Fxp) + *(2, phi0 ^ 6, phit ^ 2) + *(20, B1, xi_x ^ 2) + *(20, B2, xi_x ^ 2) + *(8, phi0 ^ 2, xi ^ 2, xi_x ^ 2) + *(-8, phit, xi, xi_x, phi0 ^ 4), exp(*(B2, u4) + *(2, B1, u4))) + *(*(8, Fy ^ 2) + *(B1h, Fyp) + *(B1p, Fyh) + *(-1, B2h, Fyp) + *(-1, B2p, Fyh) + *(-20, B1, xi_y ^ 2) + *(2, phi0 ^ 6, phih ^ 2) + *(20, B2, xi_y ^ 2) + *(8, phi0 ^ 2, xi ^ 2, xi_y ^ 2) + *(-8, phih, xi, xi_y, phi0 ^ 4), expB2u4), coshGu4) + *(*(B2h, Fxp, exp(*(B1, u4) + *(B2, u4))) + *(B2p, Fxh + Fyt, exp(*(B1, u4) + *(B2, u4))) + *(B2t, Fyp, exp(*(B1, u4) + *(B2, u4))) + *(-1, Fxp, Gt, exp(*(B2, u4) + *(2, B1, u4))) + *(-1, Fxt, Gp, exp(*(B2, u4) + *(2, B1, u4))) + *(-1, Fyh, Gp, expB2u4) + *(-1, Fyp, Gh, expB2u4) + *(-16, Fx, Fy, exp(*(B1, u4) + *(B2, u4))) + *(20, G, xi_x ^ 2, exp(*(B2, u4) + *(2, B1, u4))) + *(20, G, xi_y ^ 2, expB2u4) + *(-40, B2, xi_x, xi_y, exp(*(B1, u4) + *(B2, u4))) + *(-4, phih, phit, phi0 ^ 6, exp(*(B1, u4) + *(B2, u4))) + *(-16, xi_x, xi_y, phi0 ^ 2, xi ^ 2, exp(*(B1, u4) + *(B2, u4))) + *(8, phih, xi, xi_x, phi0 ^ 4, exp(*(B1, u4) + *(B2, u4))) + *(8, phit, xi, xi_y, phi0 ^ 4, exp(*(B1, u4) + *(B2, u4))), sinhGu4)) + *(2, u6, *(B1t, *(B2t, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-1, Gh, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(2, Gt, exp(*(B2, u4) + *(2, B1, u4)), sinhGu4)) + *(B2t, *(Gt, exp(*(B2, u4) + *(2, B1, u4)), sinhGu4) + *(-1, Gh, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-4, B2h, exp(*(B1, u4) + *(B2, u4)), sinhGu4)) + *(B1h ^ 2, coshGu4, expB2u4) + *(B1t ^ 2, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(Gh ^ 2, coshGu4, expB2u4) + *(Gt ^ 2, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(2, B2h ^ 2, coshGu4, expB2u4) + *(2, B2t ^ 2, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(B1h, Gt, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(B2h, Gh, expB2u4, sinhGu4) + *(-1, B1h, B2h, coshGu4, expB2u4) + *(-1, B2h, Gt, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-2, B1h, Gh, expB2u4, sinhGu4) + *(-2, Gh, Gt, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(2, Fx ^ 2, phi0 ^ 2, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(2, Fy ^ 2, phi0 ^ 2, coshGu4, expB2u4) + *(-56, B1, Fy, xi_y, coshGu4, expB2u4) + *(-56, B2, Fx, xi_y, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-56, B2, Fy, xi_x, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-56, Fx, G, xi_y, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-56, Fy, G, xi_x, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-4, Fx, Fy, phi0 ^ 2, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(18, phi ^ 2, phi0 ^ 6, xi_x ^ 2, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(18, phi ^ 2, phi0 ^ 6, xi_y ^ 2, coshGu4, expB2u4) + *(56, B1, Fx, xi_x, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(56, B2, Fx, xi_x, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(56, B2, Fy, xi_y, coshGu4, expB2u4) + *(56, Fx, G, xi_x, exp(*(B2, u4) + *(2, B1, u4)), sinhGu4) + *(56, Fy, G, xi_y, expB2u4, sinhGu4) + *(-36, xi_x, xi_y, phi ^ 2, phi0 ^ 6, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-16, Fx, xi_y, phi0 ^ 2, xi ^ 2, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-16, Fy, xi_x, phi0 ^ 2, xi ^ 2, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-12, Fx, phi, xi_y, phi0 ^ 4, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-12, Fy, phi, xi_x, phi0 ^ 4, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-8, Fx, phit, xi, phi0 ^ 4, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-8, Fy, phih, xi, phi0 ^ 4, coshGu4, expB2u4) + *(8, Fx, phih, xi, phi0 ^ 4, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(8, Fy, phit, xi, phi0 ^ 4, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(12, Fx, phi, xi_x, phi0 ^ 4, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(12, Fy, phi, xi_y, phi0 ^ 4, coshGu4, expB2u4) + *(16, Fx, xi_x, phi0 ^ 2, xi ^ 2, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(16, Fy, xi_y, phi0 ^ 2, xi ^ 2, coshGu4, expB2u4)) + *(4, u, *(Fxt + *(-1, Fxp, xi_x), coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(Fyh + *(-1, Fyp, xi_y), coshGu4, expB2u4) + *(-1, Fxh + Fyt + *(-1, Fxp, xi_y) + *(-1, Fyp, xi_x), exp(*(B1, u4) + *(B2, u4)), sinhGu4)) + *(4, u5, *(-5, B1h, Fy, coshGu4, expB2u4) + *(-5, B2h, Fx, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-5, B2t, Fy, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-5, Fx, Gh, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-5, Fy, Gt, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-2, B1, Fyh, coshGu4, expB2u4) + *(-2, B2, Fxh, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-2, B2, Fyt, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-2, Fxh, G, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-2, Fyt, G, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(2, B1, Fxt, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(2, B2, Fxt, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(2, B2, Fyh, coshGu4, expB2u4) + *(2, Fxt, G, exp(*(B2, u4) + *(2, B1, u4)), sinhGu4) + *(2, Fyh, G, expB2u4, sinhGu4) + *(5, B1t, Fx, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(5, B2h, Fy, coshGu4, expB2u4) + *(5, B2t, Fx, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(5, Fx, Gt, exp(*(B2, u4) + *(2, B1, u4)), sinhGu4) + *(5, Fy, Gh, expB2u4, sinhGu4) + *(B1p, Fy, xi_y, coshGu4, expB2u4) + *(B2p, Fx, xi_y, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(B2p, Fy, xi_x, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(Fx, Gp, xi_y, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(Fy, Gp, xi_x, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-1, B1p, Fx, xi_x, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-1, B2p, Fx, xi_x, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-1, B2p, Fy, xi_y, coshGu4, expB2u4) + *(-1, Fx, Gp, xi_x, exp(*(B2, u4) + *(2, B1, u4)), sinhGu4) + *(-1, Fy, Gp, xi_y, expB2u4, sinhGu4) + *(-2, B1, Fxp, xi_x, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-2, B2, Fxp, xi_x, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-2, B2, Fyp, xi_y, coshGu4, expB2u4) + *(-2, Fx, phih, phi0 ^ 4, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-2, Fxp, G, xi_x, exp(*(B2, u4) + *(2, B1, u4)), sinhGu4) + *(-2, Fy, phit, phi0 ^ 4, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-2, Fyp, G, xi_y, expB2u4, sinhGu4) + *(2, B1, Fyp, xi_y, coshGu4, expB2u4) + *(2, B2, Fxp, xi_y, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(2, B2, Fyp, xi_x, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(2, Fx, phit, phi0 ^ 4, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(2, Fxp, G, xi_y, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(2, Fy, phih, phi0 ^ 4, coshGu4, expB2u4) + *(2, Fyp, G, xi_x, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-12, phi, xi, phi0 ^ 4, xi_x ^ 2, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-12, phi, xi, phi0 ^ 4, xi_y ^ 2, coshGu4, expB2u4) + *(-6, phi, phih, xi_x, phi0 ^ 6, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-6, phi, phit, xi_y, phi0 ^ 6, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-4, Fx, xi, xi_x, phi0 ^ 2, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-4, Fy, xi, xi_y, phi0 ^ 2, coshGu4, expB2u4) + *(4, Fx, xi, xi_y, phi0 ^ 2, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(4, Fy, xi, xi_x, phi0 ^ 2, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(6, phi, phih, xi_y, phi0 ^ 6, coshGu4, expB2u4) + *(6, phi, phit, xi_x, phi0 ^ 6, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(24, phi, xi, xi_x, xi_y, phi0 ^ 4, exp(*(B1, u4) + *(B2, u4)), sinhGu4)) + *(4, u ^ 10, *(*(Fx, *(xi_x, *(16, B1 ^ 2) + *(16, G ^ 2) + *(32, B2 ^ 2) + *(16, B1, B2)) + *(9, Fx, phi ^ 2, phi0 ^ 6), exp(*(B2, u4) + *(2, B1, u4))) + *(Fy, *(xi_y, *(16, B1 ^ 2) + *(16, G ^ 2) + *(32, B2 ^ 2) + *(-16, B1, B2)) + *(9, Fy, phi ^ 2, phi0 ^ 6), expB2u4) + *(-16, B2, G, *(Fx, xi_y) + *(Fy, xi_x), exp(*(B1, u4) + *(B2, u4))), coshGu4) + *(-2, *(8, Fx, xi_y, G ^ 2 + *(2, B2 ^ 2)) + *(8, Fy, xi_x, G ^ 2 + *(2, B2 ^ 2)) + *(9, Fx, Fy, phi ^ 2, phi0 ^ 6), exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(16, G, *(B2, Fx, xi_x, exp(*(B2, u4) + *(2, B1, u4))) + *(B2, Fy, xi_y, expB2u4) + *(-2, B1, Fy, xi_y, expB2u4) + *(2, B1, Fx, xi_x, exp(*(B2, u4) + *(2, B1, u4))), sinhGu4)) + *(8, u8, *(B2, *(*(9, Fx ^ 2, exp(*(B2, u4) + *(2, B1, u4))) + *(9, Fy ^ 2, expB2u4) + *(-8, G, xi_x, xi_y, exp(*(B1, u4) + *(B2, u4))), coshGu4) + *(4, G, *(xi_x ^ 2, exp(*(B2, u4) + *(2, B1, u4))) + *(xi_y ^ 2, expB2u4), sinhGu4) + *(-18, Fx, Fy, exp(*(B1, u4) + *(B2, u4)), sinhGu4)) + *(8, B2 ^ 2, *(xi_x ^ 2, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(xi_y ^ 2, coshGu4, expB2u4) + *(-2, xi_x, xi_y, exp(*(B1, u4) + *(B2, u4)), sinhGu4)) + *(B1, *(-9, Fy ^ 2, expB2u4) + *(9, Fx ^ 2, exp(*(B2, u4) + *(2, B1, u4))) + *(-4, B2, xi_y ^ 2, expB2u4) + *(4, B2, xi_x ^ 2, exp(*(B2, u4) + *(2, B1, u4))), coshGu4) + *(4, B1 ^ 2, *(xi_x ^ 2, exp(*(B2, u4) + *(2, B1, u4))) + *(xi_y ^ 2, expB2u4), coshGu4) + *(4, G ^ 2, xi_x ^ 2, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(4, G ^ 2, xi_y ^ 2, coshGu4, expB2u4) + *(8, B1, G, *(xi_x ^ 2, exp(*(B2, u4) + *(2, B1, u4))) + *(-1, xi_y ^ 2, expB2u4), sinhGu4) + *(9, G, Fx ^ 2, exp(*(B2, u4) + *(2, B1, u4)), sinhGu4) + *(9, G, Fy ^ 2, expB2u4, sinhGu4) + *(-18, Fx, Fy, G, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-8, xi_x, xi_y, G ^ 2, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(2, Fx ^ 2, phi0 ^ 2, xi ^ 2, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(2, Fy ^ 2, phi0 ^ 2, xi ^ 2, coshGu4, expB2u4) + *(3, phi, Fx ^ 2, phi0 ^ 4, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(3, phi, Fy ^ 2, phi0 ^ 4, coshGu4, expB2u4) + *(-9, Fx, xi_y, phi ^ 2, phi0 ^ 6, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-9, Fy, xi_x, phi ^ 2, phi0 ^ 6, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-6, Fx, Fy, phi, phi0 ^ 4, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-4, Fx, Fy, phi0 ^ 2, xi ^ 2, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(9, Fx, xi_x, phi ^ 2, phi0 ^ 6, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(9, Fy, xi_y, phi ^ 2, phi0 ^ 6, coshGu4, expB2u4)) + *(8, u ^ 9, *(B1, *(B2t, Fx, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(Fy, Gt, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-1, B2h, Fy, coshGu4, expB2u4) + *(-1, Fx, Gh, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-2, Fy, Gh, expB2u4, sinhGu4) + *(2, B1h, Fy, coshGu4, expB2u4) + *(2, B1t, Fx, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(2, Fx, Gt, exp(*(B2, u4) + *(2, B1, u4)), sinhGu4)) + *(B1t, *(B2, Fx, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-1, Fy, G, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(2, Fx, G, exp(*(B2, u4) + *(2, B1, u4)), sinhGu4)) + *(B1h, Fx, G, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(B2, Fx, Gt, exp(*(B2, u4) + *(2, B1, u4)), sinhGu4) + *(B2, Fy, Gh, expB2u4, sinhGu4) + *(B2h, Fy, G, expB2u4, sinhGu4) + *(B2t, Fx, G, exp(*(B2, u4) + *(2, B1, u4)), sinhGu4) + *(-1, B1h, B2, Fy, coshGu4, expB2u4) + *(-1, B2, Fx, Gh, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-1, B2, Fy, Gt, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-1, B2h, Fx, G, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-1, B2t, Fy, G, coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-4, B2, B2h, Fx, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-4, B2, B2t, Fy, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-2, B1h, Fy, G, expB2u4, sinhGu4) + *(-2, Fx, G, Gh, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-2, Fy, G, Gt, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(2, Fx, G, Gt, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(2, Fy, G, Gh, coshGu4, expB2u4) + *(4, B2, B2h, Fy, coshGu4, expB2u4) + *(4, B2, B2t, Fx, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-6, phi, xi, Fx ^ 2, phi0 ^ 4, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-6, phi, xi, Fy ^ 2, phi0 ^ 4, coshGu4, expB2u4) + *(12, Fx, Fy, phi, xi, phi0 ^ 4, exp(*(B1, u4) + *(B2, u4)), sinhGu4)) + *(32, u ^ 12, *(G ^ 2, *(Fx ^ 2, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(Fy ^ 2, coshGu4, expB2u4) + *(-2, Fx, Fy, exp(*(B1, u4) + *(B2, u4)), sinhGu4)) + *(2, B2 ^ 2, *(Fx ^ 2, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(Fy ^ 2, coshGu4, expB2u4) + *(-2, Fx, Fy, exp(*(B1, u4) + *(B2, u4)), sinhGu4)) + *(B1, *(B2, coshGu4) + *(2, G, sinhGu4), *(Fx ^ 2, exp(*(B2, u4) + *(2, B1, u4))) + *(-1, Fy ^ 2, expB2u4)) + *(B2, G, *(Fx ^ 2, exp(*(B2, u4) + *(2, B1, u4)), sinhGu4) + *(Fy ^ 2, expB2u4, sinhGu4) + *(-2, Fx, Fy, coshGu4, exp(*(B1, u4) + *(B2, u4)))) + *(B1 ^ 2, *(Fx ^ 2, exp(*(B2, u4) + *(2, B1, u4))) + *(Fy ^ 2, expB2u4), coshGu4)) + *(-2, Fxpt, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-2, Fyph, coshGu4, expB2u4) + *(2, Fxph, exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(2, Fypt, exp(*(B1, u4) + *(B2, u4)), sinhGu4))

    nothing
end


function B2d_eq_coeff!(ABCS::Vector, vars::BdGVars, ::Inner)
    @unpack (
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

    ABCS[2] = *(-4/27, u2, (3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi))) ^ 4, expB1u4)

    ABCS[3] = *(-2/9, u, (3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi))) ^ 3, 3 + *(-3, Sp, u3) + *(6, u, xi) + *(15, S, u4) + *(phi0 ^ 2, u2, -3 + *(4, u, xi)), expB1u4)

    ABCS[4] = *(4/9, u5, *((*(Fx, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(u, *(3, St) + *(2, xi, xi_x, phi0 ^ 2) + *(9, S, u, xi_x))) ^ 2, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *((*(Fy, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(u, *(3, Sh) + *(2, xi, xi_y, phi0 ^ 2) + *(9, S, u, xi_y))) ^ 2, coshGu4, expB2u4) + *(-2, *(Fx, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(u, *(3, St) + *(2, xi, xi_x, phi0 ^ 2) + *(9, S, u, xi_x)), *(Fy, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(u, *(3, Sh) + *(2, xi, xi_y, phi0 ^ 2) + *(9, S, u, xi_y)), exp(*(B1, u4) + *(B2, u4)), sinhGu4)) + *(-2/9, u2, 3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi)), *(u4, *(Gh + *(4, G, u, xi_y + *(Fy, u2)), *(Fy, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(u, *(3, Sh) + *(2, xi, xi_y, phi0 ^ 2) + *(9, S, u, xi_y)), expB2u4) + *(Gt + *(4, G, u, xi_x + *(Fx, u2)), *(Fx, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(u, *(3, St) + *(2, xi, xi_x, phi0 ^ 2) + *(9, S, u, xi_x)), exp(*(B2, u4) + *(2, B1, u4))), sinhGu4) + *(*(Fxt, -3 + *(-3, Sp, u3) + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(u, *(3, Sb) + *(Fx, *(-2, xi_x, 3 + *(-27, S, u4) + *(-12, B2, u4) + *(3, Sp, u3) + *(6, B1, u4) + *(-36, B1, S, u8) + *(72, B2, S, u8)) + *(2, B2t, *(3, u3) + *(phi0 ^ 2, u5) + *(-9, S, u ^ 7) + *(-2, xi, phi0 ^ 2, u6)) + *(6, St, u3) + *(B1t, u3, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(2, Fxp, u, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(2, St, u ^ 7, *(-12, B2) + *(6, B1)) + *(2, xi_x, phi0 ^ 2, u2, -1 + *(4, u, xi), 1 + *(-4, B2, u4) + *(2, B1, u4))) + *(-3, Sp, xi_xx) + *(2, phi0 ^ 2, xi_x ^ 2) + *(-6, B2t, St, u4) + *(2, xi, xi_xx, phi0 ^ 2) + *(2, Fx ^ 2, u2, 3 + *(-3, Sp, u3) + *(9, S, u4) + *(12, B2, u4) + *(xi, phi0 ^ 2, u3) + *(-36, B2, S, u8) + *(2, B1, u4, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(4, B2, phi0 ^ 2, u6) + *(-8, B2, xi, phi0 ^ 2, u ^ 7)) + *(3, B1t, St, u4) + *(6, Fxp, St, u2) + *(9, S, u, xi_xx) + *(18, St, u, xi_x) + *(36, S, u2, xi_x ^ 2) + *(-72, B2, S, u6, xi_x ^ 2) + *(-24, B2, St, xi_x, u5) + *(-18, B2t, S, xi_x, u5) + *(6, u, xi, phi0 ^ 2, xi_x ^ 2) + *(9, B1t, S, xi_x, u5) + *(12, B1, St, xi_x, u5) + *(18, Fxp, S, xi_x, u3) + *(36, B1, S, u6, xi_x ^ 2) + *(-16, B2, xi, phi0 ^ 2, u5, xi_x ^ 2) + *(-4, B2t, xi, xi_x, phi0 ^ 2, u4) + *(2, B1t, xi, xi_x, phi0 ^ 2, u4) + *(4, Fxp, xi, xi_x, phi0 ^ 2, u2) + *(8, B1, xi, phi0 ^ 2, u5, xi_x ^ 2)), coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(*(Fyh, -3 + *(-3, Sp, u3) + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(-1, u, *(-3, Ss) + *(Fy, *(2, xi_y, 3 + *(-27, S, u4) + *(-12, B2, u4) + *(-6, B1, u4) + *(3, Sp, u3) + *(36, B1, S, u8) + *(72, B2, S, u8)) + *(B1h, u3, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(2, B2h, u3, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(2, Fyp, u, 3 + *(phi0 ^ 2, u2 + *(-2, xi, u3)) + *(-9, S, u4)) + *(6, Sh, u3, -1 + *(2, B1, u4) + *(4, B2, u4)) + *(2, xi_y, phi0 ^ 2, u2, -1 + *(4, u, xi), -1 + *(2, B1, u4) + *(4, B2, u4))) + *(-2, phi0 ^ 2, xi_y ^ 2) + *(3, Sp, xi_yy) + *(-36, S, u2, xi_y ^ 2) + *(-18, Sh, u, xi_y) + *(-9, S, u, xi_yy) + *(-6, Fyp, Sh, u2) + *(-2, xi, xi_yy, phi0 ^ 2) + *(2, Fy ^ 2, u2, -3 + *(-12, B2, u4) + *(-9, S, u4) + *(3, Sp, u3) + *(-1, xi, phi0 ^ 2, u3) + *(-4, B2, phi0 ^ 2, u6) + *(2, B1, u4, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(36, B2, S, u8) + *(8, B2, xi, phi0 ^ 2, u ^ 7)) + *(3, B1h, Sh, u4) + *(6, B2h, Sh, u4) + *(-18, Fyp, S, xi_y, u3) + *(-6, u, xi, phi0 ^ 2, xi_y ^ 2) + *(9, B1h, S, xi_y, u5) + *(12, B1, Sh, xi_y, u5) + *(18, B2h, S, xi_y, u5) + *(24, B2, Sh, xi_y, u5) + *(36, B1, S, u6, xi_y ^ 2) + *(72, B2, S, u6, xi_y ^ 2) + *(-4, Fyp, xi, xi_y, phi0 ^ 2, u2) + *(2, B1h, xi, xi_y, phi0 ^ 2, u4) + *(4, B2h, xi, xi_y, phi0 ^ 2, u4) + *(8, B1, xi, phi0 ^ 2, u5, xi_y ^ 2) + *(16, B2, xi, phi0 ^ 2, u5, xi_y ^ 2)), coshGu4, expB2u4) + *(*(Fxh, 3 + *(phi0 ^ 2, u2 + *(-2, xi, u3)) + *(-9, S, u4) + *(3, Sp, u3)) + *(Fyt, 3 + *(phi0 ^ 2, u2 + *(-2, xi, u3)) + *(-9, S, u4) + *(3, Sp, u3)) + *(2, u, *(-3, Sc) + *(3, Fx, xi_y) + *(3, Fy, xi_x) + *(3, Sp, xi_xy) + *(-1, Fxp, u, *(Fy, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(u, *(3, Sh) + *(2, xi, xi_y, phi0 ^ 2) + *(9, S, u, xi_y))) + *(-9, S, u, xi_xy) + *(-9, Sh, u, xi_x) + *(-9, St, u, xi_y) + *(-6, Fx, Fy, u2) + *(-3, B2h, Fx, u3) + *(-3, B2t, Fy, u3) + *(-3, Fx, Sh, u3) + *(-3, Fy, St, u3) + *(-3, Fyp, St, u2) + *(-2, xi, xi_xy, phi0 ^ 2) + *(-2, xi_x, xi_y, phi0 ^ 2) + *(3, B2h, St, u4) + *(3, B2t, Sh, u4) + *(3, Fx, Fyp, u) + *(Fx, Fyp, phi0 ^ 2, u3) + *(Fx, xi_y, phi0 ^ 2, u2) + *(Fy, xi_x, phi0 ^ 2, u2) + *(-1, B2h, Fx, phi0 ^ 2, u5) + *(-1, B2t, Fy, phi0 ^ 2, u5) + *(-36, S, xi_x, xi_y, u2) + *(-27, Fx, S, xi_y, u4) + *(-27, Fy, S, xi_x, u4) + *(-24, B2, Fx, Fy, u6) + *(-18, Fx, Fy, S, u6) + *(-12, B2, Fx, xi_y, u4) + *(-12, B2, Fy, xi_x, u4) + *(-9, Fx, Fyp, S, u5) + *(-9, Fyp, S, xi_x, u3) + *(3, Fx, Sp, xi_y, u3) + *(3, Fy, Sp, xi_x, u3) + *(6, Fx, Fy, Sp, u5) + *(9, B2h, Fx, S, u ^ 7) + *(9, B2h, S, xi_x, u5) + *(9, B2t, Fy, S, u ^ 7) + *(9, B2t, S, xi_y, u5) + *(12, B2, Fx, Sh, u ^ 7) + *(12, B2, Fy, St, u ^ 7) + *(12, B2, Sh, xi_x, u5) + *(12, B2, St, xi_y, u5) + *(-8, B2, Fx, Fy, phi0 ^ 2, u8) + *(-6, u, xi, xi_x, xi_y, phi0 ^ 2) + *(-4, B2, Fx, xi_y, phi0 ^ 2, u6) + *(-4, B2, Fy, xi_x, phi0 ^ 2, u6) + *(-4, Fx, xi, xi_y, phi0 ^ 2, u3) + *(-4, Fy, xi, xi_x, phi0 ^ 2, u3) + *(-2, Fx, Fy, xi, phi0 ^ 2, u5) + *(-2, Fx, Fyp, xi, phi0 ^ 2, u4) + *(-2, Fyp, xi, xi_x, phi0 ^ 2, u2) + *(2, B2h, Fx, xi, phi0 ^ 2, u6) + *(2, B2h, xi, xi_x, phi0 ^ 2, u4) + *(2, B2t, Fy, xi, phi0 ^ 2, u6) + *(2, B2t, xi, xi_y, phi0 ^ 2, u4) + *(72, B2, Fx, Fy, S, u ^ 10) + *(72, B2, Fx, S, xi_y, u8) + *(72, B2, Fy, S, xi_x, u8) + *(72, B2, S, xi_x, xi_y, u6) + *(16, B2, Fx, Fy, xi, phi0 ^ 2, u ^ 9) + *(16, B2, Fx, xi, xi_y, phi0 ^ 2, u ^ 7) + *(16, B2, Fy, xi, xi_x, phi0 ^ 2, u ^ 7) + *(16, B2, xi, xi_x, xi_y, phi0 ^ 2, u5)), exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-1, u4, *(Fx, Gh, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(Fy, Gt, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(Gh, u, *(3, St) + *(2, xi, xi_x, phi0 ^ 2) + *(9, S, u, xi_x)) + *(Gt, u, *(3, Sh) + *(2, xi, xi_y, phi0 ^ 2) + *(9, S, u, xi_y)) + *(4, G, u2, *(xi_y, *(3, St) + *(4, xi, xi_x, phi0 ^ 2) + *(18, S, u, xi_x)) + *(3, Sh, xi_x)) + *(4, Fx, G, u, *(xi_y, -3 + *(18, S, u4) + *(phi0 ^ 2, u2, -1 + *(4, u, xi))) + *(3, Sh, u3) + *(2, Fy, u2, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi)))) + *(4, Fy, G, u, *(xi_x, -3 + *(18, S, u4) + *(phi0 ^ 2, u2, -1 + *(4, u, xi))) + *(3, St, u3)), coshGu4, exp(*(B1, u4) + *(B2, u4)))) + *(-1/9, u, (3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi))) ^ 2, *(*(-2, u2, *(*(-1, Gs) + *(Gp, xi_yy) + *(Fyh, Gp, u2) + *(Fyp, u2, Gh + *(4, G, u, xi_y + *(Fy, u2))) + *(-36, G, Fy ^ 2, u6) + *(-20, G, u2, xi_y ^ 2) + *(-10, Fy, Gh, u3) + *(-8, Gh, u, xi_y) + *(-4, Fyh, G, u3) + *(-4, G, u, xi_yy) + *(2, B1h, Gh, u4) + *(2, B2h, Gh, u4) + *(2, Gp, Fy ^ 2, u5) + *(-56, Fy, G, xi_y, u4) + *(2, Fy, Gp, xi_y, u3) + *(8, B1, Fy, Gh, u ^ 7) + *(8, B1, Gh, xi_y, u5) + *(8, B1h, Fy, G, u ^ 7) + *(8, B1h, G, xi_y, u5) + *(8, B2, Fy, Gh, u ^ 7) + *(8, B2, Gh, xi_y, u5) + *(8, B2h, Fy, G, u ^ 7) + *(8, B2h, G, xi_y, u5) + *(32, B1, G, Fy ^ 2, u ^ 10) + *(32, B1, G, u6, xi_y ^ 2) + *(32, B2, G, Fy ^ 2, u ^ 10) + *(32, B2, G, u6, xi_y ^ 2) + *(64, B1, Fy, G, xi_y, u8) + *(64, B2, Fy, G, xi_y, u8), expB2u4) + *(*(-1, Gb) + *(Gp, xi_xx) + *(Fxp, Gt, u2) + *(Fxt, u2, Gp + *(-4, G, u)) + *(-36, G, Fx ^ 2, u6) + *(-20, G, u2, xi_x ^ 2) + *(-10, Fx, Gt, u3) + *(-8, Gt, u, xi_x) + *(-4, G, u, xi_xx) + *(-2, B1t, Gt, u4) + *(2, B2t, Gt, u4) + *(2, Gp, Fx ^ 2, u5) + *(-56, Fx, G, xi_x, u4) + *(-32, B1, G, Fx ^ 2, u ^ 10) + *(-32, B1, G, u6, xi_x ^ 2) + *(-8, B1, Fx, Gt, u ^ 7) + *(-8, B1, Gt, xi_x, u5) + *(-8, B1t, Fx, G, u ^ 7) + *(-8, B1t, G, xi_x, u5) + *(2, Fx, Gp, xi_x, u3) + *(4, Fx, Fxp, G, u5) + *(4, Fxp, G, xi_x, u3) + *(8, B2, Fx, Gt, u ^ 7) + *(8, B2, Gt, xi_x, u5) + *(8, B2t, Fx, G, u ^ 7) + *(8, B2t, G, xi_x, u5) + *(32, B2, G, Fx ^ 2, u ^ 10) + *(32, B2, G, u6, xi_x ^ 2) + *(-64, B1, Fx, G, xi_x, u8) + *(64, B2, Fx, G, xi_x, u8), exp(*(B2, u4) + *(2, B1, u4)))) + *(2, Fxph + Fypt + *(u, *(-2, Fxh, 1 + *(B2p, u3) + *(-4, B2, u4)) + *(-2, Fyt, 1 + *(B2p, u3) + *(-4, B2, u4)) + *(2, Fxp, xi_y) + *(2, Fyp, xi_x) + *(4, B2c, u) + *(-1, Fxp, Fyp, u) + *(-16, Fx, Fy, u3) + *(-6, Fx, u, xi_y) + *(-6, Fy, u, xi_x) + *(-4, B2p, u, xi_xy) + *(-2, B2h, Fxp, u3) + *(-2, B2t, Fyp, u3) + *(-2, Gh, Gt, u5) + *(2, B2h, B2t, u5) + *(4, Fx, Fyp, u2) + *(4, Fxp, Fy, u2) + *(16, B2, xi_xy, u2) + *(16, B2h, xi_x, u2) + *(16, B2t, xi_y, u2) + *(20, B2h, Fx, u4) + *(20, B2t, Fy, u4) + *(-32, Fx, Fy, G ^ 2, u ^ 11) + *(-32, Fx, xi_y, G ^ 2, u ^ 9) + *(-32, Fy, xi_x, G ^ 2, u ^ 9) + *(-32, xi_x, xi_y, G ^ 2, u ^ 7) + *(-8, B2, Fx, Fyp, u6) + *(-8, B2, Fxp, Fy, u6) + *(-8, B2, Fxp, xi_y, u4) + *(-8, B2, Fyp, xi_x, u4) + *(-8, B2p, Fx, Fy, u6) + *(-8, Fx, G, Gh, u8) + *(-8, Fy, G, Gt, u8) + *(-8, G, Gh, xi_x, u6) + *(-8, G, Gt, xi_y, u6) + *(-4, B2p, Fx, xi_y, u4) + *(-4, B2p, Fy, xi_x, u4) + *(-4, Fx, Fy, phi0 ^ 2, u5) + *(-4, Fx, phih, phi0 ^ 4, u4) + *(-4, Fy, phit, phi0 ^ 4, u4) + *(-4, phih, phit, phi0 ^ 6, u3) + *(8, B2, B2h, Fx, u8) + *(8, B2, B2h, xi_x, u6) + *(8, B2, B2t, Fy, u8) + *(8, B2, B2t, xi_y, u6) + *(32, Fx, Fy, B2 ^ 2, u ^ 11) + *(32, Fx, xi_y, B2 ^ 2, u ^ 9) + *(32, Fy, xi_x, B2 ^ 2, u ^ 9) + *(32, xi_x, xi_y, B2 ^ 2, u ^ 7) + *(80, B2, xi_x, xi_y, u3) + *(112, B2, Fx, xi_y, u5) + *(112, B2, Fy, xi_x, u5) + *(144, B2, Fx, Fy, u ^ 7) + *(-36, Fx, Fy, phi ^ 2, phi0 ^ 6, u ^ 9) + *(-36, Fx, xi_y, phi ^ 2, phi0 ^ 6, u ^ 7) + *(-36, Fy, xi_x, phi ^ 2, phi0 ^ 6, u ^ 7) + *(-36, xi_x, xi_y, phi ^ 2, phi0 ^ 6, u5) + *(-24, Fx, Fy, phi, phi0 ^ 4, u ^ 7) + *(-16, Fx, Fy, phi0 ^ 2, u ^ 7, xi ^ 2) + *(-16, Fx, xi_y, phi0 ^ 2, u5, xi ^ 2) + *(-16, Fy, xi_x, phi0 ^ 2, u5, xi ^ 2) + *(-16, xi_x, xi_y, phi0 ^ 2, u3, xi ^ 2) + *(-12, Fx, phi, phih, phi0 ^ 6, u6) + *(-12, Fx, phi, xi_y, phi0 ^ 4, u5) + *(-12, Fy, phi, phit, phi0 ^ 6, u6) + *(-12, Fy, phi, xi_x, phi0 ^ 4, u5) + *(-12, phi, phih, xi_x, phi0 ^ 6, u4) + *(-12, phi, phit, xi_y, phi0 ^ 6, u4) + *(8, Fx, phih, xi, phi0 ^ 4, u5) + *(8, Fx, xi, xi_y, phi0 ^ 2, u4) + *(8, Fy, phit, xi, phi0 ^ 4, u5) + *(8, Fy, xi, xi_x, phi0 ^ 2, u4) + *(8, phih, xi, xi_x, phi0 ^ 4, u3) + *(8, phit, xi, xi_y, phi0 ^ 4, u3) + *(16, Fx, Fy, xi, phi0 ^ 2, u6) + *(48, Fx, Fy, phi, xi, phi0 ^ 4, u8) + *(48, Fx, phi, xi, xi_y, phi0 ^ 4, u6) + *(48, Fy, phi, xi, xi_x, phi0 ^ 4, u6) + *(48, phi, xi, xi_x, xi_y, phi0 ^ 4, u4)), exp(*(B1, u4) + *(B2, u4))), sinhGu4) + *(*(*(-2, Fxpt) + *(u, *(Fxt, 4 + *(-16, B2, u4) + *(-2, B1p, u3) + *(4, B2p, u3) + *(8, B1, u4)) + *(u, Fxp ^ 2) + *(-4, B2b, u) + *(-4, Fxp, xi_x) + *(-2, B2t ^ 2, u5) + *(2, B1b, u) + *(2, B1t ^ 2, u5) + *(2, Gt ^ 2, u5) + *(16, Fx ^ 2, u3) + *(-144, B2, Fx ^ 2, u ^ 7) + *(-80, B2, u3, xi_x ^ 2) + *(-40, B2t, Fx, u4) + *(-32, B2t, xi_x, u2) + *(-32, B2 ^ 2, Fx ^ 2, u ^ 11) + *(-32, B2 ^ 2, u ^ 7, xi_x ^ 2) + *(-16, B2, xi_xx, u2) + *(-8, Fx, Fxp, u2) + *(-4, B1p, Fx ^ 2, u6) + *(-4, B1t, B2t, u5) + *(-2, B1p, u, xi_xx) + *(-2, B1t, Fxp, u3) + *(4, B2p, u, xi_xx) + *(4, B2t, Fxp, u3) + *(4, Fx ^ 2, phi0 ^ 2, u5) + *(4, phi0 ^ 6, phit ^ 2, u3) + *(8, B1, xi_xx, u2) + *(8, B2p, Fx ^ 2, u6) + *(12, Fx, u, xi_x) + *(16, B1t, xi_x, u2) + *(20, B1t, Fx, u4) + *(32, B1 ^ 2, Fx ^ 2, u ^ 11) + *(32, B1 ^ 2, u ^ 7, xi_x ^ 2) + *(32, Fx ^ 2, G ^ 2, u ^ 11) + *(32, G ^ 2, u ^ 7, xi_x ^ 2) + *(40, B1, u3, xi_x ^ 2) + *(72, B1, Fx ^ 2, u ^ 7) + *(-224, B2, Fx, xi_x, u5) + *(-64, B1, B2, Fx ^ 2, u ^ 11) + *(-64, B1, B2, u ^ 7, xi_x ^ 2) + *(-64, Fx, xi_x, B2 ^ 2, u ^ 9) + *(-16, B1, B2t, Fx, u8) + *(-16, B1, B2t, xi_x, u6) + *(-16, B1t, B2, Fx, u8) + *(-16, B1t, B2, xi_x, u6) + *(-16, B2, B2t, Fx, u8) + *(-16, B2, B2t, xi_x, u6) + *(-16, xi, Fx ^ 2, phi0 ^ 2, u6) + *(-8, B1, Fx, Fxp, u6) + *(-8, B1, Fxp, xi_x, u4) + *(-4, B1p, Fx, xi_x, u4) + *(8, B2p, Fx, xi_x, u4) + *(8, Fx, phit, phi0 ^ 4, u4) + *(16, B1, B1t, Fx, u8) + *(16, B1, B1t, xi_x, u6) + *(16, B2, Fx, Fxp, u6) + *(16, B2, Fxp, xi_x, u4) + *(16, Fx, G, Gt, u8) + *(16, G, Gt, xi_x, u6) + *(16, Fx ^ 2, phi0 ^ 2, u ^ 7, xi ^ 2) + *(16, phi0 ^ 2, u3, xi ^ 2, xi_x ^ 2) + *(24, phi, Fx ^ 2, phi0 ^ 4, u ^ 7) + *(36, Fx ^ 2, phi ^ 2, phi0 ^ 6, u ^ 9) + *(36, phi ^ 2, phi0 ^ 6, u5, xi_x ^ 2) + *(64, Fx, xi_x, B1 ^ 2, u ^ 9) + *(64, Fx, xi_x, G ^ 2, u ^ 9) + *(112, B1, Fx, xi_x, u5) + *(-128, B1, B2, Fx, xi_x, u ^ 9) + *(-48, phi, xi, Fx ^ 2, phi0 ^ 4, u8) + *(-48, phi, xi, phi0 ^ 4, u4, xi_x ^ 2) + *(-16, Fx, phit, xi, phi0 ^ 4, u5) + *(-16, Fx, xi, xi_x, phi0 ^ 2, u4) + *(-16, phit, xi, xi_x, phi0 ^ 4, u3) + *(24, Fx, phi, phit, phi0 ^ 6, u6) + *(24, Fx, phi, xi_x, phi0 ^ 4, u5) + *(24, phi, phit, xi_x, phi0 ^ 6, u4) + *(32, Fx, xi_x, phi0 ^ 2, u5, xi ^ 2) + *(72, Fx, xi_x, phi ^ 2, phi0 ^ 6, u ^ 7) + *(-96, Fx, phi, xi, xi_x, phi0 ^ 4, u6)), exp(*(B2, u4) + *(2, B1, u4))) + *(-2, Fyph, expB2u4) + *(u, *(u, Fyp ^ 2) + *(-4, B2s, u) + *(-4, Fyp, xi_y) + *(-2, B1s, u) + *(-2, B2h ^ 2, u5) + *(2, Fyh, 2 + *(B1p, u3) + *(-8, B2, u4) + *(-4, B1, u4) + *(2, B2p, u3)) + *(2, B1h ^ 2, u5) + *(2, Gh ^ 2, u5) + *(16, Fy ^ 2, u3) + *(-144, B2, Fy ^ 2, u ^ 7) + *(-80, B2, u3, xi_y ^ 2) + *(-72, B1, Fy ^ 2, u ^ 7) + *(-40, B1, u3, xi_y ^ 2) + *(-40, B2h, Fy, u4) + *(-32, B2h, xi_y, u2) + *(-32, B2 ^ 2, Fy ^ 2, u ^ 11) + *(-32, B2 ^ 2, u ^ 7, xi_y ^ 2) + *(-20, B1h, Fy, u4) + *(-16, B1h, xi_y, u2) + *(-16, B2, xi_yy, u2) + *(-8, B1, xi_yy, u2) + *(-8, Fy, Fyp, u2) + *(2, B1h, Fyp, u3) + *(2, B1p, u, xi_yy) + *(4, B1h, B2h, u5) + *(4, B1p, Fy ^ 2, u6) + *(4, B2h, Fyp, u3) + *(4, B2p, u, xi_yy) + *(4, Fy ^ 2, phi0 ^ 2, u5) + *(4, phi0 ^ 6, phih ^ 2, u3) + *(8, B2p, Fy ^ 2, u6) + *(12, Fy, u, xi_y) + *(32, B1 ^ 2, Fy ^ 2, u ^ 11) + *(32, B1 ^ 2, u ^ 7, xi_y ^ 2) + *(32, Fy ^ 2, G ^ 2, u ^ 11) + *(32, G ^ 2, u ^ 7, xi_y ^ 2) + *(-224, B2, Fy, xi_y, u5) + *(-112, B1, Fy, xi_y, u5) + *(-64, Fy, xi_y, B2 ^ 2, u ^ 9) + *(-16, B2, B2h, Fy, u8) + *(-16, B2, B2h, xi_y, u6) + *(-16, xi, Fy ^ 2, phi0 ^ 2, u6) + *(4, B1p, Fy, xi_y, u4) + *(8, B1, Fy, Fyp, u6) + *(8, B1, Fyp, xi_y, u4) + *(8, B2p, Fy, xi_y, u4) + *(8, Fy, phih, phi0 ^ 4, u4) + *(16, B1, B1h, Fy, u8) + *(16, B1, B1h, xi_y, u6) + *(16, B1, B2h, Fy, u8) + *(16, B1, B2h, xi_y, u6) + *(16, B1h, B2, Fy, u8) + *(16, B1h, B2, xi_y, u6) + *(16, B2, Fy, Fyp, u6) + *(16, B2, Fyp, xi_y, u4) + *(16, Fy, G, Gh, u8) + *(16, G, Gh, xi_y, u6) + *(16, Fy ^ 2, phi0 ^ 2, u ^ 7, xi ^ 2) + *(16, phi0 ^ 2, u3, xi ^ 2, xi_y ^ 2) + *(24, phi, Fy ^ 2, phi0 ^ 4, u ^ 7) + *(36, Fy ^ 2, phi ^ 2, phi0 ^ 6, u ^ 9) + *(36, phi ^ 2, phi0 ^ 6, u5, xi_y ^ 2) + *(64, B1, B2, Fy ^ 2, u ^ 11) + *(64, B1, B2, u ^ 7, xi_y ^ 2) + *(64, Fy, xi_y, B1 ^ 2, u ^ 9) + *(64, Fy, xi_y, G ^ 2, u ^ 9) + *(-48, phi, xi, Fy ^ 2, phi0 ^ 4, u8) + *(-48, phi, xi, phi0 ^ 4, u4, xi_y ^ 2) + *(-16, Fy, phih, xi, phi0 ^ 4, u5) + *(-16, Fy, xi, xi_y, phi0 ^ 2, u4) + *(-16, phih, xi, xi_y, phi0 ^ 4, u3) + *(24, Fy, phi, phih, phi0 ^ 6, u6) + *(24, Fy, phi, xi_y, phi0 ^ 4, u5) + *(24, phi, phih, xi_y, phi0 ^ 6, u4) + *(32, Fy, xi_y, phi0 ^ 2, u5, xi ^ 2) + *(72, Fy, xi_y, phi ^ 2, phi0 ^ 6, u ^ 7) + *(128, B1, B2, Fy, xi_y, u ^ 9) + *(-96, Fy, phi, xi, xi_y, phi0 ^ 4, u6), expB2u4) + *(2, u2, *(-2, Gc) + *(2, Gp, xi_xy) + *(B1h, Gt, u4) + *(Fxh, Gp, u2) + *(Fxp, u2, Gh + *(4, G, u, xi_y + *(Fy, u2))) + *(Fyp, Gt, u2) + *(Fyt, Gp, u2) + *(-1, B1t, Gh, u4) + *(-10, Fx, Gh, u3) + *(-10, Fy, Gt, u3) + *(-8, G, u, xi_xy) + *(-8, Gh, u, xi_x) + *(-8, Gt, u, xi_y) + *(-4, Fxh, G, u3) + *(-4, Fyt, G, u3) + *(2, B2h, Gt, u4) + *(2, B2t, Gh, u4) + *(-72, Fx, Fy, G, u6) + *(-56, Fx, G, xi_y, u4) + *(-56, Fy, G, xi_x, u4) + *(-40, G, xi_x, xi_y, u2) + *(-4, B1, Fx, Gh, u ^ 7) + *(-4, B1, Gh, xi_x, u5) + *(-4, B1t, Fy, G, u ^ 7) + *(-4, B1t, G, xi_y, u5) + *(2, Fx, Gp, xi_y, u3) + *(2, Fy, Gp, xi_x, u3) + *(4, B1, Fy, Gt, u ^ 7) + *(4, B1, Gt, xi_y, u5) + *(4, B1h, Fx, G, u ^ 7) + *(4, B1h, G, xi_x, u5) + *(4, Fx, Fy, Gp, u5) + *(4, Fx, Fyp, G, u5) + *(4, Fyp, G, xi_x, u3) + *(8, B2, Fx, Gh, u ^ 7) + *(8, B2, Fy, Gt, u ^ 7) + *(8, B2, Gh, xi_x, u5) + *(8, B2, Gt, xi_y, u5) + *(8, B2h, Fx, G, u ^ 7) + *(8, B2h, G, xi_x, u5) + *(8, B2t, Fy, G, u ^ 7) + *(8, B2t, G, xi_y, u5) + *(64, B2, Fx, Fy, G, u ^ 10) + *(64, B2, Fx, G, xi_y, u8) + *(64, B2, Fy, G, xi_x, u8) + *(64, B2, G, xi_x, xi_y, u6), exp(*(B1, u4) + *(B2, u4))), coshGu4)) + *(1/9, (3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi))) ^ 3, B2p + *(-4, B2, u), *(3, (1 + *(u, xi)) ^ 2) + *(-1, phi0 ^ 2, u2) + *(6, Sd, u4), expB1u4)

    nothing
end


# this is another coupled equation, for B1d and Gd. the notation used is
#
# ( A11 d_uu B1d + A12 d_uu Gd + B11 d_u B1d + B12 d_u Gd + C11 B1d + C12 Gd ) = -S1
# ( A21 d_uu B1d + A22 d_uu Gd + B21 d_u B1d + B22 d_u Gd + C21 B1d + C22 Gd ) = -S2

function B1dGd_eq_coeff!(AA::Matrix, BB::Matrix, CC::Matrix, SS::Vector, vars::BdGVars, ::Inner)
    @unpack (
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

    u2 = u*u
    u3 = u*u2
    u4 = u2*u2
    u5 = u4*u
    u5 = u3*u2
    u6 = u3*u3
    u8 = u4*u4

    coshGu4 = cosh(*(G, u4))
    sinhGu4 = sinh(*(G, u4))

    tanhGu4 = tanh(*(G, u4))
    sechGu4 = sech(*(G, u4))

    cosh2Gu4 = cosh(*(2, G, u4))
    sinh2Gu4 = sinh(*(2, G, u4))

    expB1u4 = exp(*(B1, u4))
    expB2u4 = exp(*(B2, u4))


    AA[1,1] = 0
    AA[1,2] = 0


    BB[1,1] = *(-4/27, u2, (3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi))) ^ 4, expB1u4)
    BB[1,2] = 0


    CC[1,1] = *(-2/27, u, (3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi))) ^ 3, 9 + *(-9, Sp, u3) + *(18, u, xi) + *(45, S, u4) + *(3, phi0 ^ 2, u2, -3 + *(4, u, xi)) + *(-2, u3, Gp + *(-4, G, u), 3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi)), tanhGu4), expB1u4)

    CC[1,2] = *(4/27, (*(3, u, 1 + *(S, u4) + *(u, xi)) + *(phi0 ^ 2, u3, -1 + *(u, xi))) ^ 4, B1p + *(-4, B1, u), expB1u4, tanhGu4)


    SS[1] = *(1/9, (3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi))) ^ 3, B1p + *(-4, B1, u), *(3, (1 + *(u, xi)) ^ 2) + *(-1, phi0 ^ 2, u2) + *(6, Sd, u4), expB1u4) + *(4/3, u5, *(-1, (*(Fy, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(u, *(3, Sh) + *(2, xi, xi_y, phi0 ^ 2) + *(9, S, u, xi_y))) ^ 2) + *((*(Fx, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(u, *(3, St) + *(2, xi, xi_x, phi0 ^ 2) + *(9, S, u, xi_x))) ^ 2, exp(*(2, B1, u4))), expB2u4, sechGu4) + *(-2/3, u2, 3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi)), *(Fyh, 3 + *(phi0 ^ 2, u2 + *(-2, xi, u3)) + *(-9, S, u4) + *(3, Sp, u3)) + *(*(Fxt, -3 + *(-3, Sp, u3) + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(u, *(3, Sb) + *(Fx, *(-1, xi_x, 6 + *(-54, S, u4) + *(6, Sp, u3) + *(12, B2, u4) + *(-72, B2, S, u8)) + *(6, St, u3 + *(2, B2, u ^ 7)) + *(B2t, u3, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(2, Fxp, u, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(2, xi_x, phi0 ^ 2, u2, 1 + *(2, B2, u4), -1 + *(4, u, xi))) + *(-3, Sp, xi_xx) + *(2, phi0 ^ 2, xi_x ^ 2) + *(2, xi, xi_xx, phi0 ^ 2) + *(2, Fx ^ 2, u2, 3 + *(-3, Sp, u3) + *(9, S, u4) + *(xi, phi0 ^ 2, u3) + *(2, B2, u4, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi)))) + *(3, B2t, St, u4) + *(6, Fxp, St, u2) + *(9, S, u, xi_xx) + *(18, St, u, xi_x) + *(36, S, u2, xi_x ^ 2) + *(6, u, xi, phi0 ^ 2, xi_x ^ 2) + *(9, B2t, S, xi_x, u5) + *(12, B2, St, xi_x, u5) + *(18, Fxp, S, xi_x, u3) + *(36, B2, S, u6, xi_x ^ 2) + *(2, B2t, xi, xi_x, phi0 ^ 2, u4) + *(4, Fxp, xi, xi_x, phi0 ^ 2, u2) + *(8, B2, xi, phi0 ^ 2, u5, xi_x ^ 2)), exp(*(2, B1, u4))) + *(-1, u, *(3, Ss) + *(Fy, *(-1, xi_y, 6 + *(-54, S, u4) + *(6, Sp, u3) + *(12, B2, u4) + *(-72, B2, S, u8)) + *(6, Sh, u3 + *(2, B2, u ^ 7)) + *(B2h, u3, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(2, Fyp, u, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(2, xi_y, phi0 ^ 2, u2, 1 + *(2, B2, u4), -1 + *(4, u, xi))) + *(-3, Sp, xi_yy) + *(2, phi0 ^ 2, xi_y ^ 2) + *(2, xi, xi_yy, phi0 ^ 2) + *(2, Fy ^ 2, u2, 3 + *(-3, Sp, u3) + *(9, S, u4) + *(xi, phi0 ^ 2, u3) + *(2, B2, u4, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi)))) + *(3, B2h, Sh, u4) + *(6, Fyp, Sh, u2) + *(9, S, u, xi_yy) + *(18, Sh, u, xi_y) + *(36, S, u2, xi_y ^ 2) + *(6, u, xi, phi0 ^ 2, xi_y ^ 2) + *(9, B2h, S, xi_y, u5) + *(12, B2, Sh, xi_y, u5) + *(18, Fyp, S, xi_y, u3) + *(36, B2, S, u6, xi_y ^ 2) + *(2, B2h, xi, xi_y, phi0 ^ 2, u4) + *(4, Fyp, xi, xi_y, phi0 ^ 2, u2) + *(8, B2, xi, phi0 ^ 2, u5, xi_y ^ 2)) + *(u4, *(Fx, Gh, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(Fy, Gt, 3 + *(phi0 ^ 2, u2 + *(-2, xi, u3)) + *(-9, S, u4)) + *(Gh, u, *(3, St) + *(2, xi, xi_x, phi0 ^ 2) + *(9, S, u, xi_x)) + *(-1, Gt, u, *(3, Sh) + *(2, xi, xi_y, phi0 ^ 2) + *(9, S, u, xi_y)) + *(12, G, u2, *(St, xi_y) + *(-1, Sh, xi_x)) + *(-4, Fx, G, u, *(xi_y, 3 + *(phi0 ^ 2, u2)) + *(3, Sh, u3)) + *(4, Fy, G, u, *(xi_x, 3 + *(phi0 ^ 2, u2)) + *(3, St, u3)), expB1u4), expB2u4, sechGu4) + *(-1/3, u, (3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi))) ^ 2, *(2, Fyph) + *(*(-2, Fxpt) + *(u, *(Fxt, 4 + *(-2, B2p, u3) + *(8, B2, u4)) + *(u, Fxp ^ 2) + *(-4, Fxp, xi_x) + *(2, B2b, u) + *(4, B2t ^ 2, u5) + *(16, Fx ^ 2, u3) + *(-8, Fx, Fxp, u2) + *(-4, B2p, Fx ^ 2, u6) + *(-2, B2p, u, xi_xx) + *(-2, B2t, Fxp, u3) + *(4, Fx ^ 2, phi0 ^ 2, u5) + *(4, phi0 ^ 6, phit ^ 2, u3) + *(8, B2, xi_xx, u2) + *(12, Fx, u, xi_x) + *(16, B2t, xi_x, u2) + *(20, B2t, Fx, u4) + *(40, B2, u3, xi_x ^ 2) + *(64, B2 ^ 2, Fx ^ 2, u ^ 11) + *(64, B2 ^ 2, u ^ 7, xi_x ^ 2) + *(72, B2, Fx ^ 2, u ^ 7) + *(-16, xi, Fx ^ 2, phi0 ^ 2, u6) + *(-8, B2, Fx, Fxp, u6) + *(-8, B2, Fxp, xi_x, u4) + *(-4, B2p, Fx, xi_x, u4) + *(8, Fx, phit, phi0 ^ 4, u4) + *(16, Fx ^ 2, phi0 ^ 2, u ^ 7, xi ^ 2) + *(16, phi0 ^ 2, u3, xi ^ 2, xi_x ^ 2) + *(24, phi, Fx ^ 2, phi0 ^ 4, u ^ 7) + *(32, B2, B2t, Fx, u8) + *(32, B2, B2t, xi_x, u6) + *(36, Fx ^ 2, phi ^ 2, phi0 ^ 6, u ^ 9) + *(36, phi ^ 2, phi0 ^ 6, u5, xi_x ^ 2) + *(112, B2, Fx, xi_x, u5) + *(128, Fx, xi_x, B2 ^ 2, u ^ 9) + *(-48, phi, xi, Fx ^ 2, phi0 ^ 4, u8) + *(-48, phi, xi, phi0 ^ 4, u4, xi_x ^ 2) + *(-16, Fx, phit, xi, phi0 ^ 4, u5) + *(-16, Fx, xi, xi_x, phi0 ^ 2, u4) + *(-16, phit, xi, xi_x, phi0 ^ 4, u3) + *(24, Fx, phi, phit, phi0 ^ 6, u6) + *(24, Fx, phi, xi_x, phi0 ^ 4, u5) + *(24, phi, phit, xi_x, phi0 ^ 6, u4) + *(32, Fx, xi_x, phi0 ^ 2, u5, xi ^ 2) + *(72, Fx, xi_x, phi ^ 2, phi0 ^ 6, u ^ 7) + *(-96, Fx, phi, xi, xi_x, phi0 ^ 4, u6)), exp(*(2, B1, u4))) + *(-1, u, *(Fyh, 4 + *(-2, B2p, u3) + *(8, B2, u4)) + *(u, Fyp ^ 2) + *(-4, Fyp, xi_y) + *(2, B2s, u) + *(4, B2h ^ 2, u5) + *(16, Fy ^ 2, u3) + *(-8, Fy, Fyp, u2) + *(-4, B2p, Fy ^ 2, u6) + *(-2, B2h, Fyp, u3) + *(-2, B2p, u, xi_yy) + *(4, Fy ^ 2, phi0 ^ 2, u5) + *(4, phi0 ^ 6, phih ^ 2, u3) + *(8, B2, xi_yy, u2) + *(12, Fy, u, xi_y) + *(16, B2h, xi_y, u2) + *(20, B2h, Fy, u4) + *(40, B2, u3, xi_y ^ 2) + *(64, B2 ^ 2, Fy ^ 2, u ^ 11) + *(64, B2 ^ 2, u ^ 7, xi_y ^ 2) + *(72, B2, Fy ^ 2, u ^ 7) + *(-16, xi, Fy ^ 2, phi0 ^ 2, u6) + *(-8, B2, Fy, Fyp, u6) + *(-8, B2, Fyp, xi_y, u4) + *(-4, B2p, Fy, xi_y, u4) + *(8, Fy, phih, phi0 ^ 4, u4) + *(16, Fy ^ 2, phi0 ^ 2, u ^ 7, xi ^ 2) + *(16, phi0 ^ 2, u3, xi ^ 2, xi_y ^ 2) + *(24, phi, Fy ^ 2, phi0 ^ 4, u ^ 7) + *(32, B2, B2h, Fy, u8) + *(32, B2, B2h, xi_y, u6) + *(36, Fy ^ 2, phi ^ 2, phi0 ^ 6, u ^ 9) + *(36, phi ^ 2, phi0 ^ 6, u5, xi_y ^ 2) + *(112, B2, Fy, xi_y, u5) + *(128, Fy, xi_y, B2 ^ 2, u ^ 9) + *(-48, phi, xi, Fy ^ 2, phi0 ^ 4, u8) + *(-48, phi, xi, phi0 ^ 4, u4, xi_y ^ 2) + *(-16, Fy, phih, xi, phi0 ^ 4, u5) + *(-16, Fy, xi, xi_y, phi0 ^ 2, u4) + *(-16, phih, xi, xi_y, phi0 ^ 4, u3) + *(24, Fy, phi, phih, phi0 ^ 6, u6) + *(24, Fy, phi, xi_y, phi0 ^ 4, u5) + *(24, phi, phih, xi_y, phi0 ^ 6, u4) + *(32, Fy, xi_y, phi0 ^ 2, u5, xi ^ 2) + *(72, Fy, xi_y, phi ^ 2, phi0 ^ 6, u ^ 7) + *(-96, Fy, phi, xi, xi_y, phi0 ^ 4, u6)) + *(2, u4, *(Fxh, Gp) + *(Fyp, Gt) + *(-1, Fxp, Gh + *(4, G, u, xi_y + *(Fy, u2))) + *(-1, Fyt, Gp) + *(B2t, Gh, u2) + *(-1, B2h, Gt, u2) + *(-4, Fxh, G, u) + *(-2, Fy, Gt, u) + *(2, Fx, Gh, u) + *(4, Fyt, G, u) + *(-4, B2, Fy, Gt, u5) + *(-4, B2, Gt, xi_y, u3) + *(-4, B2h, Fx, G, u5) + *(-4, B2h, G, xi_x, u3) + *(-2, Fy, Gp, u, xi_x) + *(2, Fx, Gp, u, xi_y) + *(4, B2, Fx, Gh, u5) + *(4, B2, Gh, xi_x, u3) + *(4, B2t, Fy, G, u5) + *(4, B2t, G, xi_y, u3) + *(4, Fx, Fyp, G, u3) + *(4, Fyp, G, u, xi_x), expB1u4), expB2u4, sechGu4)


    AA[2,1] = 0
    AA[2,2] = 0

    BB[2,1] = 0
    BB[2,2] = *(-4/27, u2, (3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi))) ^ 4, expB1u4)

    CC[2,1] = *(-2/27, (*(3, u, 1 + *(S, u4) + *(u, xi)) + *(phi0 ^ 2, u3, -1 + *(u, xi))) ^ 4, B1p + *(-4, B1, u), expB1u4, sinh2Gu4)

    CC[2,2] = *(-2/9, u, (3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi))) ^ 3, 3 + *(-3, Sp, u3) + *(6, u, xi) + *(15, S, u4) + *(phi0 ^ 2, u2, -3 + *(4, u, xi)), expB1u4)

    SS[2] = *(-1/3, u, (3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi))) ^ 2, *(*(-2, Fyph) + *(u, *(Fyh, 4 + *(-2, B2p, u3) + *(8, B2, u4)) + *(u, Fyp ^ 2) + *(-4, Fyp, xi_y) + *(2, B2s, u) + *(4, B2h ^ 2, u5) + *(16, Fy ^ 2, u3) + *(-8, Fy, Fyp, u2) + *(-4, B2p, Fy ^ 2, u6) + *(-2, B2h, Fyp, u3) + *(-2, B2p, u, xi_yy) + *(4, Fy ^ 2, phi0 ^ 2, u5) + *(4, phi0 ^ 6, phih ^ 2, u3) + *(8, B2, xi_yy, u2) + *(12, Fy, u, xi_y) + *(16, B2h, xi_y, u2) + *(20, B2h, Fy, u4) + *(40, B2, u3, xi_y ^ 2) + *(64, B2 ^ 2, Fy ^ 2, u ^ 11) + *(64, B2 ^ 2, u ^ 7, xi_y ^ 2) + *(72, B2, Fy ^ 2, u ^ 7) + *(-16, xi, Fy ^ 2, phi0 ^ 2, u6) + *(-8, B2, Fy, Fyp, u6) + *(-8, B2, Fyp, xi_y, u4) + *(-4, B2p, Fy, xi_y, u4) + *(8, Fy, phih, phi0 ^ 4, u4) + *(16, Fy ^ 2, phi0 ^ 2, u ^ 7, xi ^ 2) + *(16, phi0 ^ 2, u3, xi ^ 2, xi_y ^ 2) + *(24, phi, Fy ^ 2, phi0 ^ 4, u ^ 7) + *(32, B2, B2h, Fy, u8) + *(32, B2, B2h, xi_y, u6) + *(36, Fy ^ 2, phi ^ 2, phi0 ^ 6, u ^ 9) + *(36, phi ^ 2, phi0 ^ 6, u5, xi_y ^ 2) + *(112, B2, Fy, xi_y, u5) + *(128, Fy, xi_y, B2 ^ 2, u ^ 9) + *(-48, phi, xi, Fy ^ 2, phi0 ^ 4, u8) + *(-48, phi, xi, phi0 ^ 4, u4, xi_y ^ 2) + *(-16, Fy, phih, xi, phi0 ^ 4, u5) + *(-16, Fy, xi, xi_y, phi0 ^ 2, u4) + *(-16, phih, xi, xi_y, phi0 ^ 4, u3) + *(24, Fy, phi, phih, phi0 ^ 6, u6) + *(24, Fy, phi, xi_y, phi0 ^ 4, u5) + *(24, phi, phih, xi_y, phi0 ^ 6, u4) + *(32, Fy, xi_y, phi0 ^ 2, u5, xi ^ 2) + *(72, Fy, xi_y, phi ^ 2, phi0 ^ 6, u ^ 7) + *(-96, Fy, phi, xi, xi_y, phi0 ^ 4, u6)) + *(*(-2, Fxpt) + *(u, *(Fxt, 4 + *(-2, B2p, u3) + *(8, B2, u4)) + *(u, Fxp ^ 2) + *(-4, Fxp, xi_x) + *(2, B2b, u) + *(4, B2t ^ 2, u5) + *(16, Fx ^ 2, u3) + *(-8, Fx, Fxp, u2) + *(-4, B2p, Fx ^ 2, u6) + *(-2, B2p, u, xi_xx) + *(-2, B2t, Fxp, u3) + *(4, Fx ^ 2, phi0 ^ 2, u5) + *(4, phi0 ^ 6, phit ^ 2, u3) + *(8, B2, xi_xx, u2) + *(12, Fx, u, xi_x) + *(16, B2t, xi_x, u2) + *(20, B2t, Fx, u4) + *(40, B2, u3, xi_x ^ 2) + *(64, B2 ^ 2, Fx ^ 2, u ^ 11) + *(64, B2 ^ 2, u ^ 7, xi_x ^ 2) + *(72, B2, Fx ^ 2, u ^ 7) + *(-16, xi, Fx ^ 2, phi0 ^ 2, u6) + *(-8, B2, Fx, Fxp, u6) + *(-8, B2, Fxp, xi_x, u4) + *(-4, B2p, Fx, xi_x, u4) + *(8, Fx, phit, phi0 ^ 4, u4) + *(16, Fx ^ 2, phi0 ^ 2, u ^ 7, xi ^ 2) + *(16, phi0 ^ 2, u3, xi ^ 2, xi_x ^ 2) + *(24, phi, Fx ^ 2, phi0 ^ 4, u ^ 7) + *(32, B2, B2t, Fx, u8) + *(32, B2, B2t, xi_x, u6) + *(36, Fx ^ 2, phi ^ 2, phi0 ^ 6, u ^ 9) + *(36, phi ^ 2, phi0 ^ 6, u5, xi_x ^ 2) + *(112, B2, Fx, xi_x, u5) + *(128, Fx, xi_x, B2 ^ 2, u ^ 9) + *(-48, phi, xi, Fx ^ 2, phi0 ^ 4, u8) + *(-48, phi, xi, phi0 ^ 4, u4, xi_x ^ 2) + *(-16, Fx, phit, xi, phi0 ^ 4, u5) + *(-16, Fx, xi, xi_x, phi0 ^ 2, u4) + *(-16, phit, xi, xi_x, phi0 ^ 4, u3) + *(24, Fx, phi, phit, phi0 ^ 6, u6) + *(24, Fx, phi, xi_x, phi0 ^ 4, u5) + *(24, phi, phit, xi_x, phi0 ^ 6, u4) + *(32, Fx, xi_x, phi0 ^ 2, u5, xi ^ 2) + *(72, Fx, xi_x, phi ^ 2, phi0 ^ 6, u ^ 7) + *(-96, Fx, phi, xi, xi_x, phi0 ^ 4, u6)), exp(*(2, B1, u4))), expB2u4, sinhGu4) + *(2, Fxph + Fypt + *(-1, u, *(Fxh, 2 + *(B1p, u3) + *(-1, B2p, u3) + *(-4, B1, u4) + *(4, B2, u4)) + *(Fyt, 2 + *(-1, B1p, u3) + *(-1, B2p, u3) + *(4, B1, u4) + *(4, B2, u4)) + *(-2, Fxp, xi_y) + *(-2, Fyp, xi_x) + *(2, B2c, u) + *(B1h, B2t, u5) + *(B1t, Fyp, u3) + *(Fxp, Fyp, u) + *(-1, B1h, Fxp, u3) + *(-1, B1t, B2h, u5) + *(-1, B2h, Fxp, u3) + *(-1, B2t, Fyp, u3) + *(-4, Fx, Fyp, u2) + *(-4, Fxp, Fy, u2) + *(-2, B1t, Fy, u4) + *(-2, B2p, u, xi_xy) + *(2, B1h, Fx, u4) + *(4, B2h, B2t, u5) + *(6, Fx, u, xi_y) + *(6, Fy, u, xi_x) + *(8, B2, xi_xy, u2) + *(8, B2h, xi_x, u2) + *(8, B2t, xi_y, u2) + *(10, B2h, Fx, u4) + *(10, B2t, Fy, u4) + *(16, Fx, Fy, u3) + *(-4, B1, B2h, Fx, u8) + *(-4, B1, B2h, xi_x, u6) + *(-4, B1, Fxp, Fy, u6) + *(-4, B1, Fxp, xi_y, u4) + *(-4, B1t, B2, Fy, u8) + *(-4, B1t, B2, xi_y, u6) + *(-4, B2, Fx, Fyp, u6) + *(-4, B2, Fxp, Fy, u6) + *(-4, B2, Fxp, xi_y, u4) + *(-4, B2, Fyp, xi_x, u4) + *(-4, B2p, Fx, Fy, u6) + *(-2, B1p, Fy, xi_x, u4) + *(-2, B2p, Fx, xi_y, u4) + *(-2, B2p, Fy, xi_x, u4) + *(2, B1p, Fx, xi_y, u4) + *(4, B1, B2t, Fy, u8) + *(4, B1, B2t, xi_y, u6) + *(4, B1, Fx, Fyp, u6) + *(4, B1, Fyp, xi_x, u4) + *(4, B1h, B2, Fx, u8) + *(4, B1h, B2, xi_x, u6) + *(4, Fx, Fy, phi0 ^ 2, u5) + *(4, Fx, phih, phi0 ^ 4, u4) + *(4, Fy, phit, phi0 ^ 4, u4) + *(4, phih, phit, phi0 ^ 6, u3) + *(16, B2, B2h, Fx, u8) + *(16, B2, B2h, xi_x, u6) + *(16, B2, B2t, Fy, u8) + *(16, B2, B2t, xi_y, u6) + *(40, B2, xi_x, xi_y, u3) + *(56, B2, Fx, xi_y, u5) + *(56, B2, Fy, xi_x, u5) + *(64, Fx, Fy, B2 ^ 2, u ^ 11) + *(64, Fx, xi_y, B2 ^ 2, u ^ 9) + *(64, Fy, xi_x, B2 ^ 2, u ^ 9) + *(64, xi_x, xi_y, B2 ^ 2, u ^ 7) + *(72, B2, Fx, Fy, u ^ 7) + *(-16, Fx, Fy, xi, phi0 ^ 2, u6) + *(-8, Fx, phih, xi, phi0 ^ 4, u5) + *(-8, Fx, xi, xi_y, phi0 ^ 2, u4) + *(-8, Fy, phit, xi, phi0 ^ 4, u5) + *(-8, Fy, xi, xi_x, phi0 ^ 2, u4) + *(-8, phih, xi, xi_x, phi0 ^ 4, u3) + *(-8, phit, xi, xi_y, phi0 ^ 4, u3) + *(12, Fx, phi, phih, phi0 ^ 6, u6) + *(12, Fx, phi, xi_y, phi0 ^ 4, u5) + *(12, Fy, phi, phit, phi0 ^ 6, u6) + *(12, Fy, phi, xi_x, phi0 ^ 4, u5) + *(12, phi, phih, xi_x, phi0 ^ 6, u4) + *(12, phi, phit, xi_y, phi0 ^ 6, u4) + *(16, Fx, Fy, phi0 ^ 2, u ^ 7, xi ^ 2) + *(16, Fx, xi_y, phi0 ^ 2, u5, xi ^ 2) + *(16, Fy, xi_x, phi0 ^ 2, u5, xi ^ 2) + *(16, xi_x, xi_y, phi0 ^ 2, u3, xi ^ 2) + *(24, Fx, Fy, phi, phi0 ^ 4, u ^ 7) + *(36, Fx, Fy, phi ^ 2, phi0 ^ 6, u ^ 9) + *(36, Fx, xi_y, phi ^ 2, phi0 ^ 6, u ^ 7) + *(36, Fy, xi_x, phi ^ 2, phi0 ^ 6, u ^ 7) + *(36, xi_x, xi_y, phi ^ 2, phi0 ^ 6, u5) + *(-48, Fx, Fy, phi, xi, phi0 ^ 4, u8) + *(-48, Fx, phi, xi, xi_y, phi0 ^ 4, u6) + *(-48, Fy, phi, xi, xi_x, phi0 ^ 4, u6) + *(-48, phi, xi, xi_x, xi_y, phi0 ^ 4, u4)), coshGu4, exp(*(u4, B1 + B2)))) + *(2/3, u2, *(*(Fxh, -3 + *(-3, Sp, u3) + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(Fyt, -3 + *(-3, Sp, u3) + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(u, *(6, Sc) + *(-6, Fx, xi_y) + *(-6, Fy, xi_x) + *(-6, Sp, xi_xy) + *(-6, Fx, Fyp, u) + *(-3, B1h, Fx, u3) + *(-3, B1t, Sh, u4) + *(-3, B2h, Fx, u3) + *(-3, B2t, Fy, u3) + *(2, Fxp, u, *(Fy, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(u, *(3, Sh) + *(2, xi, xi_y, phi0 ^ 2) + *(9, S, u, xi_y))) + *(3, B1h, St, u4) + *(3, B1t, Fy, u3) + *(3, B2h, St, u4) + *(3, B2t, Sh, u4) + *(4, xi, xi_xy, phi0 ^ 2) + *(4, xi_x, xi_y, phi0 ^ 2) + *(6, Fx, Sh, u3) + *(6, Fy, St, u3) + *(6, Fyp, St, u2) + *(12, Fx, Fy, u2) + *(18, S, u, xi_xy) + *(18, Sh, u, xi_x) + *(18, St, u, xi_y) + *(B1t, Fy, phi0 ^ 2, u5) + *(-1, B1h, Fx, phi0 ^ 2, u5) + *(-1, B2h, Fx, phi0 ^ 2, u5) + *(-1, B2t, Fy, phi0 ^ 2, u5) + *(-24, B2, Fx, Fy, u6) + *(-12, B1, Fx, Sh, u ^ 7) + *(-12, B1, Fx, xi_y, u4) + *(-12, B1, Sh, xi_x, u5) + *(-12, B2, Fx, xi_y, u4) + *(-12, B2, Fy, xi_x, u4) + *(-12, Fx, Fy, Sp, u5) + *(-9, B1t, Fy, S, u ^ 7) + *(-9, B1t, S, xi_y, u5) + *(-6, Fx, Sp, xi_y, u3) + *(-6, Fy, Sp, xi_x, u3) + *(-2, Fx, Fyp, phi0 ^ 2, u3) + *(-2, Fx, xi_y, phi0 ^ 2, u2) + *(-2, Fy, xi_x, phi0 ^ 2, u2) + *(9, B1h, Fx, S, u ^ 7) + *(9, B1h, S, xi_x, u5) + *(9, B2h, Fx, S, u ^ 7) + *(9, B2h, S, xi_x, u5) + *(9, B2t, Fy, S, u ^ 7) + *(9, B2t, S, xi_y, u5) + *(12, B1, Fy, St, u ^ 7) + *(12, B1, Fy, xi_x, u4) + *(12, B1, St, xi_y, u5) + *(12, B2, Fx, Sh, u ^ 7) + *(12, B2, Fy, St, u ^ 7) + *(12, B2, Sh, xi_x, u5) + *(12, B2, St, xi_y, u5) + *(18, Fx, Fyp, S, u5) + *(18, Fyp, S, xi_x, u3) + *(36, Fx, Fy, S, u6) + *(54, Fx, S, xi_y, u4) + *(54, Fy, S, xi_x, u4) + *(72, S, xi_x, xi_y, u2) + *(-8, B2, Fx, Fy, phi0 ^ 2, u8) + *(-4, B1, Fx, xi_y, phi0 ^ 2, u6) + *(-4, B2, Fx, xi_y, phi0 ^ 2, u6) + *(-4, B2, Fy, xi_x, phi0 ^ 2, u6) + *(-2, B1t, Fy, xi, phi0 ^ 2, u6) + *(-2, B1t, xi, xi_y, phi0 ^ 2, u4) + *(2, B1h, Fx, xi, phi0 ^ 2, u6) + *(2, B1h, xi, xi_x, phi0 ^ 2, u4) + *(2, B2h, Fx, xi, phi0 ^ 2, u6) + *(2, B2h, xi, xi_x, phi0 ^ 2, u4) + *(2, B2t, Fy, xi, phi0 ^ 2, u6) + *(2, B2t, xi, xi_y, phi0 ^ 2, u4) + *(4, B1, Fy, xi_x, phi0 ^ 2, u6) + *(4, Fx, Fy, xi, phi0 ^ 2, u5) + *(4, Fx, Fyp, xi, phi0 ^ 2, u4) + *(4, Fyp, xi, xi_x, phi0 ^ 2, u2) + *(8, Fx, xi, xi_y, phi0 ^ 2, u3) + *(8, Fy, xi, xi_x, phi0 ^ 2, u3) + *(12, u, xi, xi_x, xi_y, phi0 ^ 2) + *(72, B2, Fx, Fy, S, u ^ 10) + *(72, B2, Fx, S, xi_y, u8) + *(72, B2, Fy, S, xi_x, u8) + *(72, B2, S, xi_x, xi_y, u6) + *(16, B2, Fx, Fy, xi, phi0 ^ 2, u ^ 9) + *(16, B2, Fx, xi, xi_y, phi0 ^ 2, u ^ 7) + *(16, B2, Fy, xi, xi_x, phi0 ^ 2, u ^ 7) + *(16, B2, xi, xi_x, xi_y, phi0 ^ 2, u5)), coshGu4, exp(*(u4, B1 + B2))) + *(-1, *(Fyh, -3 + *(-3, Sp, u3) + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(u, *(3, Ss) + *(Fy, *(-1, xi_y, 6 + *(-54, S, u4) + *(6, Sp, u3) + *(12, B2, u4) + *(-72, B2, S, u8)) + *(6, Sh, u3 + *(2, B2, u ^ 7)) + *(B2h, u3, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(2, Fyp, u, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(2, xi_y, phi0 ^ 2, u2, 1 + *(2, B2, u4), -1 + *(4, u, xi))) + *(-3, Sp, xi_yy) + *(2, phi0 ^ 2, xi_y ^ 2) + *(2, xi, xi_yy, phi0 ^ 2) + *(2, Fy ^ 2, u2, 3 + *(-3, Sp, u3) + *(9, S, u4) + *(xi, phi0 ^ 2, u3) + *(2, B2, u4, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi)))) + *(3, B2h, Sh, u4) + *(6, Fyp, Sh, u2) + *(9, S, u, xi_yy) + *(18, Sh, u, xi_y) + *(36, S, u2, xi_y ^ 2) + *(6, u, xi, phi0 ^ 2, xi_y ^ 2) + *(9, B2h, S, xi_y, u5) + *(12, B2, Sh, xi_y, u5) + *(18, Fyp, S, xi_y, u3) + *(36, B2, S, u6, xi_y ^ 2) + *(2, B2h, xi, xi_y, phi0 ^ 2, u4) + *(4, Fyp, xi, xi_y, phi0 ^ 2, u2) + *(8, B2, xi, phi0 ^ 2, u5, xi_y ^ 2)) + *(Fxt, -3 + *(-3, Sp, u3) + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi)), exp(*(2, B1, u4))) + *(u, *(3, Sb) + *(Fx, *(-1, xi_x, 6 + *(-54, S, u4) + *(6, Sp, u3) + *(12, B2, u4) + *(-72, B2, S, u8)) + *(6, St, u3 + *(2, B2, u ^ 7)) + *(B2t, u3, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(2, Fxp, u, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(2, xi_x, phi0 ^ 2, u2, 1 + *(2, B2, u4), -1 + *(4, u, xi))) + *(-3, Sp, xi_xx) + *(2, phi0 ^ 2, xi_x ^ 2) + *(2, xi, xi_xx, phi0 ^ 2) + *(2, Fx ^ 2, u2, 3 + *(-3, Sp, u3) + *(9, S, u4) + *(xi, phi0 ^ 2, u3) + *(2, B2, u4, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi)))) + *(3, B2t, St, u4) + *(6, Fxp, St, u2) + *(9, S, u, xi_xx) + *(18, St, u, xi_x) + *(36, S, u2, xi_x ^ 2) + *(6, u, xi, phi0 ^ 2, xi_x ^ 2) + *(9, B2t, S, xi_x, u5) + *(12, B2, St, xi_x, u5) + *(18, Fxp, S, xi_x, u3) + *(36, B2, S, u6, xi_x ^ 2) + *(2, B2t, xi, xi_x, phi0 ^ 2, u4) + *(4, Fxp, xi, xi_x, phi0 ^ 2, u2) + *(8, B2, xi, phi0 ^ 2, u5, xi_x ^ 2), exp(*(2, B1, u4))), expB2u4, sinhGu4), 3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi))) + *(4/3, u5, *((*(Fy, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(u, *(3, Sh) + *(2, xi, xi_y, phi0 ^ 2) + *(9, S, u, xi_y))) ^ 2, sinhGu4) + *((*(Fx, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(u, *(3, St) + *(2, xi, xi_x, phi0 ^ 2) + *(9, S, u, xi_x))) ^ 2, exp(*(2, B1, u4)), sinhGu4) + *(-2, *(Fx, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(u, *(3, St) + *(2, xi, xi_x, phi0 ^ 2) + *(9, S, u, xi_x)), *(Fy, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(u, *(3, Sh) + *(2, xi, xi_y, phi0 ^ 2) + *(9, S, u, xi_y)), coshGu4, expB1u4), expB2u4) + *(1/9, (3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi))) ^ 3, Gp + *(-4, G, u), *(3, (1 + *(u, xi)) ^ 2) + *(-1, phi0 ^ 2, u2) + *(6, Sd, u4), expB1u4)

    nothing
end

function phid_eq_coeff!(ABCS::Vector, vars::phidVars, ::Inner)
    @unpack (
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


function A_eq_coeff!(ABCS::Vector, vars::AVars, ::Inner)
    @unpack (
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

    u2 = u*u
    u3 = u*u2
    u4 = u2*u2
    u5 = u4*u
    u5 = u3*u2
    u6 = u3*u3
    u7 = u4*u3
    u8 = u4*u4

    coshGu4 = cosh(*(G, u4))
    sinhGu4 = sinh(*(G, u4))

    coshGu4sq = cosh(*(G, u4)) ^ 2
    sinhGu4sq = sinh(*(G, u4)) ^ 2

    expB1u4 = exp(*(B1, u4))
    expB2u4 = exp(*(B2, u4))

    ABCS[1] = *(2/27, (*(3, u, 1 + *(S, u4) + *(u, xi)) + *(phi0 ^ 2, u3, -1 + *(u, xi))) ^ 4, expB1u4)

    ABCS[2] = *(4/9, u3, (3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi))) ^ 4, expB1u4)

    ABCS[3] = *(4/9, u2, (3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi))) ^ 4, expB1u4)

    ABCS[4] = *(-4/9, *(486, xi ^ 2) + *(-405, S ^ 2, u6) + *(-243, S ^ 3, u ^ 10) + *(-81, S, u2) + *(81, u2, xi ^ 4) + *(324, u, xi ^ 3) + *(-1377, S ^ 2, u8, xi ^ 2) + *(-1296, xi, S ^ 2, u7) + *(-972, S, u4, xi ^ 2) + *(-810, S, u5, xi ^ 3) + *(-486, S, xi, u3) + *(-486, xi, S ^ 3, u ^ 11) + *(-486, S ^ 2, u ^ 9, xi ^ 3) + *(-243, S, u6, xi ^ 4) + *(-243, S ^ 3, u ^ 12, xi ^ 2) + *(9, phi0 ^ 4, u4, *(7, xi ^ 2) + *(-11, S, u2) + *(-7, S ^ 2, u6) + *(-5, u2, xi ^ 4) + *(-4, u3, xi ^ 5) + *(6, u, xi ^ 3) + *(-7, S, u6, xi ^ 4) + *(-2, S, u5, xi ^ 3) + *(2, S, xi, u3) + *(8, xi, S ^ 2, u7) + *(22, S, u4, xi ^ 2) + *(Sp, u, -1 + *(u, xi), -3 + *(u2, xi ^ 2) + *(u3, xi ^ 3) + *(-3, u, xi) + *(-2, S, u4) + *(2, Sd, u4, -1 + *(u, xi))) + *(-2, Sd, u2, -1 + *(u, xi), -1 + *(u, xi) + *(4, u2, xi ^ 2) + *(S, u4, -5 + *(7, u, xi)))) + *(27, phi0 ^ 2, 1 + *(S, u4) + *(u, xi), -2 + *(-1, u2, xi ^ 2) + *(-3, u3, xi ^ 3) + *(-2, u, xi) + *(-2, u4, xi ^ 4) + *(3, S ^ 2, u8) + *(9, S, u4) + *(Sp, u3, -3 + *(-1, S, u4) + *(-3, u, xi) + *(2, u2, xi ^ 2) + *(2, u3, xi ^ 3) + *(4, Sd, u4, -1 + *(u, xi))) + *(-9, S, u6, xi ^ 2) + *(-8, S, u7, xi ^ 3) + *(-2, Sd, u4, 1 + *(-1, u, xi) + *(2, u2, xi ^ 2) + *(S, u4, -7 + *(8, u, xi))) + *(9, S, xi, u5)) + *(phi0 ^ 8, u6, (-1 + *(u, xi)) ^ 2, -1 + *(2, u, xi)) + *(-162, Sd, u2, (1 + *(S, u4) + *(u, xi)) ^ 2, -1 + *(3, S, u4)) + *(-3, phi0 ^ 6, u4, -1 + *(u, xi), 2 + *(u3, xi ^ 3) + *(-7, u2, xi ^ 2) + *(-2, u, xi) + *(2, u4, xi ^ 4) + *(5, S, u4) + *(Sp, u3, -1 + *(u, xi)) + *(-7, S, xi, u5) + *(2, Sd, u4, 1 + *(-3, u, xi) + *(2, u2, xi ^ 2))) + *(81, Sp, u, (1 + *(S, u4) + *(u, xi)) ^ 2, (1 + *(u, xi)) ^ 2 + *(2, Sd, u4)), expB1u4) + *(2/81, *(5832, xi ^ 2) + *(-324, UU(phi, potential), u2) + *(972, u2, xi ^ 4) + *(3888, u, xi ^ 3) + *(-2916, B2, B2d, u6) + *(-1944, UU(phi, potential), u4, xi ^ 2) + *(-1296, UU(phi, potential), xi, u3) + *(-1296, UU(phi, potential), u5, xi ^ 3) + *(-972, G, Gd, u6) + *(-486, u6, (S + *(S, u, xi)) ^ 2, -12 + *(4, UU(phi, potential), u4) + *(-9, B2d, u7, B2p + *(-4, B2, u)) + *(-3, Gd, Gp, u7) + *(12, G, Gd, u8) + *(-3, B1d, u7, coshGu4sq, B1p + *(-4, B1, u))) + *(-324, UU(phi, potential), u6, xi ^ 4) + *(-81, S ^ 4, u ^ 14, -12 + *(4, UU(phi, potential), u4) + *(-9, B2d, u7, B2p + *(-4, B2, u)) + *(-3, Gd, Gp, u7) + *(12, G, Gd, u8) + *(-3, B1d, u7, coshGu4sq, B1p + *(-4, B1, u))) + *(-54, u, phi0 ^ 4, *(9, phip, (1 + *(S, u4) + *(u, xi)) ^ 4) + *(u, (1 + *(S, u4) + *(u, xi)) ^ 2, *(-1 + *(u, xi), -12 + *(-4, UU(phi, potential), u4) + *(12, u, xi) + *(48, u2, xi ^ 2) + *(-12, G, Gd, u8) + *(3, Gd, Gp, u7) + *(4, UU(phi, potential), xi, u5) + *(24, S, u4, -1 + *(2, u, xi)) + *(-12, B1, B1d, u8, coshGu4sq) + *(-9, B2d, u7, -1 + *(u, xi), B2p + *(-4, B2, u)) + *(-3, Gd, Gp, xi, u8) + *(3, B1d, B1p, u7, coshGu4sq) + *(12, G, Gd, xi, u ^ 9) + *(-3, B1d, B1p, xi, u8, coshGu4sq) + *(12, B1, B1d, xi, u ^ 9, coshGu4sq)) + *(-27, phi, (1 + *(S, u4) + *(u, xi)) ^ 2) + *(-18, phid, (1 + *(S, u4) + *(u, xi)) ^ 2, -1 + *(2, u, xi)))) + *(108, phi0 ^ 2, (1 + *(S, u4) + *(u, xi)) ^ 3, -3 + *(-18, u2, xi ^ 2) + *(3, u, xi) + *(4, UU(phi, potential), u4) + *(-4, UU(phi, potential), xi, u5) + *(-3, Gd, Gp, u7) + *(9, S, u4, 1 + *(-2, u, xi)) + *(12, G, Gd, u8) + *(-12, G, Gd, xi, u ^ 9) + *(-3, B1d, B1p, u7, coshGu4sq) + *(3, Gd, Gp, xi, u8) + *(9, B2d, u7, -1 + *(u, xi), B2p + *(-4, B2, u)) + *(12, B1, B1d, u8, coshGu4sq) + *(-12, B1, B1d, xi, u ^ 9, coshGu4sq) + *(3, B1d, B1p, xi, u8, coshGu4sq)) + *(243, Gd, Gp, u5) + *(729, B2d, B2p, u5) + *(phi0 ^ 8, u5, -1 + *(u, xi), *(-1, u, *(-1 + *(u, xi), *(-1 + *(u, xi), -132 + *(-4, UU(phi, potential), u4) + *(132, u, xi) + *(288, u2, xi ^ 2) + *(-12, G, Gd, u8) + *(3, Gd, Gp, u7) + *(4, UU(phi, potential), xi, u5) + *(144, S, u4, -1 + *(2, u, xi)) + *(-12, B1, B1d, u8, coshGu4sq) + *(-9, B2d, u7, -1 + *(u, xi), B2p + *(-4, B2, u)) + *(-3, Gd, Gp, xi, u8) + *(3, B1d, B1p, u7, coshGu4sq) + *(12, G, Gd, xi, u ^ 9) + *(-3, B1d, B1p, xi, u8, coshGu4sq) + *(12, B1, B1d, xi, u ^ 9, coshGu4sq)) + *(-648, phid, (1 + *(S, u4) + *(u, xi)) ^ 2, -1 + *(2, u, xi))) + *(972, phi, (1 + *(S, u4) + *(u, xi)) ^ 2, 1 + *(-1, u, xi) + *(4, phid, 1 + *(S, u4) + *(u, xi)))) + *(324, phip, (1 + *(S, u4) + *(u, xi)) ^ 2, 1 + *(-1, u, xi) + *(4, phid, 1 + *(S, u4) + *(u, xi)))) + *(-17496, B2, B2d, u8, xi ^ 2) + *(-11664, B2, B2d, xi, u7) + *(-11664, B2, B2d, u ^ 9, xi ^ 3) + *(-5832, G, Gd, u8, xi ^ 2) + *(-3888, G, Gd, xi, u7) + *(-3888, G, Gd, u ^ 9, xi ^ 3) + *(-2916, B2, B2d, u ^ 10, xi ^ 4) + *(-972, B1, B1d, u6, coshGu4sq) + *(-972, G, Gd, u ^ 10, xi ^ 4) + *(-324, S, u2, (1 + *(u, xi)) ^ 3, -12 + *(4, UU(phi, potential), u4) + *(-9, B2d, u7, B2p + *(-4, B2, u)) + *(-3, Gd, Gp, u7) + *(12, G, Gd, u8) + *(-3, B1d, u7, coshGu4sq, B1p + *(-4, B1, u))) + *(-324, S ^ 3, u ^ 10, 1 + *(u, xi), -12 + *(4, UU(phi, potential), u4) + *(-9, B2d, u7, B2p + *(-4, B2, u)) + *(-3, Gd, Gp, u7) + *(12, G, Gd, u8) + *(-3, B1d, u7, coshGu4sq, B1p + *(-4, B1, u))) + *(-6, phi0 ^ 12, u ^ 9, (-1 + *(u, xi)) ^ 3, *(-1, phip) + *(phi, u, 3 + *(-3, u, xi) + *(72, phid, 1 + *(S, u4) + *(u, xi))) + *(phip, u, xi) + *(-24, phid, phip, 1 + *(S, u4) + *(u, xi)) + *(-2, phid, u, 1 + *(-3, u, xi) + *(2, u2, xi ^ 2))) + *(12, phi0 ^ 6, u3, *(-1, u, *(-1 + *(u, xi), *(-1 + *(u, xi), -42 + *(-4, UU(phi, potential), u4) + *(42, u, xi) + *(108, u2, xi ^ 2) + *(-12, G, Gd, u8) + *(3, Gd, Gp, u7) + *(4, UU(phi, potential), xi, u5) + *(54, S, u4, -1 + *(2, u, xi)) + *(-12, B1, B1d, u8, coshGu4sq) + *(-9, B2d, u7, -1 + *(u, xi), B2p + *(-4, B2, u)) + *(-3, Gd, Gp, xi, u8) + *(3, B1d, B1p, u7, coshGu4sq) + *(12, G, Gd, xi, u ^ 9) + *(-3, B1d, B1p, xi, u8, coshGu4sq) + *(12, B1, B1d, xi, u ^ 9, coshGu4sq)) + *(-108, phid, (1 + *(S, u4) + *(u, xi)) ^ 2, -1 + *(2, u, xi))) + *(81, phi, (1 + *(S, u4) + *(u, xi)) ^ 2, 2 + *(-2, u, xi) + *(3, phid, 1 + *(S, u4) + *(u, xi)))) + *(27, phip, (1 + *(S, u4) + *(u, xi)) ^ 2, 2 + *(-2, u, xi) + *(3, phid, 1 + *(S, u4) + *(u, xi))), 1 + *(S, u4) + *(u, xi)) + *(12, phi0 ^ 10, u7, (-1 + *(u, xi)) ^ 2, *(u, 1 + *(-1, u, xi) + *(12, phid, 1 + *(S, u4) + *(u, xi)), 1 + *(-3, u, xi) + *(2, u2, xi ^ 2)) + *(6, phip, 1 + *(S, u4) + *(u, xi), 1 + *(-1, u, xi) + *(9, phid, 1 + *(S, u4) + *(u, xi))) + *(-18, phi, u, 1 + *(S, u4) + *(u, xi), 1 + *(-1, u, xi) + *(9, phid, 1 + *(S, u4) + *(u, xi)))) + *(243, B1d, B1p, u5, coshGu4sq) + *(243, Gd, Gp, u ^ 9, xi ^ 4) + *(729, B2d, B2p, u ^ 9, xi ^ 4) + *(972, Gd, Gp, xi, u6) + *(972, Gd, Gp, u8, xi ^ 3) + *(1458, Gd, Gp, u7, xi ^ 2) + *(2916, B2d, B2p, xi, u6) + *(2916, B2d, B2p, u8, xi ^ 3) + *(4374, B2d, B2p, u7, xi ^ 2) + *(-5832, B1, B1d, u8, xi ^ 2, coshGu4sq) + *(-3888, B1, B1d, xi, u7, coshGu4sq) + *(-3888, B1, B1d, u ^ 9, xi ^ 3, coshGu4sq) + *(-972, B1, B1d, u ^ 10, xi ^ 4, coshGu4sq) + *(12, phid, phi0 ^ 14, u ^ 11, (-1 + *(u, xi)) ^ 4, phip + *(-3, phi, u)) + *(243, B1d, B1p, u ^ 9, xi ^ 4, coshGu4sq) + *(972, B1d, B1p, xi, u6, coshGu4sq) + *(972, B1d, B1p, u8, xi ^ 3, coshGu4sq) + *(1458, B1d, B1p, u7, xi ^ 2, coshGu4sq), expB1u4) + *(4/3, u6, *((*(Fx, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(u, *(3, St) + *(2, xi, xi_x, phi0 ^ 2) + *(9, S, u, xi_x))) ^ 2, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *((*(Fy, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(u, *(3, Sh) + *(2, xi, xi_y, phi0 ^ 2) + *(9, S, u, xi_y))) ^ 2, coshGu4, expB2u4) + *(-2, *(Fx, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(u, *(3, St) + *(2, xi, xi_x, phi0 ^ 2) + *(9, S, u, xi_x)), *(Fy, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(u, *(3, Sh) + *(2, xi, xi_y, phi0 ^ 2) + *(9, S, u, xi_y)), exp(*(B1, u4) + *(B2, u4)), sinhGu4)) + *(4/27, *(81, S ^ 4, u ^ 14) + *(81, xi ^ 2, 6 + *(u2, xi ^ 2) + *(4, u, xi)) + *(486, u6, (S + *(S, u, xi)) ^ 2) + *(phi0 ^ 8, u6, (-1 + *(u, xi)) ^ 4) + *(108, phi0 ^ 2, (1 + *(S, u4) + *(u, xi)) ^ 3, -1 + *(u, xi)) + *(324, S, u2, (1 + *(u, xi)) ^ 3) + *(324, S ^ 3, u ^ 10, 1 + *(u, xi)) + *(12, phi0 ^ 6, u4, (-1 + *(u, xi)) ^ 3, 1 + *(S, u4) + *(u, xi)) + *(54, phi0 ^ 4, u2, (-1 + *(u, xi)) ^ 2, (1 + *(S, u4) + *(u, xi)) ^ 2), expB1u4) + *(3, u4, (1 + *(S, u4) + *(u, xi) + *(1/3, phi0 ^ 2, u2, -1 + *(u, xi))) ^ 2, *((Fxp + *(-2, Fx, u)) ^ 2, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *((Fyp + *(-2, Fy, u)) ^ 2, coshGu4, expB2u4) + *(-2, B1b + *(4, u, xi_x + *(Fx, u2), *(2, B1t) + *(5, B1, u, xi_x + *(Fx, u2))), coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-2, B2s + *(4, u, xi_y + *(Fy, u2), *(2, B2h) + *(5, B2, u, xi_y + *(Fy, u2))), coshGu4, expB2u4) + *(-2, B2b + *(4, u, xi_x + *(Fx, u2), *(2, B2t) + *(5, B2, u, xi_x + *(Fx, u2))), coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-2, Gs + *(4, u, xi_y + *(Fy, u2), *(2, Gh) + *(5, G, u, xi_y + *(Fy, u2))), expB2u4, sinhGu4) + *(-2, Gb + *(4, u, xi_x + *(Fx, u2), *(2, Gt) + *(5, G, u, xi_x + *(Fx, u2))), exp(*(B2, u4) + *(2, B1, u4)), sinhGu4) + *(2, B1s + *(4, u, xi_y + *(Fy, u2), *(2, B1h) + *(5, B1, u, xi_y + *(Fy, u2))), coshGu4, expB2u4) + *(4, B2c + *(4, u, *(B2h, xi_x + *(Fx, u2)) + *(B2t + *(5, B2, u, xi_x + *(Fx, u2)), xi_y + *(Fy, u2))), exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(4, Gc + *(4, u, *(Gh, xi_x) + *(Gt, xi_y) + *(Fx, u2, Gh + *(5, G, u, xi_y + *(Fy, u2))) + *(Fy, u2, Gt + *(5, G, u, xi_x)) + *(5, G, u, xi_x, xi_y)), coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-4, u4, (B2h + *(4, B2, u, xi_y + *(Fy, u2))) ^ 2, coshGu4, expB2u4) + *(-4, u4, (B2t + *(4, B2, u, xi_x + *(Fx, u2))) ^ 2, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-2, u4, (B1h + *(4, B1, u, xi_y + *(Fy, u2))) ^ 2, coshGu4, expB2u4) + *(-2, u4, (B1t + *(4, B1, u, xi_x + *(Fx, u2))) ^ 2, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-2, u4, (Gh + *(4, G, u, xi_y + *(Fy, u2))) ^ 2, coshGu4, expB2u4) + *(-2, u4, (Gt + *(4, G, u, xi_x + *(Fx, u2))) ^ 2, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-2, B1p + *(-4, B1, u), xi_yy + *(Fyh, u2) + *(2, Fy ^ 2, u5) + *(2, Fy, xi_y, u3), coshGu4, expB2u4) + *(-2, B2p + *(-4, B2, u), xi_xy + *(Fxh, u2) + *(2, Fx, u3, xi_y + *(Fy, u2)), exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-2, B2p + *(-4, B2, u), xi_xy + *(Fyt, u2) + *(2, Fx, Fy, u5) + *(2, Fy, xi_x, u3), exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-2, Gp + *(-4, G, u), xi_xy + *(Fxh, u2) + *(2, Fx, u3, xi_y + *(Fy, u2)), coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-2, Gp + *(-4, G, u), xi_xy + *(Fyt, u2) + *(2, Fx, Fy, u5) + *(2, Fy, xi_x, u3), coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(2, B1p + *(-4, B1, u), xi_xx + *(Fxt, u2) + *(2, Fx ^ 2, u5) + *(2, Fx, xi_x, u3), coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(2, B2p + *(-4, B2, u), xi_xx + *(Fxt, u2) + *(2, Fx ^ 2, u5) + *(2, Fx, xi_x, u3), coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(2, B2p + *(-4, B2, u), xi_yy + *(Fyh, u2) + *(2, Fy ^ 2, u5) + *(2, Fy, xi_y, u3), coshGu4, expB2u4) + *(2, Fxp + *(-2, Fx, u), *(-1, Fyp) + *(2, Fy, u), exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(2, Gp + *(-4, G, u), xi_xx + *(Fxt, u2) + *(2, Fx ^ 2, u5) + *(2, Fx, xi_x, u3), exp(*(B2, u4) + *(2, B1, u4)), sinhGu4) + *(2, Gp + *(-4, G, u), xi_yy + *(Fyh, u2) + *(2, Fy ^ 2, u5) + *(2, Fy, xi_y, u3), expB2u4, sinhGu4) + *(-4, phi0 ^ 2, u2, (*(Fx, u + *(-2, xi, u2)) + *(phi0 ^ 2, phit + *(3, phi, u, xi_x + *(Fx, u2))) + *(-2, xi, xi_x)) ^ 2, coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-4, phi0 ^ 2, u2, (*(Fy, u + *(-2, xi, u2)) + *(phi0 ^ 2, phih + *(3, phi, u, xi_y + *(Fy, u2))) + *(-2, xi, xi_y)) ^ 2, coshGu4, expB2u4) + *(-4, u4, B1t + *(4, B1, u, xi_x + *(Fx, u2)), Gt + *(4, G, u, xi_x + *(Fx, u2)), exp(*(B2, u4) + *(2, B1, u4)), sinhGu4) + *(-2, u4, B1h + *(4, B1, u, xi_y + *(Fy, u2)), Gt + *(4, G, u, xi_x + *(Fx, u2)), coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(-2, u4, B1t + *(4, B1, u, xi_x + *(Fx, u2)), B2t + *(4, B2, u, xi_x + *(Fx, u2)), coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(-2, u4, B2h + *(4, B2, u, xi_y + *(Fy, u2)), Gh + *(4, G, u, xi_y + *(Fy, u2)), expB2u4, sinhGu4) + *(-2, u4, B2t + *(4, B2, u, xi_x + *(Fx, u2)), Gt + *(4, G, u, xi_x + *(Fx, u2)), exp(*(B2, u4) + *(2, B1, u4)), sinhGu4) + *(2, u4, B1h + *(4, B1, u, xi_y + *(Fy, u2)), B2h + *(4, B2, u, xi_y + *(Fy, u2)), coshGu4, expB2u4) + *(2, u4, B1t + *(4, B1, u, xi_x + *(Fx, u2)), Gh + *(4, G, u, xi_y + *(Fy, u2)), coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(2, u4, B2h + *(4, B2, u, xi_y + *(Fy, u2)), Gt + *(4, G, u, xi_x + *(Fx, u2)), coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(2, u4, B2t + *(4, B2, u, xi_x + *(Fx, u2)), Gh + *(4, G, u, xi_y + *(Fy, u2)), coshGu4, exp(*(B1, u4) + *(B2, u4))) + *(4, u4, B1h + *(4, B1, u, xi_y + *(Fy, u2)), Gh + *(4, G, u, xi_y + *(Fy, u2)), expB2u4, sinhGu4) + *(4, u4, Gh + *(4, G, u, xi_y + *(Fy, u2)), Gt + *(4, G, u, xi_x + *(Fx, u2)), exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(8, u4, B2h + *(4, B2, u, xi_y + *(Fy, u2)), B2t + *(4, B2, u, xi_x + *(Fx, u2)), exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(8, phi0 ^ 2, u2, *(Fx, u + *(-2, xi, u2)) + *(phi0 ^ 2, phit + *(3, phi, u, xi_x + *(Fx, u2))) + *(-2, xi, xi_x), *(Fy, u + *(-2, xi, u2)) + *(phi0 ^ 2, phih + *(3, phi, u, xi_y + *(Fy, u2))) + *(-2, xi, xi_y), exp(*(B1, u4) + *(B2, u4)), sinhGu4)) + *(-8/3, u3, 3 + *(3, S, u4) + *(3, u, xi) + *(phi0 ^ 2, u2, -1 + *(u, xi)), *(u4, *(Gh + *(4, G, u, xi_y + *(Fy, u2)), *(Fy, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(u, *(3, Sh) + *(2, xi, xi_y, phi0 ^ 2) + *(9, S, u, xi_y)), expB2u4) + *(Gt + *(4, G, u, xi_x + *(Fx, u2)), *(Fx, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(u, *(3, St) + *(2, xi, xi_x, phi0 ^ 2) + *(9, S, u, xi_x)), exp(*(B2, u4) + *(2, B1, u4))), sinhGu4) + *(*(Fxt, -3 + *(-3, Sp, u3) + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(u, *(3, Sb) + *(Fx, *(2, St, *(9, u3) + *(u7, *(6, B1) + *(6, B2))) + *(2, xi_x, -3 + *(-6, B1, u4) + *(-6, B2, u4) + *(-3, Sp, u3) + *(45, S, u4) + *(phi0 ^ 2, u2, -1 + *(8, u, xi) + *(2, B1, u4, -1 + *(4, u, xi)) + *(2, B2, u4, -1 + *(4, u, xi))) + *(36, B1, S, u8) + *(36, B2, S, u8)) + *(B1t, u3, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(B2t, u3, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi)))) + *(-3, Sp, xi_xx) + *(2, phi0 ^ 2, xi_x ^ 2) + *(2, xi, xi_xx, phi0 ^ 2) + *(2, Fx ^ 2, u2, -3 + *(-6, B1, u4) + *(-6, B2, u4) + *(-3, Sp, u3) + *(27, S, u4) + *(phi0 ^ 2, u2, -2 + *(5, u, xi) + *(2, B1, u4, -1 + *(2, u, xi)) + *(2, B2, u4, -1 + *(2, u, xi))) + *(18, B1, S, u8) + *(18, B2, S, u8)) + *(3, B1t, St, u4) + *(3, B2t, St, u4) + *(9, S, u, xi_xx) + *(18, St, u, xi_x) + *(36, S, u2, xi_x ^ 2) + *(6, u, xi, phi0 ^ 2, xi_x ^ 2) + *(9, B1t, S, xi_x, u5) + *(9, B2t, S, xi_x, u5) + *(12, B1, St, xi_x, u5) + *(12, B2, St, xi_x, u5) + *(36, B1, S, u6, xi_x ^ 2) + *(36, B2, S, u6, xi_x ^ 2) + *(2, B1t, xi, xi_x, phi0 ^ 2, u4) + *(2, B2t, xi, xi_x, phi0 ^ 2, u4) + *(8, B1, xi, phi0 ^ 2, u5, xi_x ^ 2) + *(8, B2, xi, phi0 ^ 2, u5, xi_x ^ 2)), coshGu4, exp(*(B2, u4) + *(2, B1, u4))) + *(*(Fyh, -3 + *(-3, Sp, u3) + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(u, *(3, Ss) + *(Fy, *(B1h, *(3, u3) + *(phi0 ^ 2, u5) + *(-9, S, u7) + *(-2, xi, phi0 ^ 2, u6)) + *(-2, Sh, *(-9, u3) + *(u7, *(-6, B2) + *(6, B1))) + *(-2, xi_y, 3 + *(-45, S, u4) + *(-6, B1, u4) + *(3, Sp, u3) + *(6, B2, u4) + *(phi0 ^ 2, u2, 1 + *(-8, u, xi) + *(2, B1, u4, -1 + *(4, u, xi)) + *(2, B2, u4, 1 + *(-4, u, xi))) + *(-36, B2, S, u8) + *(36, B1, S, u8)) + *(B2h, u3, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi)))) + *(-3, Sp, xi_yy) + *(2, phi0 ^ 2, xi_y ^ 2) + *(-3, B1h, Sh, u4) + *(-2, Fy ^ 2, u2, 3 + *(-27, S, u4) + *(-6, B1, u4) + *(3, Sp, u3) + *(6, B2, u4) + *(phi0 ^ 2, u2, 2 + *(-5, u, xi) + *(2, B1, u4, -1 + *(2, u, xi)) + *(2, B2, u4, 1 + *(-2, u, xi))) + *(-18, B2, S, u8) + *(18, B1, S, u8)) + *(2, xi, xi_yy, phi0 ^ 2) + *(3, B2h, Sh, u4) + *(9, S, u, xi_yy) + *(18, Sh, u, xi_y) + *(36, S, u2, xi_y ^ 2) + *(-36, B1, S, u6, xi_y ^ 2) + *(-12, B1, Sh, xi_y, u5) + *(-9, B1h, S, xi_y, u5) + *(6, u, xi, phi0 ^ 2, xi_y ^ 2) + *(9, B2h, S, xi_y, u5) + *(12, B2, Sh, xi_y, u5) + *(36, B2, S, u6, xi_y ^ 2) + *(-8, B1, xi, phi0 ^ 2, u5, xi_y ^ 2) + *(-2, B1h, xi, xi_y, phi0 ^ 2, u4) + *(2, B2h, xi, xi_y, phi0 ^ 2, u4) + *(8, B2, xi, phi0 ^ 2, u5, xi_y ^ 2)), coshGu4, expB2u4) + *(-1, *(Fxh, -3 + *(-3, Sp, u3) + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(Fyt, -3 + *(-3, Sp, u3) + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(u, *(6, Sc) + *(Fx, *(2, Sh, *(9, u3) + *(6, B2, u7)) + *(2, xi_y, -3 + *(-6, B2, u4) + *(-3, Sp, u3) + *(45, S, u4) + *(phi0 ^ 2, u2, -1 + *(8, u, xi) + *(2, B2, u4, -1 + *(4, u, xi))) + *(36, B2, S, u8)) + *(4, Fy, *(-3, u2, 1 + *(Sp, u3) + *(-9, S, u4) + *(2, B2, u4) + *(-6, B2, S, u8)) + *(phi0 ^ 2, u4, -2 + *(5, u, xi) + *(2, B2, u4, -1 + *(2, u, xi)))) + *(B2h, u3, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi)))) + *(-6, Fy, xi_x) + *(-6, Sp, xi_xy) + *(-3, B2t, Fy, u3) + *(3, B2h, St, u4) + *(3, B2t, Sh, u4) + *(4, xi, xi_xy, phi0 ^ 2) + *(4, xi_x, xi_y, phi0 ^ 2) + *(18, Fy, St, u3) + *(18, S, u, xi_xy) + *(18, Sh, u, xi_x) + *(18, St, u, xi_y) + *(-1, B2t, Fy, phi0 ^ 2, u5) + *(-12, B2, Fy, xi_x, u4) + *(-6, Fy, Sp, xi_x, u3) + *(-2, Fy, xi_x, phi0 ^ 2, u2) + *(9, B2h, S, xi_x, u5) + *(9, B2t, Fy, S, u7) + *(9, B2t, S, xi_y, u5) + *(12, B2, Fy, St, u7) + *(12, B2, Sh, xi_x, u5) + *(12, B2, St, xi_y, u5) + *(72, S, xi_x, xi_y, u2) + *(90, Fy, S, xi_x, u4) + *(-4, B2, Fy, xi_x, phi0 ^ 2, u6) + *(2, B2h, xi, xi_x, phi0 ^ 2, u4) + *(2, B2t, Fy, xi, phi0 ^ 2, u6) + *(2, B2t, xi, xi_y, phi0 ^ 2, u4) + *(12, u, xi, xi_x, xi_y, phi0 ^ 2) + *(16, Fy, xi, xi_x, phi0 ^ 2, u3) + *(72, B2, Fy, S, xi_x, u8) + *(72, B2, S, xi_x, xi_y, u6) + *(16, B2, Fy, xi, xi_x, phi0 ^ 2, u7) + *(16, B2, xi, xi_x, xi_y, phi0 ^ 2, u5)), exp(*(B1, u4) + *(B2, u4)), sinhGu4) + *(-1, u4, *(Fx, Gh, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(Fy, Gt, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi))) + *(Gh, u, *(3, St) + *(2, xi, xi_x, phi0 ^ 2) + *(9, S, u, xi_x)) + *(Gt, u, *(3, Sh) + *(2, xi, xi_y, phi0 ^ 2) + *(9, S, u, xi_y)) + *(4, G, u2, *(xi_y, *(3, St) + *(4, xi, xi_x, phi0 ^ 2) + *(18, S, u, xi_x)) + *(3, Sh, xi_x)) + *(4, Fx, G, u, *(xi_y, -3 + *(18, S, u4) + *(phi0 ^ 2, u2, -1 + *(4, u, xi))) + *(3, Sh, u3) + *(2, Fy, u2, -3 + *(9, S, u4) + *(phi0 ^ 2, u2, -1 + *(2, u, xi)))) + *(4, Fy, G, u, *(xi_x, -3 + *(18, S, u4) + *(phi0 ^ 2, u2, -1 + *(4, u, xi))) + *(3, St, u3)), coshGu4, exp(*(B1, u4) + *(B2, u4))))

    nothing
end


#= tilde, hat, etc, definitions

We use these macros as shorthand notation. For instance

  @tilde_outer("B1")
  @hat_outer("B2")

should expand to

  B1t = B1_x - (Fx + xi_x) * B1p
  B2h = B2_y - (Fy + xi_y) * B2p

etc.

=#
macro tilde_outer(fname::String)
    ft    = Symbol(fname, "t")
    f_x   = Symbol(fname, "_x")
    fp    = Symbol(fname, "p")
    return esc( :($ft = $f_x - (Fx + xi_x) * $fp) )
end
macro hat_outer(fname::String)
    fh    = Symbol(fname, "h")
    f_y   = Symbol(fname, "_y")
    fp    = Symbol(fname, "p")
    return esc( :($fh = $f_y - (Fy + xi_y) * $fp) )
end
macro bar_outer(fname::String)
    fb    = Symbol(fname, "b")
    f_xx  = Symbol(fname, "_xx")
    fpp   = Symbol(fname, "pp")
    fp_x  = Symbol(fname, "p_x")
    return esc( :($fb = $f_xx + (Fx + xi_x) * ( -2*($fp_x) + (Fx + xi_x) * ($fpp) )) )
end
macro star_outer(fname::String)
    fs    = Symbol(fname, "s")
    f_yy  = Symbol(fname, "_yy")
    fpp   = Symbol(fname, "pp")
    fp_y  = Symbol(fname, "p_y")
    return esc( :($fs = $f_yy + (Fy + xi_y) * ( -2*($fp_y) + (Fy + xi_y) * ($fpp) )) )
end
macro cross_outer(fname::String)
    fc    = Symbol(fname, "c")
    f_xy  = Symbol(fname, "_xy")
    fpp   = Symbol(fname, "pp")
    fp_x  = Symbol(fname, "p_x")
    fp_y  = Symbol(fname, "p_y")
    return esc( :($fc = $f_xy  - (Fx + xi_x) * ($fp_y) -
                  (Fy + xi_y) * ( $fp_x - (Fx + xi_x) * ($fpp) ) ) )
end


# assuming
# (A d_uu + B d_u + C Id) f = -S

function S_eq_coeff!(ABCS::Vector, vars::Tuple, ::Outer)
    (phi0, u, xi, B1, B1p, B2, B2p, G, Gp, phi, phip) = vars

    ABCS[1] = *(6, u ^ 4)
    ABCS[2] = *(12, u ^ 3)
    ABCS[3] = Gp ^ 2 + *(3, B2p ^ 2) + *(4, phip ^ 2) + *(B1p ^ 2, cosh(G) ^ 2)
    ABCS[4] = 0

    nothing
end

# this is a coupled equation for Fx and Fy. the notation used is
#
# ( A11 d_uu Fx + A12 d_uu Fy + B11 d_u Fx + B12 d_u Fy + C11 Fx + C12 Fy ) = -S1
# ( A21 d_uu Fx + A22 d_uu Fy + B21 d_u Fx + B22 d_u Fy + C21 Fx + C22 Fy ) = -S2

function Fxy_eq_coeff!(AA::Matrix, BB::Matrix, CC::Matrix, SS::Vector, vars::Tuple, ::Outer)
    (
        phi0, u, xi, xi_x, xi_y,
        B1     ,    B2     ,    G      ,    phi    ,    S      ,
        B1p    ,    B2p    ,    Gp     ,    phip   ,    Sp     ,
        B1pp   ,    B2pp   ,    Gpp    ,                Spp    ,
        B1_x   ,    B2_x   ,    G_x    ,    phi_x  ,    S_x    ,
        B1_y   ,    B2_y   ,    G_y    ,    phi_y  ,    S_y    ,
        B1p_x  ,    B2p_x  ,    Gp_x   ,                Sp_x   ,
        B1p_y  ,    B2p_y  ,    Gp_y   ,                Sp_y
    ) = vars

    expB1   = exp(B1)
    sinhG   = sinh(G)
    sinh2G  = sinh(*(2, G))
    coshG   = cosh(G)
    cosh2G  = cosh(*(2, G))
    coshGsq = cosh(G)^2


    AA[1,1] = *(2, S ^ 2, u ^ 4, expB1)
    AA[1,2] = 0
    AA[2,1] = 0
    AA[2,2] = *(2, S ^ 2, u ^ 4)

    BB[1,1] = *(2, S, u ^ 2, *(-1, Sp) + *(-1, B2p, S) + *(2, S, u) + *(-1, B1p, S, coshGsq), expB1)
    BB[1,2] = *(2, S ^ 2, u ^ 2, Gp + *(B1p, coshG, sinhG))
    BB[2,1] = *(2, S ^ 2, u ^ 2, Gp + *(-1, B1p, coshG, sinhG), expB1)
    BB[2,2] = *(S, u ^ 2, *(-2, Sp) + *(S, B1p + *(-2, B2p) + *(4, u)) + *(B1p, S, cosh2G))

    CC[1,1] = *(-2, *(4, Sp ^ 2) + *(-1, B2pp, S ^ 2) + *(-1, Gp ^ 2, S ^ 2) + *(-4, S, Spp) + *(-4, phip ^ 2, S ^ 2) + *(-3, B2p ^ 2, S ^ 2) + *(-1, S, coshGsq, *(B1pp, S) + *(S, B1p ^ 2) + *(3, B1p, Sp)) + *(-3, B2p, S, Sp) + *(-1, B1p, Gp, S ^ 2, sinh2G), expB1)

    CC[1,2] = *(-1, S, *(*(B1pp, S) + *(-1, S, B1p ^ 2) + *(3, B1p, Sp), sinh2G) + *(2, Gpp, S) + *(6, Gp, Sp) + *(-2, B1p, Gp, S) + *(2, B1p, Gp, S, cosh2G))

    CC[2,1] = *(S, *(*(B1pp, S) + *(S, B1p ^ 2) + *(3, B1p, Sp), sinh2G) + *(-6, Gp, Sp) + *(-2, Gpp, S) + *(-2, B1p, Gp, S) + *(2, B1p, Gp, S, cosh2G), expB1)

    CC[2,2] = *(-8, Sp ^ 2) + *(2, B2pp, S ^ 2) + *(2, Gp ^ 2, S ^ 2) + *(6, B2p ^ 2, S ^ 2) + *(8, S, Spp) + *(8, phip ^ 2, S ^ 2) + *(2, S, coshGsq, *(S, B1p ^ 2) + *(-1, B1pp, S) + *(-3, B1p, Sp)) + *(6, B2p, S, Sp) + *(-2, B1p, Gp, S ^ 2, sinh2G)

    SS[1] = *(S ^ 2, *(2, Gp_y) + *(B1p_y + *(xi_y, B1p ^ 2) + *(-1, B1_y, B1p) + *(-1, B1pp, xi_y), sinh2G) + *(-2, B1_y, Gp) + *(-2, Gpp, xi_y) + *(2, B1p, Gp, xi_y) + *(2, B1p, G_y + *(-1, Gp, xi_y), cosh2G)) + *(*(S, *(-8, Sp_x) + *(-1, B2p + *(B1p, coshGsq), *(6, S_x) + *(-6, Sp, xi_x)) + *(8, Spp, xi_x)) + *(S ^ 2, *(-2, B2p_x) + *(-8, phip, phi_x + *(-1, phip, xi_x)) + *(-6, B2p, B2_x + *(-1, B2p, xi_x)) + *(-2, coshGsq, B1p_x + *(-1, B1pp, xi_x)) + *(-2, G_x + *(-1, Gp, xi_x), Gp + *(B1p, sinh2G)) + *(2, B2pp, xi_x) + *(2, B1p, coshGsq, *(-1, B1_x) + *(B1p, xi_x))) + *(8, Sp, S_x + *(-1, Sp, xi_x)), expB1) + *(S, Gp + *(B1p, coshG, sinhG), *(6, S_y) + *(-6, Sp, xi_y))

    SS[2] = *(S, *(-8, Sp_y) + *(-1, B2p + *(-1, B1p, coshGsq), *(6, S_y) + *(-6, Sp, xi_y)) + *(8, Spp, xi_y)) + *(2, S ^ 2, *(-1, B2p_y) + *(B2pp, xi_y) + *(xi_y, Gp ^ 2) + *(coshGsq, B1p_y + *(xi_y, B1p ^ 2) + *(-1, B1_y, B1p) + *(-1, B1pp, xi_y)) + *(-1, G_y, Gp) + *(-4, phi_y, phip) + *(-3, B2_y, B2p) + *(3, xi_y, B2p ^ 2) + *(4, xi_y, phip ^ 2) + *(B1p, G_y + *(-1, Gp, xi_y), sinh2G)) + *(8, Sp, S_y + *(-1, Sp, xi_y)) + *(S, *(S, *(2, Gp_x) + *(-1, B1p_x + *(B1_x, B1p) + *(-1, xi_x, B1pp + B1p ^ 2), sinh2G) + *(-2, xi_x, Gpp + *(B1p, Gp)) + *(2, B1_x, Gp) + *(-2, B1p, G_x + *(-1, Gp, xi_x), cosh2G)) + *(Gp + *(-1, B1p, coshG, sinhG), *(6, S_x) + *(-6, Sp, xi_x)), expB1)

    nothing
end


function Sd_eq_coeff!(ABCS::Vector, vars::Tuple, ::Outer)
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

    @tilde_outer("B1")
    @tilde_outer("B2")
    @tilde_outer("G")
    @tilde_outer("phi")
    @tilde_outer("S")
    @tilde_outer("Fx")
    @tilde_outer("Fy")

    @hat_outer("B1")
    @hat_outer("B2")
    @hat_outer("G")
    @hat_outer("phi")
    @hat_outer("S")
    @hat_outer("Fx")
    @hat_outer("Fy")

    @bar_outer("B1")
    @bar_outer("B2")
    @bar_outer("G")
    @bar_outer("phi")
    @bar_outer("S")

    @star_outer("B1")
    @star_outer("B2")
    @star_outer("G")
    @star_outer("phi")
    @star_outer("S")

    @tilde_outer("Fxp")
    @tilde_outer("Fyp")

    @hat_outer("Fxp")
    @hat_outer("Fyp")

    @cross_outer("B2")
    @cross_outer("G")
    @cross_outer("S")

    VV  = -3 - 3/2 * phi0 * u * (1 - 2*u*xi) + u*u*u*u * UU(phi, potential)

    # VVp = -3 * phi0*u * (1 - u*xi + phi0*phi0 * u*u * phi) + u*u*u * Up(phi)

    expB1   = exp(B1)
    expB2   = exp(B2)
    sinh2G  = sinh(*(2, G))
    cosh2G  = cosh(*(2, G))
    coshGsq = cosh(G)^2
    coshG   = cosh(G)
    sinhG   = sinh(G)


    ABCS[1] = 0

    ABCS[2] = *(-12, S ^ 3, u ^ 2, expB1)

    ABCS[3] = *(24, Sp, S ^ 2, expB1)

    ABCS[4] = *(*(S, *(*(-8, Gh, St) + *(-8, Gt, Sh), coshG) + *(*(-16, Sc) + *(2, Sh, Fxp + *(-4, B2t)) + *(2, Sp, *(4, Fxh) + *(4, Fyt) + *(8, xi_xy)) + *(2, St, Fyp + *(-4, B2h)), sinhG)) + *(S ^ 2, *(*(-4, Gc) + *(-2, Gh, B1t + B2t + *(-1, Fxp)) + *(2, Gp, Fxh + Fyt + *(2, xi_xy)) + *(2, Gt, B1h + Fyp + *(-1, B2h)), coshG) + *(*(-4, B2c) + *(2, Fxph) + *(2, Fypt) + *(-8, phih, phit) + *(-4, Gh, Gt) + *(2, B2h, Fxp + *(-4, B2t)) + *(2, B2p, Fxh + Fyt + *(2, xi_xy)) + *(2, Fyp, B2t + *(-1, Fxp)), sinhG)) + *(8, Sh, St, sinhG), exp(B1 + B2)) + *(*(S, *(*(8, Ss) + *(-2, Sh, Fyp + *(-4, B2h) + *(4, B1h)) + *(-2, Sp, *(4, Fyh) + *(4, xi_yy)), coshG) + *(8, Gh, Sh, sinhG)) + *(S ^ 2, *(*(2, Gs) + *(-2, Gp, Fyh + xi_yy) + *(2, Gh, B2h + *(-1, Fyp) + *(-2, B1h)), sinhG) + *(Fyp ^ 2 + *(-2, B1s) + *(-2, Fyph) + *(2, B2s) + *(2, B1h ^ 2) + *(2, Gh ^ 2) + *(4, B2h ^ 2) + *(4, phih ^ 2) + *(Fyp, *(-2, B2h) + *(2, B1h)) + *(-2, B1h, B2h) + *(2, B1p + *(-1, B2p), Fyh + xi_yy), coshG)) + *(-4, Sh ^ 2, coshG), expB2) + *(*(S, *(*(8, Sb) + *(-8, Fxt, Sp) + *(-8, Sp, xi_xx) + *(-2, Fxp, St) + *(8, B1t, St) + *(8, B2t, St), coshG) + *(8, Gt, St, sinhG)) + *(S ^ 2, *(*(2, Gb) + *(-2, Gp, Fxt + xi_xx) + *(2, Gt, B2t + *(-1, Fxp) + *(2, B1t)), sinhG) + *(Fxp ^ 2 + *(-2, Fxpt) + *(2, B1b) + *(2, B2b) + *(2, B1t ^ 2) + *(2, Gt ^ 2) + *(4, B2t ^ 2) + *(4, phit ^ 2) + *(-1, Fxp, *(2, B1t) + *(2, B2t)) + *(-2, B1p + B2p, Fxt + xi_xx) + *(2, B1t, B2t), coshG)) + *(-4, St ^ 2, coshG), exp(B2 + *(2, B1))) + *(8, VV, S ^ 4, expB1)

    nothing
end


function B2d_eq_coeff!(ABCS::Vector, vars::Tuple, ::Outer)
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

    @tilde_outer("B1")
    @tilde_outer("B2")
    @tilde_outer("G")
    @tilde_outer("phi")
    @tilde_outer("S")
    @tilde_outer("Fx")
    @tilde_outer("Fy")

    @hat_outer("B1")
    @hat_outer("B2")
    @hat_outer("G")
    @hat_outer("phi")
    @hat_outer("S")
    @hat_outer("Fx")
    @hat_outer("Fy")

    @bar_outer("B1")
    @bar_outer("B2")
    @bar_outer("G")
    @bar_outer("phi")
    @bar_outer("S")

    @star_outer("B1")
    @star_outer("B2")
    @star_outer("G")
    @star_outer("phi")
    @star_outer("S")

    @tilde_outer("Fxp")
    @tilde_outer("Fyp")

    @hat_outer("Fxp")
    @hat_outer("Fyp")

    @cross_outer("B2")
    @cross_outer("G")
    @cross_outer("S")


    expB1   = exp(B1)
    expB2   = exp(B2)
    sinh2G  = sinh(*(2, G))
    cosh2G  = cosh(*(2, G))
    coshGsq = cosh(G)^2
    coshG   = cosh(G)
    sinhG   = sinh(G)


    ABCS[1] = 0

    ABCS[2] = *(-12, S ^ 4, u ^ 2, expB1)

    ABCS[3] = *(18, Sp, S ^ 3, expB1)

    ABCS[4] = *(*(S, *(*(2, Gh, St) + *(2, Gt, Sh), coshG) + *(-1, *(-4, Sc) + *(-4, Fxp, Sh) + *(-4, Fyp, St) + *(2, Fxh, Sp) + *(2, Fyt, Sp) + *(4, B2h, St) + *(4, B2t, Sh) + *(4, Sp, xi_xy), sinhG)) + *(S ^ 2, *(*(4, Gc) + *(-2, Gp, Fxh + Fyt + *(2, xi_xy)) + *(-2, Gt, B1h + Fyp + *(2, B2h)) + *(2, Gh, B1t + *(-1, Fxp) + *(-2, B2t)), coshG) + *(*(-8, B2c) + *(-2, Fxph) + *(-2, Fypt) + *(2, Fyp, Fxp + *(2, B2t)) + *(4, B2h, Fxp + *(-1, B2t)) + *(4, B2p, Fxh + Fyt + *(2, xi_xy)) + *(4, Gh, Gt) + *(8, phih, phit), sinhG)) + *(-8, Sh, St, sinhG), exp(B1 + B2)) + *(*(S, *(*(-2, Ss) + *(2, Sh, B1h + *(-2, Fyp) + *(2, B2h)) + *(2, Sp, Fyh + xi_yy), coshG) + *(-2, Gh, Sh, sinhG)) + *(S ^ 2, *(*(-2, Gs) + *(2, Gh, Fyp + *(2, B1h) + *(2, B2h)) + *(2, Gp, Fyh + xi_yy), sinhG) + *(*(-1, Fyp ^ 2) + *(-4, phih ^ 2) + *(-2, B1h ^ 2) + *(-2, Gh ^ 2) + *(2, B1s) + *(2, Fyph) + *(2, B2h ^ 2) + *(4, B2s) + *(-1, Fyp, *(2, B1h) + *(4, B2h)) + *(-1, Fyh + xi_yy, *(2, B1p) + *(4, B2p)) + *(-4, B1h, B2h), coshG)) + *(4, Sh ^ 2, coshG), expB2) + *(*(S, *(*(-2, Sb) + *(-4, Fxp, St) + *(-2, B1t, St) + *(2, Fxt, Sp) + *(2, Sp, xi_xx) + *(4, B2t, St), coshG) + *(-2, Gt, St, sinhG)) + *(S ^ 2, *(*(-2, Gb) + *(2, Gp, Fxt + xi_xx) + *(2, Gt, Fxp + *(-2, B1t) + *(2, B2t)), sinhG) + *(*(-1, Fxp ^ 2) + *(-4, phit ^ 2) + *(-2, B1b) + *(-2, B1t ^ 2) + *(-2, Gt ^ 2) + *(2, Fxpt) + *(2, B2t ^ 2) + *(4, B2b) + *(Fxp, *(-4, B2t) + *(2, B1t)) + *(2, B1p + *(-2, B2p), Fxt + xi_xx) + *(4, B1t, B2t), coshG)) + *(4, St ^ 2, coshG), exp(B2 + *(2, B1))) + *(18, B2p, Sd, S ^ 3, expB1)

    nothing
end


# this is another coupled equation, for B1d and Gd. the notation used is
#
# ( A11 d_uu B1d + A12 d_uu Gd + B11 d_u B1d + B12 d_u Gd + C11 B1d + C12 Gd ) = -S1
# ( A21 d_uu B1d + A22 d_uu Gd + B21 d_u B1d + B22 d_u Gd + C21 B1d + C22 Gd ) = -S2

function B1dGd_eq_coeff!(AA::Matrix, BB::Matrix, CC::Matrix, SS::Vector, vars::Tuple, ::Outer)
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

    @tilde_outer("B1")
    @tilde_outer("B2")
    @tilde_outer("G")
    @tilde_outer("phi")
    @tilde_outer("S")
    @tilde_outer("Fx")
    @tilde_outer("Fy")

    @hat_outer("B1")
    @hat_outer("B2")
    @hat_outer("G")
    @hat_outer("phi")
    @hat_outer("S")
    @hat_outer("Fx")
    @hat_outer("Fy")

    @bar_outer("B1")
    @bar_outer("B2")
    @bar_outer("G")
    @bar_outer("phi")
    @bar_outer("S")

    @star_outer("B1")
    @star_outer("B2")
    @star_outer("G")
    @star_outer("phi")
    @star_outer("S")

    @tilde_outer("Fxp")
    @tilde_outer("Fyp")

    @hat_outer("Fxp")
    @hat_outer("Fyp")

    @cross_outer("B2")
    @cross_outer("G")
    @cross_outer("S")


    expB1   = exp(B1)
    expB2   = exp(B2)
    sinh2G  = sinh(*(2, G))
    cosh2G  = cosh(*(2, G))
    coshGsq = cosh(G)^2
    coshG   = cosh(G)
    tanhG   = tanh(G)
    sinhG   = sinh(G)
    sechG   = sech(G)


    AA[1,1] = 0
    AA[1,2] = 0
    AA[2,1] = 0
    AA[2,2] = 0


    BB[1,1] = *(-12, S ^ 4, u ^ 2, expB1)

    BB[1,2] = 0

    BB[2,1] = 0

    BB[2,2] = *(-12, S ^ 4, u ^ 2, expB1)


    CC[1,1] = *(6, S ^ 3, *(3, Sp) + *(2, Gp, S, tanhG), expB1)

    CC[1,2] = *(12, B1p, S ^ 4, expB1, tanhG)

    CC[2,1] = *(-6, B1p, S ^ 4, expB1, sinh2G)

    CC[2,2] = *(18, Sp, S ^ 3, expB1)


    SS[1] = *(3, *(-4, Sh ^ 2) + *(S, *(2, Ss) + *(-2, Sp, Fyh + xi_yy) + *(2, Sh, B2h + *(2, Fyp))) + *(S ^ 2, Fyp ^ 2 + *(-2, Fyph) + *(2, B2s) + *(4, B2h ^ 2) + *(4, phih ^ 2) + *(-2, B2h, Fyp) + *(-2, B2p, Fyh + xi_yy)), expB2, sechG) + *(3, *(4, St ^ 2) + *(-1, S ^ 2, Fxp ^ 2 + *(-2, Fxpt) + *(2, B2b) + *(4, B2t ^ 2) + *(4, phit ^ 2) + *(-2, B2p, Fxt + xi_xx) + *(-2, B2t, Fxp)) + *(2, S, *(-1, Sb) + *(Fxt, Sp) + *(Sp, xi_xx) + *(-1, B2t, St) + *(-2, Fxp, St)), exp(B2 + *(2, B1)), sechG) + *(6, S, *(Gt, Sh + *(S, B2h + *(-1, Fyp))) + *(-1, Gh, St) + *(Fxp, Gh, S) + *(Fyt, Gp, S) + *(-1, B2t, Gh, S) + *(-1, Fxh, Gp, S), exp(B1 + B2), sechG) + *(18, B1p, Sd, S ^ 3, expB1)

    SS[2] = *(3, *(4, Sh ^ 2) + *(-1, S, *(2, Ss) + *(-2, Sp, Fyh + xi_yy) + *(2, Sh, B2h + *(2, Fyp))) + *(-1, S ^ 2, Fyp ^ 2 + *(-2, Fyph) + *(2, B2s) + *(4, B2h ^ 2) + *(4, phih ^ 2) + *(-2, B2h, Fyp) + *(-2, B2p, Fyh + xi_yy)), expB2, sinhG) + *(3, *(4, St ^ 2) + *(-1, S ^ 2, Fxp ^ 2 + *(-2, Fxpt) + *(2, B2b) + *(4, B2t ^ 2) + *(4, phit ^ 2) + *(-2, B2p, Fxt + xi_xx) + *(-2, B2t, Fxp)) + *(2, S, *(-1, Sb) + *(Fxt, Sp) + *(Sp, xi_xx) + *(-1, B2t, St) + *(-2, Fxp, St)), exp(B2 + *(2, B1)), sinhG) + *(6, *(S, *(2, Sc) + *(Sh, B2t + *(-1, B1t) + *(2, Fxp)) + *(St, B1h + B2h + *(2, Fyp)) + *(-1, Sp, Fxh + Fyt + *(2, xi_xy))) + *(S ^ 2, *(-1, Fxph) + *(-1, Fypt) + *(2, B2c) + *(B1p, Fxh + *(-1, Fyt)) + *(B1t, Fyp) + *(B1h + *(-1, Fyp), B2t + *(-1, Fxp)) + *(-1, B2h, B1t + Fxp + *(-4, B2t)) + *(-1, B2p, Fxh + Fyt + *(2, xi_xy)) + *(4, phih, phit)) + *(-4, Sh, St), coshG, exp(B1 + B2)) + *(18, Gp, Sd, S ^ 3, expB1)

    nothing
end


function phid_eq_coeff!(ABCS::Vector, vars::Tuple, ::Outer)
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

    @tilde_outer("B1")
    @tilde_outer("B2")
    @tilde_outer("G")
    @tilde_outer("phi")
    @tilde_outer("S")
    @tilde_outer("Fx")
    @tilde_outer("Fy")

    @hat_outer("B1")
    @hat_outer("B2")
    @hat_outer("G")
    @hat_outer("phi")
    @hat_outer("S")
    @hat_outer("Fx")
    @hat_outer("Fy")

    # @bar_outer("B1")
    # @bar_outer("B2")
    # @bar_outer("G")
    @bar_outer("phi")
    # @bar_outer("S")

    # @star_outer("B1")
    # @star_outer("B2")
    # @star_outer("G")
    @star_outer("phi")
    # @star_outer("S")

    # @tilde_outer("Fxp")
    # @tilde_outer("Fyp")

    # @hat_outer("Fxp")
    # @hat_outer("Fyp")

    # @cross_outer("B2")
    # @cross_outer("G")
    # @cross_outer("S")
    @cross_outer("phi")

    # VV  = -3 - 3/2 * phi0 * u * (1 - 2*u*xi) + u*u*u*u * UU(phi, potential)

    VVp = -3 * phi0*u * (1 - u*xi + phi0*phi0 * u*u * phi) + u*u*u * UUp(phi, potential)



    expB1   = exp(B1)
    expB2   = exp(B2)
    sinh2G  = sinh(*(2, G))
    cosh2G  = cosh(*(2, G))
    coshGsq = cosh(G)^2
    coshG   = cosh(G)
    sinhG   = sinh(G)


    ABCS[1] = 0

    ABCS[2] = *(-8, S ^ 3, u ^ 2, expB1)

    ABCS[3] = *(12, Sp, S ^ 2, expB1)

    ABCS[4] = *(*(*(4, phih, Sh + *(-1, S, B1h + Fyp + *(-1, B2h))) + *(4, S, phis + *(-1, phip, Fyh + xi_yy)), coshG) + *(4, Gh, phih, S, sinhG), expB2) + *(*(4, S, *(phib + *(phit, B1t + B2t + *(-1, Fxp)) + *(-1, phip, Fxt + xi_xx), coshG) + *(Gt, phit, sinhG)) + *(4, phit, St, coshG), exp(B2 + *(2, B1))) + *(*(-1, *(4, phih, St) + *(4, phit, Sh), sinhG) + *(S, *(-4, Gh, phit) + *(-4, Gt, phih), coshG) + *(S, *(-8, phic) + *(4, phih, Fxp + *(-1, B2t)) + *(4, phip, Fxh + Fyt + *(2, xi_xy)) + *(4, phit, Fyp + *(-1, B2h)), sinhG), exp(B1 + B2)) + *(-4, S ^ 2, *(S, VVp) + *(-3, phip, Sd), expB1)

    nothing
end


function A_eq_coeff!(ABCS::Vector, vars::Tuple, ::Outer)
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

    @tilde_outer("B1")
    @tilde_outer("B2")
    @tilde_outer("G")
    @tilde_outer("phi")
    @tilde_outer("S")
    @tilde_outer("Fx")
    @tilde_outer("Fy")

    @hat_outer("B1")
    @hat_outer("B2")
    @hat_outer("G")
    @hat_outer("phi")
    @hat_outer("S")
    @hat_outer("Fx")
    @hat_outer("Fy")

    @bar_outer("B1")
    @bar_outer("B2")
    @bar_outer("G")
    @bar_outer("phi")
    @bar_outer("S")

    @star_outer("B1")
    @star_outer("B2")
    @star_outer("G")
    @star_outer("phi")
    @star_outer("S")

    @tilde_outer("Fxp")
    @tilde_outer("Fyp")

    @hat_outer("Fxp")
    @hat_outer("Fyp")

    @cross_outer("B2")
    @cross_outer("G")
    @cross_outer("S")
    @cross_outer("phi")

    VV  = -3 - 3/2 * phi0 * u * (1 - 2*u*xi) + u*u*u*u * UU(phi, potential)

    # VVp = -3 * phi0*u * (1 - u*xi + phi0*phi0 * u*u * phi) + u*u*u * Up(phi)

    expB1   = exp(B1)
    expB2   = exp(B2)
    sinh2G  = sinh(*(2, G))
    cosh2G  = cosh(*(2, G))
    coshGsq = cosh(G)^2
    coshG   = cosh(G)
    sinhG   = sinh(G)


    ABCS[1] = *(6, S ^ 4, u ^ 4, expB1)

    ABCS[2] = *(12, S ^ 4, u ^ 3, expB1)

    ABCS[3] = 0

    ABCS[4] = *(*(2, S ^ 4, *(-4, VV) + *(3, Gd, Gp) + *(9, B2d, B2p) + *(12, phid, phip) + *(3, B1d, B1p, coshGsq)) + *(-72, Sd, Sp, S ^ 2), expB1) + *(*(S, *(*(-24, Ss) + *(24, Sh, B1h + *(-1, B2h)) + *(24, Sp, Fyh + xi_yy), coshG) + *(-24, Gh, Sh, sinhG)) + *(S ^ 2, *(*(-6, Gs) + *(-6, B2h, Gh) + *(6, Fyh, Gp) + *(6, Gp, xi_yy) + *(12, B1h, Gh), sinhG) + *(*(-12, B2h ^ 2) + *(-12, phih ^ 2) + *(-6, B2s) + *(-6, B1h ^ 2) + *(-6, Gh ^ 2) + *(3, Fyp ^ 2) + *(6, B1s) + *(-6, B1p + *(-1, B2p), Fyh + xi_yy) + *(6, B1h, B2h), coshG)) + *(12, Sh ^ 2, coshG), expB2) + *(*(S, *(*(-24, Sb) + *(-24, B1t, St) + *(-24, B2t, St) + *(24, Fxt, Sp) + *(24, Sp, xi_xx), coshG) + *(-24, Gt, St, sinhG)) + *(S ^ 2, *(*(-12, B2t ^ 2) + *(-12, phit ^ 2) + *(-6, B1b) + *(-6, B2b) + *(-6, B1t ^ 2) + *(-6, Gt ^ 2) + *(3, Fxp ^ 2) + *(-6, B1t, B2t) + *(6, B1p + B2p, Fxt + xi_xx), coshG) + *(-1, *(6, Gb) + *(-6, Gp, Fxt + xi_xx) + *(6, Gt, B2t + *(2, B1t)), sinhG)) + *(12, St ^ 2, coshG), exp(B2 + *(2, B1))) + *(*(S ^ 2, *(*(12, Gc) + *(-6, Gp, Fxh + Fyt + *(2, xi_xy)) + *(6, Gh, B1t + B2t) + *(6, Gt, B2h + *(-1, B1h)), coshG) + *(*(12, B2c) + *(-6, B2p, Fxh + Fyt + *(2, xi_xy)) + *(-6, Fxp, Fyp) + *(12, Gh, Gt) + *(24, B2h, B2t) + *(24, phih, phit), sinhG)) + *(24, S, *(*(Gh, St) + *(Gt, Sh), coshG) + *(*(2, Sc) + *(B2h, St) + *(B2t, Sh) + *(-1, Sp, Fxh + Fyt + *(2, xi_xy)), sinhG)) + *(-24, Sh, St, sinhG), exp(B1 + B2))

    nothing
end

function xi_t_eq_coeff(vars::Tuple, ::Outer)
    (
        kappa, xi, xi_x, xi_y, xi_xx, xi_yy, xi_xy,
        B1   , B2   , G   ,  S    , Fx    , Fy    , Sd ,  B1d  , B2d  , Gd,  phid, A, Ad, Fxd, Fyd,
        B1p  , B2p  , Gp  ,  Sp   , Fxp   , Fyp   , Sdp,  B1dp , B2dp , Gdp,       Ap, Adp, Fxdp, Fydp,
        B1pp , B2pp , Gpp ,  Spp  , Fxpp  , Fypp  ,                                App,
        B1_x , B2_x , G_x ,  S_x  , Fx_x  , Fy_x  , Sd_x, B1d_x, B2d_x, Gd_x,      A_x, Ad_x, Fxd_x, Fyd_x,
        B1_y , B2_y , G_y ,  S_y  , Fx_y  , Fy_y  , Sd_y, B1d_y, B2d_y, Gd_y,      A_y, Ad_y, Fxd_y, Fyd_y,
        B1p_x, B2p_x, Gp_x,  Sp_x , Fxp_x , Fyp_x ,                                Ap_x,
        B1p_y, B2p_y, Gp_y,  Sp_y , Fxp_y , Fyp_y ,                                Ap_y,
                                                                                   A_xx,
                                                                                   A_yy,
                                                                                   A_xy,
    ) = vars

    @tilde_outer("B1")
    @tilde_outer("B2")
    @tilde_outer("G")
    @tilde_outer("S")
    @tilde_outer("Fx")
    @tilde_outer("Fy")
    @tilde_outer("Sd")
    @tilde_outer("B1d")
    @tilde_outer("B2d")
    @tilde_outer("Gd")
    @tilde_outer("A")

    @hat_outer("B1")
    @hat_outer("B2")
    @hat_outer("G")
    @hat_outer("S")
    @hat_outer("Fx")
    @hat_outer("Fy")
    @hat_outer("Sd")
    @hat_outer("B1d")
    @hat_outer("B2d")
    @hat_outer("Gd")
    @hat_outer("A")

    @bar_outer("A")

    @star_outer("A")

    @tilde_outer("B1p")
    @tilde_outer("B2p")
    @tilde_outer("Gp")
    @tilde_outer("Sp")
    @tilde_outer("Fxp")
    @tilde_outer("Fyp")

    @hat_outer("B1p")
    @hat_outer("B2p")
    @hat_outer("Gp")
    @hat_outer("Sp")
    @hat_outer("Fxp")
    @hat_outer("Fyp")

    @cross_outer("A")

    ## FIXME
    Fydh = 0
    Fyd  = 0
    Fxdt = 0
    Fydt = 0
    Fxdp = 0
    Fydp = 0
    Fxd  = 0
    Fxdh = 0

##    expB1   = exp(B1)
##    expB2   = exp(B2)
##    sinh2G  = sinh(*(2, G))
##    cosh2G  = cosh(*(2, G))
##    coshGsq = cosh(G)^2
##    coshG   = cosh(G)
##    sinhG   = sinh(G)

    κ = kappa


    x0 = 2*B1
    x1 = exp(B2 + x0)
    x2 = cosh(G)
    x3 = 2*x2
    x4 = exp(B1 + B2)
    x5 = sinh(G)
    x6 = x4*x5
    x7 = 4*S
    x8 = exp(B2)
    x9 = S*x8
    x10 = x3*x9
    x11 = exp(B1)
    x12 = S*x11
    x13 = Gt*x12
    x14 = -x13
    x15 = 2*Fx
    x16 = Gp*x15
    x17 = x12*x16
    x18 = 2*xi_x
    x19 = Gp*x18
    x20 = x12*x19
    x21 = Fyp*S
    x22 = Fy*Sp
    x23 = Sp*xi_y
    x24 = 2*B2p
    x25 = S*x24
    x26 = x25*xi_y
    x27 = B2h*S
    x28 = Sh + x27
    x29 = Fy*x25 + x21 - x22 - x23 + x26 + x28
    x30 = 2*x5
    x31 = Fy + xi_y
    x32 = Gp*x31
    x33 = S*(Gh + 2*x32)
    x34 = 2*B1p
    x35 = Fx*S
    x36 = S*x34
    x37 = B1t*S
    x38 = B2t*S
    x39 = St + x38
    x40 = x37 + x39
    x41 = Fxp*S
    x42 = Fx*Sp
    x43 = Sp*xi_x
    x44 = x25*xi_x
    x45 = x41 - x42 - x43 + x44
    x46 = -Sh
    x47 = B1h*S
    x48 = x36*xi_y
    x49 = 2*Sd
    x50 = B1p*S
    x51 = 2*Gh
    x52 = 2*Fy
    x53 = Gp*x52
    x54 = Gph*S
    x55 = 2*Fyh
    x56 = Gp*S
    x57 = 2*x21
    x58 = 2*xi_yy
    x59 = 2*xi_y
    x60 = Sh*x59
    x61 = Sph*x11
    x62 = 2*Fxh
    x63 = Sp*x11
    x64 = 2*Sh
    x65 = Fxp*x11
    x66 = 2*x12
    x67 = Spt*x11
    x68 = Fy^2
    x69 = Gp*Sp
    x70 = 2*Fyp
    x71 = x11*x70
    x72 = xi_y^2
    x73 = Sp*x72
    x74 = 4*xi_xy
    x75 = x47*x59
    x76 = x27*x59
    x77 = Fy*x21
    x78 = 4*Gp
    x79 = Gp*x59
    x80 = Fy*Gpp
    x81 = 4*xi_y
    x82 = S*x81
    x83 = Gp*x81
    x84 = 2*S
    x85 = Gpp*x84
    x86 = 2*B2h
    x87 = x11*x86
    x88 = 2*x41
    x89 = x11*x88
    x90 = B2ph*x12
    x91 = B2pt*x12
    x92 = 2*x22
    x93 = B2t*x11
    x94 = 2*x23
    x95 = Spp*x52
    x96 = Fx*x11
    x97 = Fypp*x12
    x98 = Spp*x59
    x99 = Fxpp*x12
    x100 = x11*xi_x
    x101 = B2pp*Fy
    x102 = x11*x35
    x103 = 4*x102
    x104 = B2pp*x81
    x105 = 4*xi_x
    x106 = x105*x12
    x107 = x12*xi_x
    x108 = 2*Gt
    x109 = Gp*x12
    x110 = x108*x109
    x111 = Gh*xi_x
    x112 = exp(x0)
    x113 = St*x112
    x114 = S*x112
    x115 = Gpt*x114
    x116 = x108*x112
    x117 = 2*x112
    x118 = x117*x56
    x119 = Fx^2
    x120 = x112*x119
    x121 = xi_x^2
    x122 = Sp*x121
    x123 = x112*x37
    x124 = x112*x38
    x125 = Fx*Gp
    x126 = 4*x112
    x127 = x126*x41
    x128 = Gp^2
    x129 = Fy*x128
    x130 = x128*x81
    x131 = x112*x35
    x132 = x105*x131
    x133 = Gp*xi_x
    x134 = x112*x121
    x135 = Fx + xi_x
    x136 = B2p^2
    x137 = 4*x12
    x138 = x114*x135
    x139 = x23 + x28 + x57
    x140 = Fy*x43
    x141 = x43*xi_y
    x142 = Fxh*S
    x143 = Fy*St
    x144 = Fyt*S
    x145 = Sh*xi_x
    x146 = St*xi_y
    x147 = 2*xi_xy
    x148 = Fy*x38 + S*x147 + x142 + x143 + x144 + x145 + x146 + x27*xi_x + x38*xi_y
    x149 = x31^2
    x150 = Fy*x41
    x151 = x41*xi_y
    x152 = x5*x8
    x153 = 2*B1h
    x154 = B1ph*S
    x155 = B2ph*S
    x156 = Fypp*S
    x157 = B1p*x81
    x158 = Fy*x82
    x159 = Gh*x56
    x160 = B1pp*x84
    x161 = B2pp*x84
    x162 = x11*x51
    x163 = x11*x16
    x164 = Gph*x12
    x165 = 2*Gp*x11
    x166 = Gpt*x12
    x167 = x108*x11
    x168 = Spt*x112
    x169 = Fxp*x117
    x170 = x112*x84
    x171 = Sp*x117
    x172 = x12*x34
    x173 = Gt*xi_y
    x174 = x11*x38
    x175 = x11*x42
    x176 = 4*x11*x21
    x177 = Gpp*x81
    x178 = x11*x78
    x179 = x11*x43
    x180 = x128*x84
    x181 = B1p*x113
    x182 = 4*Fx
    x183 = x126*x50
    x184 = B1pt*x114
    x185 = B1t*x117
    x186 = B2pt*x114
    x187 = B2t*x117
    x188 = Fxpp*x114
    x189 = Spp*x112
    x190 = x125*x137
    x191 = B1p*Fy
    x192 = x133*x137
    x193 = x133*x81
    x194 = B1p*x123
    x195 = B1p*x124
    x196 = 6*B1p*x112*x41
    x197 = x108*x114
    x198 = B1p^2
    x199 = x198*x7
    x200 = x135^2
    x201 = 6*x50
    x202 = x2*x8
    x203 = 2*Fyd
    x204 = -A*Fyp + x203
    x205 = Sp*x31
    x206 = A*Gp
    x207 = Gd - x206/2
    x208 = 3*Sp
    x209 = A*B2p
    x210 = x149*x202
    x211 = 3*A*Spp/2 + 3*Ap*Sp/2 - 3*Sdp
    x212 = 2*Fxd
    x213 = -A*Fxp + x212
    x214 = Sp*x135
    x215 = x207*x5
    x216 = x1*x200
    x217 = 2*B2d
    x218 = -x209 + x217
    x219 = Fyp*x31
    x220 = Fyh + xi_yy
    x221 = S*(x219 + x220)
    x222 = 2*Gd - x206
    x223 = Fxp*x31
    x224 = Fxh + x223 + xi_xy
    x225 = S*x2
    x226 = x1*x2*x200
    x227 = 2*B1d
    x228 = B1p + B2p
    x229 = -A*x228 + x217 + x227
    x230 = 3*A
    x231 = x2/2
    x232 = A/2
    x233 = Ah/2 + Ap*x31/2
    x234 = S*x4
    x235 = A*B1p
    x236 = x227 - x235
    x237 = S*Sd
    x238 = Fyt + xi_xy
    x239 = Fxt + xi_xx
    x240 = Fx*Fxp + Fxp*xi_x + x239
    x241 = x11*x2
    x242 = x8*(x240*x241 - x5*(Fx*Fyp + Fyp*xi_x + x238))
    x243 = A*Sp
    x244 = Sd - x243/2
    x245 = x11*(6*x237 - 2*x242) - x202*(2*x219 + x55 + x58) + x6*(x147 + 2*x223 + x62)
    x246 = Sh + x205
    x247 = x246*x5
    x248 = St + x214
    x249 = x241*x248
    x250 = B2p*x31
    x251 = -B2h - x250
    x252 = B1t + B2t
    x253 = x11*(x2*(x135*x228 + x252) + x5*(Gt + x125 + x133))
    x254 = x2*(-Gh - x32) + x251*x5 + x253
    x255 = S*x254 - x247 + x249
    x256 = x255*x4
    x257 = 2*x214
    x258 = St - x257
    x259 = x11*x5
    x260 = x258*x259
    x261 = B1h + B1p*x31
    x262 = B2h + x250
    x263 = Gh + x32
    x264 = x263*x5
    x265 = B2p*x135 + B2t
    x266 = x259*x265
    x267 = Gp*x135 + Gt
    x268 = x241*x267
    x269 = x2*x261 - x2*x262 - x264 + x266 + x268
    x270 = S*x269 + x2*(-x205 + x46) + x260
    x271 = x270*x8
    x272 = x4*(x15 + x18)
    x273 = x8*(x52 + x59)
    x274 = S^3
    x275 = S^2
    x276 = 3*Sd
    x277 = x275*(Fyp*x135 + x238)
    x278 = x275*x4
    x279 = Ap*x135/2 + At/2
    x280 = At + x212
    x281 = S*x5
    x282 = Ah + x203
    x283 = S*(-B1h + B2h) + Sh
    x284 = x2*x207
    x285 = B1d - x235/2
    x286 = x207*x259
    x287 = Sdh + Sdp*x31 - Sp*x233 - x232*(Sph + Spp*x31)
    x288 = Sdp*x135
    x289 = Spp*x135
    x290 = B2dh + B2dp*x31 - B2p*x233 - x232*(B2ph + B2pp*x31)
    x291 = Gdh + Gdp*x31 - Gp*x233 - x232*(Gph + Gpp*x31)
    x292 = B2dp*x135 + B2dt - B2p*x279 - x232*(B2pp*x135 + B2pt)
    x293 = Gdp*x135 + Gdt - Gp*x279 - x232*(Gpp*x135 + Gpt)
    x294 = x207*x241

    axx = -S*x1*x3

    axy = x6*x7

    ayy = -x10

    bx = x4*(x3*(-x11*(x35*(x24 + x34) + x36*xi_x + x40 + x45) + x33) + x30*(x14 - x17 - x20 + x29))

    by = x8*(x3*(Fy*(Sp - x25 + x36) + x13 + x17 + x20 - x21 + x23 - x26 - x27 + x46 + x47 + x48) + x30*(x11*(x24*x35 + x39 + x45) - x33))

    cc = 6*x12*(S*Sdp + Sd*x50 + Sp*x49) - x152*(-B2h*x89 - Fxph*x66 + Fxt*x118 - Fy*x110 - Fypt*x66 - 2*Fyt*x63 - Gh*x17 + Gh*x57 - Gp*x111*x66 + Gp*x112*x122 + Gp*x60 + Gp*x73 - Gp*x75 + Gp*x76 + Gpp*x132 + Sh*x53 - St*x71 + x100*x95 + x100*x98 - x101*x103 - x101*x106 - x102*x104 - x102*x130 - x103*x129 - x104*x107 - x106*x129 - x107*x130 - x110*xi_y + x112*x19*x42 + x113*x16 + x113*x19 + x115*x15 + x115*x18 + x116*x41 + x116*x42 + x116*x43 + x118*xi_xx + x120*x69 + x120*x85 + x123*x16 + x123*x19 + x124*x16 + x124*x19 + x125*x127 + x127*x133 + x134*x85 - x135*x136*x137*x31 - x15*x61 - x15*x90 - x15*x97 - x18*x61 - x18*x90 - x18*x97 + x21*x83 + x22*x51 + x22*x79 + x23*x51 + x24*(-x11*(Fx*(x139 + x22) + x140 + x141 + x148 + x18*x21 + x41*x52 + x41*x59) + x138*(Gt + x16 + x19) + x31*x33) + x27*x53 + x34*(-x11*(Fx*x29 + Fy*x44 - x140 - x141 + x148 + x150 + x151 + x21*xi_x + x44*xi_y) + x138*(x108 + 3*x125 + 3*x133) - x149*x56) - x38*x71 + x42*x71 - x42*x87 + x43*x71 - x43*x87 - x47*x53 + x52*x54 - x52*x67 - x52*x91 - x52*x99 + x54*x59 + x55*x56 + x56*x58 - x59*x67 - x59*x91 - x59*x99 - x62*x63 - x63*x74 - x64*x65 + x65*x92 + x65*x94 + x68*x69 + x68*x85 + x72*x85 + x77*x78 + x80*x82 - x92*x93 - x93*x94 + x95*x96 + x96*x98) - x202*(-B1h*x57 - B1p*x12*x193 + B1pp*x132 - B1pp*x158 + B2p*(S*x55 + S*x58 + x112*(Fxt*x84 + St*x18 + x105*x41 + x119*(Sp + x201) + x121*x201 + x122 + x15*(x201*xi_x + x40 + x43 + x88) + x18*x37 + x18*x38 + x84*xi_xx) + x21*x81 - x36*x72 + x52*(x139 + x14 - x190 - x192 - x47 - x48) + x60 - x66*(Fx*(Gh + x83) + x111 + x173 + x193) + x68*(Sp - x36) + x73 - x75 + x76) + B2pp*x132 + B2pp*x158 + Fx*x196 + Fxpt*x170 + Fxt*x171 + Fxt*x183 - Fy*x13*x34 + Fyp*x64 + Fyph*x84 - Gh*x102*x34 - Gh*x89 - Sh*x163 - Sp*x34*x68 + Sp*x55 + Sp*x58 + Sph*x52 + Sph*x59 - Spp*x68 - Spp*x72 + St*x169 - x102*x177 - x103*x80 + x105*x181 + x105*x194 + x105*x195 - x106*x80 - x107*x177 - x109*x74 - x11*x19*x27 - x111*x172 - x119*x189 - x12*x125*x157 + x120*x160 + x120*x161 + x120*x180 + x120*x199 + x121*x128*x170 - x121*x189 - x125*x176 + x125*x197 + x128*x132 + x128*x158 + 8*x131*x198*xi_x - x133*x176 + x133*x197 + x134*x160 + x134*x161 + x134*x199 + x136*x84*(x112*x200 + x149) - x142*x165 - x143*x165 - x144*x165 - x145*x165 - x146*x165 - x15*x164 + x15*x168 + x15*x184 + x15*x186 + x15*x188 - x15*x189*xi_x - x150*x178 - x151*x178 - x153*x22 - x153*x23 - x154*x52 - x154*x59 + x155*x52 + x155*x59 + x156*x52 + x156*x59 - x157*x22 + x159*x52 + x159*x59 - x160*x68 - x160*x72 + x161*x68 + x161*x72 - x162*x42 - x162*x43 - x163*x27 - x164*x18 - x166*x52 - x166*x59 - x167*x21 - x167*x22 - x167*x23 + x168*x18 + x169*x37 + x169*x38 - x169*x42 - x169*x43 + x171*xi_xx - x172*x173 - x174*x53 - x174*x79 - x175*x53 - x175*x79 - x179*x53 - x179*x79 + x18*x184 + x18*x186 + x18*x188 + x180*x68 + x180*x72 + x181*x182 + x182*x194 + x182*x195 + x183*xi_xx + x185*x42 + x185*x43 + x187*x42 + x187*x43 - x190*x191 - x191*x192 + x196*xi_x - x21*x34*xi_y - x22*x70 + x22*x86 - x23*x70 + x23*x86 + x27*x70 - x34*x73 - x34*x77 - x95*xi_y)

    S = S*x224*x229*x6 + Sp*x210*(3*B2d - 3*x209/2) + Sp*x216*x231*(12*B1d + 6*B2d - x230*(B2p + x34)) + 3*x1*x2*x213*x214 - x10*(Fydh + Fydp*x31 - Fyp*x233 - x232*(Fyph + Fypp*x31)) - x12*x236*(-3*x237 + x242) - x135*x229*x256 + x149*x152*x207*x208 - x152*x221*x222 + 3*x202*x204*x205 - x202*x218*x221 + x204*x271 + x208*x215*x216 - x210*x211 - x211*x226 - x213*x256 + x218*x271*x31 + x222*x224*x225*x4 + x234*x30*(Fxdh + Fxdp*x31 - Fxp*x233 - x232*(Fxph + Fxpp*x31)) + x244*x245 - x272*(S*(x11*(x2*x292 + x2*(B1dp*x135 + B1dt - B1p*x279 - x232*(B1pp*x135 + B1pt)) + x215*x265 + x215*(B1p*x135 + B1t) + x267*x284 + x293*x5) - x2*x291 - x207*x264 + x222*x231*x251 + x236*x253/2 - x290*x5) + x241*(Sdt - Sp*x279 - x232*(Spt + x289) + x288) + x244*x254 - x246*x284 + x248*x286 + x249*x285 - x287*x5) + x273*(S*(-x2*x290 + x2*(B1dh + B1dp*x31 - B1p*x233 - x232*(B1ph + B1pp*x31)) + x215*x261 - x215*x262 + x241*x293 + x259*x292 - x263*x284 + x265*x294 + x266*x285 + x267*x286 + x268*x285 - x291*x5) - x2*x287 - x207*x247 + x244*x269 + x258*x294 + x259*(-A*Spt + 2*A*x289 + Ap*x257 - At*Sp - 6*Fxd*Sp + 3*Fxp*x243 + 2*Sdt - 4*x288)/2 + x260*x285) + κ*(S*x245 + x208*x210 + x208*x226 - x255*x272 + x270*x273) + x11*(S*(Gh*x281*x282 + x2*(Ah*x283 + S*(-Ap*x220 + As + 2*Fydh) + x203*x283))*exp(-B1 + B2) - Sdp*x230*x274 + x152*x218*x277 - x2*x229*x240*x278 + x202*x222*x277 - x222*x275*x6*(Fxp*x135 + x239) + x234*(Gt*x280*x281 + St*x2*x280 + x225*(Ab - Ap*Fxt - Ap*xi_xx + At*x252 + B1t*x212 + B2t*x212 + 2*Fxdt)) - x274*(-Ap*x276 + S*(B1d^2*x2^2 + 3*B2d^2 + Gd^2 + 4*phid^2)) + x275*x276*(-x243 + x49) + x275*x30*x8*(Fydp*x135 + Fydt - Fyp*x279 - x232*(Fypp*x135 + Fypt)) - x278*x3*(Fxdp*x135 + Fxdt - Fxp*x279 - x232*(Fxpp*x135 + Fxpt)) - x9*(S*(x2*(Ah*Gt + At*Gh + Fyd*x108 + Gh*x212) + x5*(2*Ac + Ah*B2t - Ap*Fxh - Ap*Fyt - Ap*x147 + At*B2h + B2h*x212 + B2t*x203 + 2*Fxdh + 2*Fydt)) + x5*(At*Sh + Sh*x212 + St*x282)) + (Fxh - Fyt)^2*exp(2*B2))./S


    return axx, ayy, axy, bx, by, cc, S
end

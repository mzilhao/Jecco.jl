
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

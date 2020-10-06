
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

    VV = -3 -3/2 * phi^2 + phi^4*UU(phi,potential)

   # VVp = -3 * phi + phi^4 * UUp(phi,potential) + 4*phi^3 * UU(phi,potential) 

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

   # VV = -3 -3/2 * phi^2 + phi^4*UU(phi,potential)

    VVp = -3 * phi + phi^4 * UUp(phi,potential) + 4*phi^3 * UU(phi,potential)



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

    VV = -3 -3/2 * phi^2 + phi^4*UU(phi,potential)

   # VVp = -3 * phi + phi^4 * UUp(phi,potential) + 4*phi^3 * UU(phi,potential)


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
        B1   , B2   , G   , phi  , S    , Fx    , Fy    , Sd ,  B1d  , B2d  , Gd,  phid, A   ,
        B1p  , B2p  , Gp  , phip , Sp   , Fxp   , Fyp   , Sdp,  B1dp , B2dp , Gdp,       Ap  ,
        B1pp , B2pp , Gpp ,        Spp  , Fxpp  , Fypp  ,                                App ,
        B1_x , B2_x , G_x , phi_x, S_x  , Fx_x  , Fy_x  , Sd_x, B1d_x, B2d_x, Gd_x,      A_x ,
	B1_y , B2_y , G_y , phi_y, S_y  , Fx_y  , Fy_y  , Sd_y, B1d_y, B2d_y, Gd_y,      A_y ,
        B1p_x, B2p_x, Gp_x,        Sp_x , Fxp_x , Fyp_x ,                                Ap_x,
        B1p_y, B2p_y, Gp_y,        Sp_y , Fxp_y , Fyp_y ,                                Ap_y,
                                                  Fy_xx ,                                A_xx,
                                          Fx_yy ,                                        A_yy,
                                          Fx_xy , Fy_xy ,                                A_xy,
    ) = vars

    @tilde_outer("B1")
    @tilde_outer("B2")
    @tilde_outer("G")
    @tilde_outer("phi")
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
    @hat_outer("phi")
    @hat_outer("S")
    @hat_outer("Fx")
    @hat_outer("Fy")
    @hat_outer("Sd")
    @hat_outer("B1d")
    @hat_outer("B2d")
    @hat_outer("Gd")
    @hat_outer("A")

    @bar_outer("A")
    @bar_outer("Fy")

    @star_outer("A")
    @star_outer("Fx")

    @tilde_outer("B1p")
    @tilde_outer("B2p")
    @tilde_outer("Gp")
    @tilde_outer("Sp")
    @tilde_outer("Fxp")
    @tilde_outer("Fyp")
    @tilde_outer("Ap")

    @hat_outer("B1p")
    @hat_outer("B2p")
    @hat_outer("Gp")
    @hat_outer("Sp")
    @hat_outer("Fxp")
    @hat_outer("Fyp")
    @hat_outer("Ap")

    @cross_outer("A")
    @cross_outer("Fx")
    @cross_outer("Fy")

    x0 = exp(B1 + B2)
    x1 = cosh(G)
    x2 = S*x1
    x3 = sinh(G)
    x4 = exp(B2)
    x5 = x3*x4
    x6 = 2*x5
    x7 = -B1
    x8 = exp(B2 + x7)
    x9 = B2h*S
    x10 = Fyp*S
    x11 = Sp*xi_y
    x12 = 2*B2p
    x13 = S*x12
    x14 = x13*xi_y
    x15 = exp(B1)
    x16 = S*x15
    x17 = Gt*x16
    x18 = 2*Fx
    x19 = Gp*x16
    x20 = x18*x19
    x21 = 2*xi_x
    x22 = x19*x21
    x23 = Fy + xi_y
    x24 = Gp*x23
    x25 = S*(Gh + 2*x24)
    x26 = B1t*S
    x27 = 2*B1p
    x28 = x12 + x27
    x29 = Fx*S
    x30 = S*x27
    x31 = B2t*S
    x32 = Fxp*S
    x33 = -Fx*Sp - Sp*xi_x + St + x13*xi_x + x31 + x32
    x34 = -Sh
    x35 = B1h*S
    x36 = 6*Sp
    x37 = Fyp*x23
    x38 = x1*x4
    x39 = x37*x38
    x40 = x23^2
    x41 = x38*x40
    x42 = 3*Spp
    x43 = 3*B2p
    x44 = Sp*x41
    x45 = 3*Sp
    x46 = x40*x45*x5
    x47 = Fx + xi_x
    x48 = Fxp*x47
    x49 = exp(2*B1 + B2)
    x50 = x1*x49
    x51 = x48*x50
    x52 = x47^2
    x53 = x50*x52
    x54 = Gp*x3
    x55 = x45*x49*x52
    x56 = Sp*x53
    x57 = B1p + B2p
    x58 = x15*(x1*(B1t + B2t + x47*x57) + x3*(Fx*Gp + Gp*xi_x + Gt))
    x59 = B2p*x23
    x60 = -B2h - x59
    x61 = x1*(-Gh - x24) + x3*x60
    x62 = x58 + x61
    x63 = Sp*x23
    x64 = Sh + x63
    x65 = x3*x64
    x66 = Sp*x47
    x67 = St + x66
    x68 = x1*x15
    x69 = x67*x68
    x70 = -x65 + x69
    x71 = S*x62 + x70
    x72 = x0*x71
    x73 = 2*xi_xy
    x74 = Fxp*x23
    x75 = x0*x3
    x76 = 6*S
    x77 = Fyp*x47
    x78 = x0*x1
    x79 = x15*(Sd*x76 + x5*(2*Fyt + x73 + 2*x77) - x78*(2*Fxt + 2*x48 + 2*xi_xx))
    x80 = -x38*(2*Fyh + 2*x37 + 2*xi_yy) + x75*(2*Fxh + x73 + 2*x74) + x79
    x81 = St - 2*x66
    x82 = x15*x3
    x83 = x81*x82
    x84 = B1h + B1p*x23
    x85 = B2h + x59
    x86 = x1*x85
    x87 = Gh + x24
    x88 = x3*x87
    x89 = B2p*x47 + B2t
    x90 = x82*x89
    x91 = Gp*x47 + Gt
    x92 = x68*x91
    x93 = x1*x84 - x86 - x88 + x90 + x92
    x94 = S*x93 + x1*(x34 - x63) + x83
    x95 = x4*x94
    x96 = Fyp*x95
    x97 = x23*x95
    x98 = Fyh + xi_yy
    x99 = x37 + x98
    x100 = x38*x99
    x101 = 2*Gp
    x102 = x5*x99
    x103 = Fypp*x23
    x104 = Fxpp*x23
    x105 = Fxh + xi_xy
    x106 = x105 + x74
    x107 = x106*x78
    x108 = x106*x75
    x109 = Fypp*x47 + Fypt
    x110 = Fxpp*x47 + Fxpt
    x111 = 2*x78
    x112 = Fyt + xi_xy
    x113 = x112 + x77
    x114 = x113*x5
    x115 = x113*x38
    x116 = Fxt + xi_xx
    x117 = x116 + x48
    x118 = x117*x75
    x119 = Sph + Spp*x23
    x120 = Gp*x1
    x121 = Spp*x47
    x122 = Spt + x121
    x123 = x67*x82
    x124 = B2ph + B2pp*x23
    x125 = Gph + Gpp*x23
    x126 = B1pp*x47 + B1pt
    x127 = B2pp*x47 + B2pt
    x128 = Gpp*x47 + Gpt
    x129 = B1p*x47 + B1t
    x130 = x0*(x18 + x21)
    x131 = x68*x81
    x132 = Fxp*Sp
    x133 = B1ph + B1pp*x23
    x134 = x1*x89
    x135 = Gp*x15
    x136 = x3*x91
    x137 = x4*(2*Fy + 2*xi_y)
    x138 = S*x80 + x137*x94 + x41*x45 + x45*x53
    x139 = -x130*x71 + x138
    x140 = exp(x7)
    x141 = x140/2
    x142 = A*B1p
    x143 = B1d - x142/2
    x144 = x15*(x1*x129 + x134 + x136)
    x145 = x144 + x61
    x146 = S*x145 + x70
    x147 = A*Sp
    x148 = 3*x147
    x149 = A*Gp
    x150 = Gd - x149/2
    x151 = 3*B2d
    x152 = A*B2p
    x153 = x151 - 3*x152/2
    x154 = A/2
    x155 = Ap*Sp
    x156 = 3*Sdp - 3*Spp*x154 - 3*x155/2
    x157 = x150*x3
    x158 = Sd - x147/2
    x159 = x0*x146
    x160 = 2*B2d - x152
    x161 = 2*B1d
    x162 = -x142 + x160 + x161
    x163 = x1*x150
    x164 = Ah + Ap*x23
    x165 = x164/2
    x166 = Sdh + Sdp*x23 - Sp*x165 - x119*x154
    x167 = Ap*x47 + At
    x168 = x167/2
    x169 = Sdp*x47 + Sdt - Sp*x168 - x122*x154
    x170 = B2dh + B2dp*x23 - B2p*x165 - x124*x154
    x171 = Gdh + Gdp*x23 - Gp*x165 - x125*x154
    x172 = B2dp*x47 + B2dt - B2p*x168 - x127*x154
    x173 = Gdp*x47 + Gdt - Gp*x168 - x128*x154
    x174 = 3*A/2
    x175 = x15*x150
    x176 = 2*Gd
    x177 = -x149 + x176
    x178 = 16*S*Sd
    x179 = 3*Sd
    x180 = x1^2
    x181 = 3*B1d*x180
    x182 = S^2
    x183 = 4*x182
    x184 = 4*S
    x185 = 4*Sh
    x186 = 8*x9
    x187 = 4*x32
    x188 = 4*x10
    x189 = Fxs*x184 - Fyc*x184 - x105*x185 + x105*x186 + x112*x185 - x112*x186 + x112*x188 - x187*x98
    x190 = 4*St
    x191 = 8*x31
    x192 = Fxc*x184 - Fyb*x184 - x105*x187 - x105*x190 + x105*x191 + x112*x190 - x112*x191 + x116*x188
    x193 = 6*B2d
    x194 = 8*phid
    x195 = 2*B1dh
    x196 = Ah*B1p
    x197 = x161*x180
    x198 = 2*G
    x199 = sinh(x198)
    x200 = B1d*x199
    x201 = 2*x200
    x202 = S^3
    x203 = 2*x202
    x204 = 4*Gd
    x205 = 2*B1dt
    x206 = S*x199
    x207 = At*B1p
    x208 = B1d*cosh(x198)
    x209 = Sh*x178/2 - x15*x182*(-At*S*x101 - 12*Gd*St - Gdt*x184 + Gt*x184*x208 + 6*St*x200 + x201*x26 - x204*x26 + x205*x206 + x206*x207)/2 + x183*(Ah*Sp - Fyp*x179 - 4*Sdh - Sh*x151 + Sh*x181)/2 - x189*x5/2 + x192*x78/2 - x203*(Ah*B2p - Ap*Fyp + 2*Aph + B1h*x197 + 2*B2dh + B2h*x193 + Gh*x176 - Gh*x201 + phih*x194 - x180*x195 - x180*x196)/2
    x210 = 1/x202
    x211 = x210*x23
    x212 = x141*(6*Sh*x182*(x176 + x200) - x15*(-St*x178 + x183*(-At*Sp + Fxp*x179 + 4*Sdt + St*x151 + St*x181) + x203*(-Ap*Fxp + 2*Apt + At*B2p + B1t*x197 + 2*B2dt + B2t*x193 + Gt*x176 + Gt*x201 + phit*x194 + x180*x205 + x180*x207)) - x189*x38 + x192*x75 - x202*(-Ah*x101 + B1h*x201 + B1h*x204 - 4*Gdh - 4*Gh*x208 - x195*x199 - x196*x199))
    x213 = 3*S
    x214 = Ah*x1
    x215 = Ah*x3
    x216 = x1*x182
    x217 = At*x1
    x218 = x182*x217
    x219 = At*x3
    x220 = Ap*x3
    x221 = x210*x47

    axx = -x0*x2

    axy = S*x6

    ayy = -x2*x8

    bx = x4*(x1*(-x15*(x26 + x28*x29 + x30*xi_x + x33) + x25) + x3*(-Fy*Sp + Fy*x13 + Sh + x10 - x11 + x14 - x17 - x20 - x22 + x9))

    by = x8*(x1*(Fy*(Sp - x13 + x30) - x10 + x11 - x14 + x17 + x20 + x22 + x30*xi_y + x34 + x35 - x9) + x3*(x15*(x12*x29 + x33) - x25))

    cc = x141*(-B1p*x139 - 2*Fxp*x72 + Gp*x46 + S*(B1p*x79 - x100*x12 - x101*x102 + x101*x107 + x108*x28 + x15*(Sd*x36 + Sdp*x76 + x101*x115 - x101*x118 + x109*x6 - x110*x111 - x111*x117*x57 + x114*x12) - x38*(2*Fyph + 2*x103) + x75*(2*Fxph + 2*x104)) + Sp*x80 + x12*x97 - x130*(B1p*x69 + Gp*x123 + S*(B1p*x58 - Gp*x86 - Gp*x88 - x1*x125 - x124*x3 + x15*(x1*x126 + x1*x127 + x120*x91 + x128*x3 + x129*x54 + x54*x89)) + Sp*x62 - x119*x3 - x120*x64 + x122*x68) + x137*(B1p*x83 + Gp*x131 - Gp*x65 + S*(B1p*x90 + B1p*x92 - x1*x124 + x1*x133 - x120*x87 - x125*x3 + x127*x82 + x128*x68 + x134*x135 + x135*x136 + x54*x84 - x54*x85) + Sp*x93 - x1*x119 + x82*(Spt - 2*x121 - 3*x132)) - x28*x47*x72 + x36*x39 + x36*x51 + x41*x42 + x42*x53 + x43*x44 + x54*x55 + x56*(6*B1p + x43) + 2*x96)

    SS = kappa*x139*x141 - x141*x143*(-x130*x146 + x138) + x141*(A*Fxp*x159 - A*x96 + S*(-x100*x160 - x102*x177 + x107*x177 + x108*x162 + x143*x79 + x15*(2*Sd*(-3*x147/2 + x179) + x114*x160 + x115*x177 - x117*x162*x78 - x118*x177 + x213*(-A*Sdp + x140*x210*(-S*x4*(Ap*x2*x98 - As*x2 - Gh*S*x215 - Sh*x214 + x214*x35 - x214*x9) - x0*(S*x3*(Ah*St + At*Sh) + x182*(2*Ac*x3 + B2h*x219 + B2t*x215 + Gh*x217 + Gt*x214 - x105*x220 - x112*x220)) - x15*x202*(-Ap*x179 + B1d^2*S*x180 + B2d^2*x213 + Gd^2*S + phid^2*x184) - x49*(-Ab*x216 + Ap*x116*x216 - At*St*x2 - B1t*x218 - B2t*x218 - Gt*x182*x219) + (Fxh - Fyt)^2*exp(B1 + 2*B2))/3) + x5*(-A*x109 - Fyp*x167 + x209*x221) - x78*(-A*x110 - Fxp*x167 + x212*x221)) - x38*(-A*(Fyph + x103) - Fyp*x164 + x209*x211) + x75*(-A*(Fxph + x104) - Fxp*x164 + x211*x212)) - x130*(S*(-x1*x171 + x143*x144 + x15*(x1*x172 + x1*(B1dp*x47 + B1dt - B1p*x168 - x126*x154) + x129*x157 + x157*x89 + x163*x91 + x173*x3) - x150*x88 + x163*x60 - x170*x3) + x123*x150 + x143*x69 + x145*x158 - x163*x64 - x166*x3 + x169*x68) + x137*(S*(-x1*x170 + x1*(B1dh + B1dp*x23 - B1p*x165 - x133*x154) + x134*x175 + x136*x175 + x143*x90 + x143*x92 + x157*x84 - x157*x85 - x163*x87 - x171*x3 + x172*x82 + x173*x68) - x1*x166 + x131*x150 + x143*x83 - x150*x65 + x158*x93 + x82*(x132*x174 + x169 - x47*(3*Sdp - Spp*x174 - 3*x155/2))) - x148*x39 - x148*x51 + x150*x46 + x153*x44 + x156*x41 + x156*x53 + x157*x55 + x158*x80 - x159*x162*x47 + x160*x97 + x56*(6*B1d - 3*x142 + x153))

    return axx, ayy, axy, bx, by, cc, SS
end

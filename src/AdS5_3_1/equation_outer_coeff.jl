
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

    κ = kappa


    x0 = 2*B1
    x1 = exp(B2 + x0)
    x2 = cosh(G)
    x3 = 2*x2
    x4 = S*x3
    x5 = exp(B1 + B2)
    x6 = sinh(G)
    x7 = x5*x6
    x8 = 4*S
    x9 = exp(B2)
    x10 = B2h*S
    x11 = Fyp*S
    x12 = Sp*xi_y
    x13 = 2*B2p
    x14 = S*x13
    x15 = x14*xi_y
    x16 = exp(B1)
    x17 = Gt*x16
    x18 = S*x17
    x19 = 2*Fx
    x20 = Gp*S*x16
    x21 = x19*x20
    x22 = 2*xi_x
    x23 = x20*x22
    x24 = 2*x6
    x25 = Fy + xi_y
    x26 = Gp*x25
    x27 = S*(Gh + 2*x26)
    x28 = B1t*S
    x29 = 2*B1p
    x30 = x13 + x29
    x31 = Fx*S
    x32 = S*x29
    x33 = B2t*S
    x34 = Sp*xi_x
    x35 = Fxp*S
    x36 = St + x35
    x37 = -Fx*Sp + x14*xi_x + x33 - x34 + x36
    x38 = -Sh
    x39 = -x10 + x18 + x38
    x40 = 6*Sp
    x41 = Fyp*x25
    x42 = x2*x9
    x43 = x41*x42
    x44 = 3*Spp
    x45 = x25^2
    x46 = x42*x45
    x47 = 3*B2p
    x48 = Sp*x46
    x49 = 3*Sp
    x50 = x6*x9
    x51 = x45*x50
    x52 = Fx + xi_x
    x53 = Fxp*x52
    x54 = x1*x2
    x55 = x53*x54
    x56 = x52^2
    x57 = x54*x56
    x58 = Gp*x6
    x59 = x1*x56
    x60 = Sp*x57
    x61 = Sp*x25
    x62 = Sh + x61
    x63 = x6*x62
    x64 = Sp*x52
    x65 = St + x64
    x66 = x16*x2
    x67 = x65*x66
    x68 = -B2h
    x69 = B2p*x25
    x70 = x68 - x69
    x71 = B1p + B2p
    x72 = x16*(x2*(B1t + B2t + x52*x71) + x6*(Fx*Gp + Gp*xi_x + Gt))
    x73 = x2*(-Gh - x26) + x6*x70 + x72
    x74 = S*x73 - x63 + x67
    x75 = x5*x74
    x76 = Fxp*x75
    x77 = x42*(2*Fyh + 2*x41 + 2*xi_yy)
    x78 = 2*Fxh
    x79 = 2*xi_xy
    x80 = Fxp*x25
    x81 = 2*x80
    x82 = 6*S
    x83 = Sd*x82
    x84 = Fyp*x52
    x85 = x50*(2*Fyt + x79 + 2*x84)
    x86 = x2*x5
    x87 = x86*(2*Fxt + 2*x53 + 2*xi_xx)
    x88 = x16*(x83 + x85 - x87)
    x89 = x7*(x78 + x79 + x81) - x77 + x88
    x90 = St - 2*x64
    x91 = x16*x6
    x92 = x90*x91
    x93 = B1h + B1p*x25
    x94 = B2h + x69
    x95 = x2*x94
    x96 = Gh + x26
    x97 = x6*x96
    x98 = B2p*x52 + B2t
    x99 = x91*x98
    x100 = Gp*x52 + Gt
    x101 = x100*x66
    x102 = x101 + x2*x93 - x95 - x97 + x99
    x103 = S*x102 + x2*(x38 - x61) + x92
    x104 = x103*x9
    x105 = Fyp*x104
    x106 = x104*x25
    x107 = x52*x75
    x108 = Fyh + xi_yy
    x109 = x108 + x41
    x110 = x24*x9
    x111 = x109*x110
    x112 = Fypp*x25
    x113 = Fxpp*x25
    x114 = Fxh + x80 + xi_xy
    x115 = x3*x5
    x116 = x114*x7
    x117 = 6*Sd
    x118 = Fypp*x52 + Fypt
    x119 = Fxpp*x52 + Fxpt
    x120 = Fyt + x84 + xi_xy
    x121 = x120*x50
    x122 = Fxt + x53 + xi_xx
    x123 = Sph + Spp*x25
    x124 = Gp*x2
    x125 = Spp*x52
    x126 = Spt + x125
    x127 = x16*x58
    x128 = B2ph + B2pp*x25
    x129 = Gph + Gpp*x25
    x130 = B1pp*x52 + B1pt
    x131 = B2pp*x52 + B2pt
    x132 = Gpp*x52 + Gpt
    x133 = B1p*x52 + B1t
    x134 = x5*(x19 + x22)
    x135 = Gp*x66
    x136 = 3*Fxp
    x137 = B1ph + B1pp*x25
    x138 = x9*(2*Fy + 2*xi_y)
    x139 = A*Sp
    x140 = 3*x139
    x141 = 6*B2d
    x142 = 2*Gd
    x143 = A*Gp
    x144 = x142 - x143
    x145 = 3*x144/2
    x146 = Gd - x143/2
    x147 = 2*Sdp
    x148 = A*Spp
    x149 = Ap*Sp
    x150 = 3*x147/2 - 3*x148/2 - 3*x149/2
    x151 = S*x42
    x152 = 4*B2d
    x153 = -A*x13 + x152
    x154 = x153/2
    x155 = Sp*x59
    x156 = 4*B1d
    x157 = -2*A*x71 + x152 + x156
    x158 = x157/2
    x159 = 3*A
    x160 = x2/2
    x161 = 2*B1d
    x162 = A*B1p
    x163 = x161/2 - x162/2
    x164 = 2*Sd - x139
    x165 = S^3
    x166 = 2*x165
    x167 = A*x166
    x168 = Ah + Ap*x25
    x169 = x166*x168
    x170 = S^2
    x171 = 3*B2d
    x172 = 3*Sd
    x173 = x2^2
    x174 = 3*x173
    x175 = 4*Fxh
    x176 = 4*Fyt
    x177 = x175*(Sh - 2*x10) - x176*(S*(-2*B2h + Fyp) + Sh) + x8*(Fxp*x108 - Fxs + Fyc - Fyp*xi_xy)
    x178 = 4*x11
    x179 = Fxh - Fyt
    x180 = Fxc*x8 + Fxt*x178 - Fyb*x8 + St*x176 - x175*x36 + x178*xi_xx + 8*x179*x33 - 4*x35*xi_xy
    x181 = 8*phid
    x182 = 2*G
    x183 = sinh(x182)
    x184 = x161*x183
    x185 = Ah*B1p + 2*B1dh - B1h*x161
    x186 = 2*Gp
    x187 = 4*Gd
    x188 = x156*cosh(x182)
    x189 = B1d*St
    x190 = 2*B1dt
    x191 = At*B1p
    x192 = 16*S*Sd*Sh + x16*x170*(At*S*x186 + 12*Gd*St + Gdt*x8 - Gt*S*x188 - x183*(S*x190 + S*x191 + x161*x28 + 6*x189) + x187*x28) - x166*(Ah*B2p - Ap*Fyp + 2*Aph + 2*B2dh + B2h*x141 + Gh*x142 - Gh*x184 + phih*x181 - x173*x185) + 4*x170*(Ah*Sp + B1d*Sh*x174 - Fyp*x172 - 4*Sdh - Sh*x171) + x177*x50 + x180*x86
    x193 = 1/(2*x170)
    x194 = -At*Sp
    x195 = 2*S
    x196 = x16*x195
    x197 = exp(-B1)
    x198 = x197*(-12*Sh*x170*(B1d*x2*x6 + Gd) + x165*(-Ah*x186 + B1h*x187 - 4*Gdh - Gh*x188 - x183*x185) - x177*x42 - x180*x7 + x196*(-8*Sd*St + x170*(-Ap*Fxp + 2*Apt + At*B2p + 2*B2dt + B2t*x141 + Gt*x142 + Gt*x184 + phit*x181 + x173*(B1t*x161 + x190 + x191)) + x195*(Sd*x136 + 4*Sdt + St*x171 + x174*x189 + x194)))
    x199 = x146*x2
    x200 = B1d - x162/2
    x201 = x146*x91
    x202 = A/2
    x203 = x168/2
    x204 = Sdh + Sdp*x25 - Sp*x203 - x123*x202
    x205 = Ap*x52 + At
    x206 = x205/2
    x207 = Sd - x139/2
    x208 = B2dh + B2dp*x25 - B2p*x203 - x128*x202
    x209 = Gdh + Gdp*x25 - Gp*x203 - x129*x202
    x210 = x146*x6
    x211 = B2dp*x52 + B2dt - B2p*x206 - x131*x202
    x212 = Gdp*x52 + Gdt - Gp*x206 - x132*x202
    x213 = x146*x66
    x214 = 4*x146
    x215 = x122*x165
    x216 = Ah*S
    x217 = Ah*x16
    x218 = Ap*S
    x219 = x16*x218
    x220 = At*x16
    x221 = exp(x0)
    x222 = At*x221
    x223 = x218*x221
    x224 = x166*x205

    axx = -x1*x4

    axy = x7*x8

    ayy = -x4*x9

    bx = x5*(x24*(-Fy*Sp + Fy*x14 + Sh + x10 + x11 - x12 + x15 - x18 - x21 - x23) + x3*(-x16*(x28 + x30*x31 + x32*xi_x + x37) + x27))

    by = x9*(x24*(x16*(x13*x31 + x37) - x27) + x3*(B1h*S + Fy*(Sp - x14 + x32) - x11 + x12 - x15 + x21 + x23 + x32*xi_y + x39))

    cc = Gp*x49*x51 + S*(B1p*x88 - Gp*x111 + Gp*x114*x115 - x109*x13*x42 + x116*x30 + x16*(Gp*x120*x3*x9 - Gp*x122*x24*x5 + Sdp*x82 + Sp*x117 + x110*x118 - x115*x119 - x115*x122*x71 + x121*x13) - x42*(2*Fyph + 2*x112) + x7*(2*Fxph + 2*x113)) + Sp*x89 + 2*x105 + x106*x13 - x107*x30 - x134*(B1p*x67 + S*(B1p*x72 - Gp*x95 - Gp*x97 - x128*x6 - x129*x2 + x16*(x100*x124 + x130*x2 + x131*x2 + x132*x6 + x133*x58 + x58*x98)) + Sp*x73 - x123*x6 - x124*x62 + x126*x66 + x127*x65) + x138*(B1p*x92 - Gp*x63 + S*(B1p*x101 + B1p*x99 + x100*x127 - x124*x96 - x128*x2 - x129*x6 + x131*x91 + x132*x66 + x135*x98 + x137*x2 + x58*x93 - x58*x94) + Sp*x102 - x123*x2 + x135*x90 + x91*(-Sp*x136 + Spt - 2*x125)) + x40*x43 + x40*x55 + x44*x46 + x44*x57 + x47*x48 + x49*x58*x59 + x60*(6*B1p + x47) - 2*x76

    S = -A*x105 + A*x76 - S*x111*x146 + S*x116*x158 + S*x163*x88 + Sp*x145*x51 + x106*x154 - x107*x158 - x109*x151*x154 + x114*x146*x4*x5 - x134*(S*(x144*x160*x70 - x146*x97 + x16*(x100*x199 + x133*x210 + x2*x211 + x2*(B1dp*x52 + B1dt - B1p*x206 - x130*x202) + x210*x98 + x212*x6) + x163*x72 - x2*x209 - x208*x6) - x199*x62 + x200*x67 + x201*x65 - x204*x6 + x207*x73 + x66*(Sdp*x52 + Sdt - Sp*x206 - x126*x202)) + x138*(S*(x100*x201 + x101*x200 - x199*x96 - x2*x208 + x2*(B1dh + B1dp*x25 - B1p*x203 - x137*x202) + x200*x99 - x209*x6 + x210*x93 - x210*x94 + x211*x91 + x212*x66 + x213*x98) + x102*x207 - x146*x63 - x2*x204 + x200*x92 + x213*x90 + x91*(-A*Spt + 2*Ap*x34 - 4*Sdp*xi_x + 2*Sdt + x136*x139 + x148*x22 + x19*(-x147 + x148 + x149) + x194)/2) - x140*x43 - x140*x55 + x145*x155*x6 + x150*x46 + x150*x57 + x155*x160*(12*B1d + x141 - x159*(B2p + x29)) + x16*x193*(x117*x164*x165 + x120*x165*x214*x42 + x121*x153*x165 - x157*x215*x86 - x195*x197*(B1d^2*S^4*x16*x173 - S*x50*(-Ac*x196 + Fxh*x219 + Fyt*x219 + Gh*x216 - St*x217 - x217*x33 + x219*x79 + x220*x39) + x151*(-Ab*S*x221 - Ah*Sh - As*S + Fxt*x223 + Fyh*x218 + Gh*S*x220 - St*x222 + x216*(B1h + x17 + x68) + x218*xi_yy - x222*x28 - x222*x33 + x223*xi_xx) + x16*(x165*(-Ap*x172 + S*(3*B2d^2 + Gd^2 + 4*phid^2) + Sdp*x159) - x179^2*exp(2*B2))) - x214*x215*x7 - x50*(Fyp*x224 + x118*x167 - x192*x52) + x86*(Fxp*x224 + x119*x167 + x198*x52)) - x164*(x16*(-x83 - x85 + x87) + x7*(-x78 - x79 - x81) + x77)/2 + x193*x42*(Fyp*x169 + x167*(Fyph + x112) - x192*x25) - x193*x7*(Fxp*x169 + x167*(Fxph + x113) + x198*x25) + x48*(-A*x47 + x141)/2 + κ*(S*x89 + x103*x138 - x134*x74 + 3*x48 + 3*x60)


    return axx, ayy, axy, bx, by, cc, S
end

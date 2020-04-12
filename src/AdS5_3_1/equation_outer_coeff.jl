
# Vf(phi)  = VV(phi)
# Vfp(phi) = ∂(VV)(phi)

# FIXME
V(phi) = -3.0
Vp(phi) = 0.0


# assuming
# (A d_uu + B d_u + C Id) f = -S

function S_outer_eq_coeff!(ABCS::Vector, vars::AllVarsOuter)
    u   = vars.u

    B1p  = vars.B1p
    B2p  = vars.B2p

    G   = vars.G
    Gp  = vars.Gp

    phip = vars.phip

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

function Fxy_outer_eq_coeff!(AA::Matrix, BB::Matrix, CC::Matrix, SS::Vector, vars::FxyVars)
    u     = vars.u

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


function Sd_outer_eq_coeff!(ABCS::Vector, vars::AllVarsOuter)
    u   = vars.u

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

    ABCS[4] = *(*(S, *(*(-8, Gh, St) + *(-8, Gt, Sh), coshG) + *(*(-16, Sc) + *(2, Sh, Fxp + *(-4, B2t)) + *(2, Sp, *(4, Fxh) + *(4, Fyt) + *(8, xi_xy)) + *(2, St, Fyp + *(-4, B2h)), sinhG)) + *(S ^ 2, *(*(-4, Gc) + *(-2, Gh, B1t + B2t + *(-1, Fxp)) + *(2, Gp, Fxh + Fyt + *(2, xi_xy)) + *(2, Gt, B1h + Fyp + *(-1, B2h)), coshG) + *(*(-4, B2c) + *(2, Fxph) + *(2, Fypt) + *(-8, phih, phit) + *(-4, Gh, Gt) + *(2, B2h, Fxp + *(-4, B2t)) + *(2, B2p, Fxh + Fyt + *(2, xi_xy)) + *(2, Fyp, B2t + *(-1, Fxp)), sinhG)) + *(8, Sh, St, sinhG), exp(B1 + B2)) + *(*(S, *(*(8, Ss) + *(-2, Sh, Fyp + *(-4, B2h) + *(4, B1h)) + *(-2, Sp, *(4, Fyh) + *(4, xi_yy)), coshG) + *(8, Gh, Sh, sinhG)) + *(S ^ 2, *(*(2, Gs) + *(-2, Gp, Fyh + xi_yy) + *(2, Gh, B2h + *(-1, Fyp) + *(-2, B1h)), sinhG) + *(Fyp ^ 2 + *(-2, B1s) + *(-2, Fyph) + *(2, B2s) + *(2, B1h ^ 2) + *(2, Gh ^ 2) + *(4, B2h ^ 2) + *(4, phih ^ 2) + *(Fyp, *(-2, B2h) + *(2, B1h)) + *(-2, B1h, B2h) + *(2, B1p + *(-1, B2p), Fyh + xi_yy), coshG)) + *(-4, Sh ^ 2, coshG), expB2) + *(*(S, *(*(8, Sb) + *(-8, Fxt, Sp) + *(-8, Sp, xi_xx) + *(-2, Fxp, St) + *(8, B1t, St) + *(8, B2t, St), coshG) + *(8, Gt, St, sinhG)) + *(S ^ 2, *(*(2, Gb) + *(-2, Gp, Fxt + xi_xx) + *(2, Gt, B2t + *(-1, Fxp) + *(2, B1t)), sinhG) + *(Fxp ^ 2 + *(-2, Fxpt) + *(2, B1b) + *(2, B2b) + *(2, B1t ^ 2) + *(2, Gt ^ 2) + *(4, B2t ^ 2) + *(4, phit ^ 2) + *(-1, Fxp, *(2, B1t) + *(2, B2t)) + *(-2, B1p + B2p, Fxt + xi_xx) + *(2, B1t, B2t), coshG)) + *(-4, St ^ 2, coshG), exp(B2 + *(2, B1))) + *(8, V(phi), S ^ 4, expB1)

    nothing
end


function B2d_outer_eq_coeff!(ABCS::Vector, vars::AllVarsOuter)
    u   = vars.u

    xi_xx  = vars.xi_xx
    xi_xy  = vars.xi_xy
    xi_yy  = vars.xi_yy

    B1     = vars.B1
    B2     = vars.B2
    G      = vars.G
    phi    = vars.phi
    S      = vars.S
    Fx     = vars.Fx
    Fy     = vars.Fy
    Sd     = vars.Sd

    B1p    = vars.B1p
    B2p    = vars.B2p
    Gp     = vars.Gp
    phip   = vars.phip
    Sp     = vars.Sp
    Fxp    = vars.Fxp
    Fyp    = vars.Fyp

    B1t    = vars.B1t
    B2t    = vars.B2t
    Gt     = vars.Gt
    phit   = vars.phit
    St     = vars.St
    Fxt    = vars.Fxt
    Fyt    = vars.Fyt

    B1h    = vars.B1h
    B2h    = vars.B2h
    Gh     = vars.Gh
    phih   = vars.phih
    Sh     = vars.Sh
    Fxh    = vars.Fxh
    Fyh    = vars.Fyh

    B1b   = vars.B1b
    B2b   = vars.B2b
    Gb    = vars.Gb
    phib  = vars.phib
    Sb    = vars.Sb
    # Fxb   = vars.Fxb
    # Fyb   = vars.Fyb

    B1s   = vars.B1s
    B2s   = vars.B2s
    Gs    = vars.Gs
    phis  = vars.phis
    Ss    = vars.Ss
    # Fxs   = vars.Fxs
    # Fys   = vars.Fys

    B1pt   = vars.B1pt
    B2pt   = vars.B2pt
    Gpt    = vars.Gpt
    phipt  = vars.phipt
    Spt    = vars.Spt
    Fxpt   = vars.Fxpt
    Fypt   = vars.Fypt

    B1ph   = vars.B1ph
    B2ph   = vars.B2ph
    Gph    = vars.Gph
    phiph  = vars.phiph
    Sph    = vars.Sph
    Fxph   = vars.Fxph
    Fyph   = vars.Fyph

    B2c   = vars.B2c
    Gc    = vars.Gc
    Sc    = vars.Sc


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

function B1dGd_outer_eq_coeff!(AA::Matrix, BB::Matrix, CC::Matrix, SS::Vector, vars::AllVarsOuter)
    u    = vars.u

    xi_xx  = vars.xi_xx
    xi_xy  = vars.xi_xy
    xi_yy  = vars.xi_yy

    B1     = vars.B1
    B2     = vars.B2
    G      = vars.G
    phi    = vars.phi
    S      = vars.S
    Fx     = vars.Fx
    Fy     = vars.Fy
    Sd     = vars.Sd

    B1p    = vars.B1p
    B2p    = vars.B2p
    Gp     = vars.Gp
    phip   = vars.phip
    Sp     = vars.Sp
    Fxp    = vars.Fxp
    Fyp    = vars.Fyp

    B1t    = vars.B1t
    B2t    = vars.B2t
    Gt     = vars.Gt
    phit   = vars.phit
    St     = vars.St
    Fxt    = vars.Fxt
    Fyt    = vars.Fyt

    B1h    = vars.B1h
    B2h    = vars.B2h
    Gh     = vars.Gh
    phih   = vars.phih
    Sh     = vars.Sh
    Fxh    = vars.Fxh
    Fyh    = vars.Fyh

    B1b   = vars.B1b
    B2b   = vars.B2b
    Gb    = vars.Gb
    phib  = vars.phib
    Sb    = vars.Sb
    # Fxb   = vars.Fxb
    # Fyb   = vars.Fyb

    B1s   = vars.B1s
    B2s   = vars.B2s
    Gs    = vars.Gs
    phis  = vars.phis
    Ss    = vars.Ss
    # Fxs   = vars.Fxs
    # Fys   = vars.Fys

    B1pt   = vars.B1pt
    B2pt   = vars.B2pt
    Gpt    = vars.Gpt
    phipt  = vars.phipt
    Spt    = vars.Spt
    Fxpt   = vars.Fxpt
    Fypt   = vars.Fypt

    B1ph   = vars.B1ph
    B2ph   = vars.B2ph
    Gph    = vars.Gph
    phiph  = vars.phiph
    Sph    = vars.Sph
    Fxph   = vars.Fxph
    Fyph   = vars.Fyph

    B2c   = vars.B2c
    Gc    = vars.Gc
    Sc    = vars.Sc


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


function phid_outer_eq_coeff!(ABCS::Vector, vars::AllVarsOuter)
    u   = vars.u

    xi_xx  = vars.xi_xx
    xi_xy  = vars.xi_xy
    xi_yy  = vars.xi_yy

    B1     = vars.B1
    B2     = vars.B2
    G      = vars.G
    phi    = vars.phi
    S      = vars.S
    Fx     = vars.Fx
    Fy     = vars.Fy
    Sd     = vars.Sd

    B1p    = vars.B1p
    B2p    = vars.B2p
    Gp     = vars.Gp
    phip   = vars.phip
    Sp     = vars.Sp
    Fxp    = vars.Fxp
    Fyp    = vars.Fyp

    B1t    = vars.B1t
    B2t    = vars.B2t
    Gt     = vars.Gt
    phit   = vars.phit
    St     = vars.St
    Fxt    = vars.Fxt
    Fyt    = vars.Fyt

    B1h    = vars.B1h
    B2h    = vars.B2h
    Gh     = vars.Gh
    phih   = vars.phih
    Sh     = vars.Sh
    Fxh    = vars.Fxh
    Fyh    = vars.Fyh

    B1b   = vars.B1b
    B2b   = vars.B2b
    Gb    = vars.Gb
    phib  = vars.phib
    Sb    = vars.Sb
    # Fxb   = vars.Fxb
    # Fyb   = vars.Fyb

    B1s   = vars.B1s
    B2s   = vars.B2s
    Gs    = vars.Gs
    phis  = vars.phis
    Ss    = vars.Ss
    # Fxs   = vars.Fxs
    # Fys   = vars.Fys

    B1pt   = vars.B1pt
    B2pt   = vars.B2pt
    Gpt    = vars.Gpt
    phipt  = vars.phipt
    Spt    = vars.Spt
    Fxpt   = vars.Fxpt
    Fypt   = vars.Fypt

    B1ph   = vars.B1ph
    B2ph   = vars.B2ph
    Gph    = vars.Gph
    phiph  = vars.phiph
    Sph    = vars.Sph
    Fxph   = vars.Fxph
    Fyph   = vars.Fyph

    B2c   = vars.B2c
    Gc    = vars.Gc
    Sc    = vars.Sc
    phic  = vars.phic


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

    ABCS[4] = *(*(*(4, phih, Sh + *(-1, S, B1h + Fyp + *(-1, B2h))) + *(4, S, phis + *(-1, phip, Fyh + xi_yy)), coshG) + *(4, Gh, phih, S, sinhG), expB2) + *(*(4, S, *(phib + *(phit, B1t + B2t + *(-1, Fxp)) + *(-1, phip, Fxt + xi_xx), coshG) + *(Gt, phit, sinhG)) + *(4, phit, St, coshG), exp(B2 + *(2, B1))) + *(*(-1, *(4, phih, St) + *(4, phit, Sh), sinhG) + *(S, *(-4, Gh, phit) + *(-4, Gt, phih), coshG) + *(S, *(-8, phic) + *(4, phih, Fxp + *(-1, B2t)) + *(4, phip, Fxh + Fyt + *(2, xi_xy)) + *(4, phit, Fyp + *(-1, B2h)), sinhG), exp(B1 + B2)) + *(-4, S ^ 2, *(S, Vp(phi)) + *(-3, phip, Sd), expB1)

    nothing
end


function A_outer_eq_coeff!(ABCS::Vector, vars::AllVarsOuter)
    u   = vars.u

    xi_xx  = vars.xi_xx
    xi_xy  = vars.xi_xy
    xi_yy  = vars.xi_yy

    B1     = vars.B1
    B2     = vars.B2
    G      = vars.G
    phi    = vars.phi
    S      = vars.S
    Fx     = vars.Fx
    Fy     = vars.Fy
    Sd     = vars.Sd

    B1d    = vars.B1d
    B2d    = vars.B2d
    Gd     = vars.Gd
    phid   = vars.phid


    B1p    = vars.B1p
    B2p    = vars.B2p
    Gp     = vars.Gp
    phip   = vars.phip
    Sp     = vars.Sp
    Fxp    = vars.Fxp
    Fyp    = vars.Fyp

    B1t    = vars.B1t
    B2t    = vars.B2t
    Gt     = vars.Gt
    phit   = vars.phit
    St     = vars.St
    Fxt    = vars.Fxt
    Fyt    = vars.Fyt

    B1h    = vars.B1h
    B2h    = vars.B2h
    Gh     = vars.Gh
    phih   = vars.phih
    Sh     = vars.Sh
    Fxh    = vars.Fxh
    Fyh    = vars.Fyh

    B1b   = vars.B1b
    B2b   = vars.B2b
    Gb    = vars.Gb
    phib  = vars.phib
    Sb    = vars.Sb
    # Fxb   = vars.Fxb
    # Fyb   = vars.Fyb

    B1s   = vars.B1s
    B2s   = vars.B2s
    Gs    = vars.Gs
    phis  = vars.phis
    Ss    = vars.Ss
    # Fxs   = vars.Fxs
    # Fys   = vars.Fys

    B1pt   = vars.B1pt
    B2pt   = vars.B2pt
    Gpt    = vars.Gpt
    phipt  = vars.phipt
    Spt    = vars.Spt
    Fxpt   = vars.Fxpt
    Fypt   = vars.Fypt

    B1ph   = vars.B1ph
    B2ph   = vars.B2ph
    Gph    = vars.Gph
    phiph  = vars.phiph
    Sph    = vars.Sph
    Fxph   = vars.Fxph
    Fyph   = vars.Fyph

    B2c   = vars.B2c
    Gc    = vars.Gc
    Sc    = vars.Sc
    phic  = vars.phic


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

    ABCS[4] = *(*(2, S ^ 4, *(-4, V(phi)) + *(3, Gd, Gp) + *(9, B2d, B2p) + *(12, phid, phip) + *(3, B1d, B1p, coshGsq)) + *(-72, Sd, Sp, S ^ 2), expB1) + *(*(S, *(*(-24, Ss) + *(24, Sh, B1h + *(-1, B2h)) + *(24, Sp, Fyh + xi_yy), coshG) + *(-24, Gh, Sh, sinhG)) + *(S ^ 2, *(*(-6, Gs) + *(-6, B2h, Gh) + *(6, Fyh, Gp) + *(6, Gp, xi_yy) + *(12, B1h, Gh), sinhG) + *(*(-12, B2h ^ 2) + *(-12, phih ^ 2) + *(-6, B2s) + *(-6, B1h ^ 2) + *(-6, Gh ^ 2) + *(3, Fyp ^ 2) + *(6, B1s) + *(-6, B1p + *(-1, B2p), Fyh + xi_yy) + *(6, B1h, B2h), coshG)) + *(12, Sh ^ 2, coshG), expB2) + *(*(S, *(*(-24, Sb) + *(-24, B1t, St) + *(-24, B2t, St) + *(24, Fxt, Sp) + *(24, Sp, xi_xx), coshG) + *(-24, Gt, St, sinhG)) + *(S ^ 2, *(*(-12, B2t ^ 2) + *(-12, phit ^ 2) + *(-6, B1b) + *(-6, B2b) + *(-6, B1t ^ 2) + *(-6, Gt ^ 2) + *(3, Fxp ^ 2) + *(-6, B1t, B2t) + *(6, B1p + B2p, Fxt + xi_xx), coshG) + *(-1, *(6, Gb) + *(-6, Gp, Fxt + xi_xx) + *(6, Gt, B2t + *(2, B1t)), sinhG)) + *(12, St ^ 2, coshG), exp(B2 + *(2, B1))) + *(*(S ^ 2, *(*(12, Gc) + *(-6, Gp, Fxh + Fyt + *(2, xi_xy)) + *(6, Gh, B1t + B2t) + *(6, Gt, B2h + *(-1, B1h)), coshG) + *(*(12, B2c) + *(-6, B2p, Fxh + Fyt + *(2, xi_xy)) + *(-6, Fxp, Fyp) + *(12, Gh, Gt) + *(24, B2h, B2t) + *(24, phih, phit), sinhG)) + *(24, S, *(*(Gh, St) + *(Gt, Sh), coshG) + *(*(2, Sc) + *(B2h, St) + *(B2t, Sh) + *(-1, Sp, Fxh + Fyt + *(2, xi_xy)), sinhG)) + *(-24, Sh, St, sinhG), exp(B1 + B2))

    nothing
end

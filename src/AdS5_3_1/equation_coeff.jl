
# Vf(phi)  = VV(phi)
# Vfp(phi) = ∂(VV)(phi)

# FIXME
V(phi) = -3.0

# assuming
# (A d_uu + B d_u + C Id) f = -S

function S_outer_eq_coeff!(ABCS::Vector, vars::AllVars)
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


#= Notation

for any function f we're using the following convention: _x denotes partial
derivative with respect to x and

fp  = f_r  = -u^2 f_u
fpp = f_rr = 2u^3 f_u + u^4 f_uu

=#

mutable struct FxyVars{T}
    u        :: T

    B1       :: T
    B1p      :: T
    B1_x     :: T
    B1_y     :: T
    B1pp     :: T
    B1p_x    :: T
    B1p_y    :: T

    B2       :: T
    B2p      :: T
    B2_x     :: T
    B2_y     :: T
    B2pp     :: T
    B2p_x    :: T
    B2p_y    :: T

    G        :: T
    Gp       :: T
    Gpp      :: T
    G_x      :: T
    G_y      :: T
    Gp_x     :: T
    Gp_y     :: T

    phi      :: T
    phip     :: T
    phi_x    :: T
    phi_y    :: T

    S        :: T
    Sp       :: T
    S_x      :: T
    S_y      :: T
    Spp      :: T
    Sp_x     :: T
    Sp_y     :: T
end
function FxyVars{T}() where {T<:AbstractFloat}
    N = 1 + 4*7 + 4
    array = zeros(N)
    FxyVars{T}(array...)
end


# this is a coupled equation for Fx and Fy. the notation used is
#
# ( A11 d_uu Fx + A12 d_uu Fy + B11 d_u Fx + B12 d_u Fy + C11 Fx + C12 Fy ) = -S1
# ( A21 d_uu Fx + A22 d_uu Fy + B21 d_u Fx + B22 d_u Fy + C21 Fx + C22 Fy ) = -S2

function Fxy_outer_eq_coeff!(AA::Matrix, BB::Matrix, CC::Matrix, SS::Vector, vars::FxyVars)
    u    = vars.u

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
    sinh2G  = sinh(*(2, G))
    cosh2G  = cosh(*(2, G))
    coshGsq = cosh(G)^2


    AA[1,1] = *(2, S ^ 2, u ^ 4, expB1)
    AA[1,2] = 0
    AA[2,1] = 0
    AA[2,2] = *(2, S ^ 2, u ^ 4)

    BB[1,1] = *(2, S, u ^ 2, *(-1, Sp) + *(-1, B2p, S) + *(2, S, u) + *(-1, B1p, S, coshGsq), expB1)
    BB[1,2] = *(2, S ^ 2, u ^ 2, Gp + *(1 / 2, B1p, sinh2G))
    BB[2,1] = *(2, S ^ 2, u ^ 2, Gp + *(-1/2, B1p, sinh2G), expB1)
    BB[2,2] = *(S, u ^ 2, *(-2, Sp) + *(S, B1p + *(-2, B2p) + *(4, u)) + *(B1p, S, cosh2G))


    CC[1,1] = *(-2, *(4, Sp ^ 2) + *(-1, B2pp, S ^ 2) + *(-1, Gp ^ 2, S ^ 2) + *(-4, S, Spp) + *(-4, phip ^ 2, S ^ 2) + *(-3, B2p ^ 2, S ^ 2) + *(-1, S, coshGsq, *(B1pp, S) + *(S, B1p ^ 2) + *(3, B1p, Sp)) + *(-3, B2p, S, Sp) + *(-1, B1p, Gp, S ^ 2, sinh2G), expB1)

    CC[1,2] = *(-1, S, *(*(B1pp, S) + *(-1, S, B1p ^ 2) + *(3, B1p, Sp), sinh2G) + *(2, Gpp, S) + *(6, Gp, Sp) + *(-2, B1p, Gp, S) + *(2, B1p, Gp, S, cosh2G))

    CC[2,1] = *(S, *(*(B1pp, S) + *(S, B1p ^ 2) + *(3, B1p, Sp), sinh2G) + *(-6, Gp, Sp) + *(-2, Gpp, S) + *(-2, B1p, Gp, S) + *(2, B1p, Gp, S, cosh2G), expB1)

    CC[2,2] = *(-8, Sp ^ 2) + *(2, B2pp, S ^ 2) + *(2, Gp ^ 2, S ^ 2) + *(6, B2p ^ 2, S ^ 2) + *(8, S, Spp) + *(8, phip ^ 2, S ^ 2) + *(2, S, coshGsq, *(S, B1p ^ 2) + *(-1, B1pp, S) + *(-3, B1p, Sp)) + *(6, B2p, S, Sp) + *(-2, B1p, Gp, S ^ 2, sinh2G)

    SS[1] = *(2, Gp_y, S ^ 2) + *(B1p_y, S ^ 2, sinh2G) + *(-8, S, Sp_x, expB1) + *(-2, B1_y, Gp, S ^ 2) + *(-2, B2p_x, S ^ 2, expB1) + *(6, Gp, S, S_y) + *(8, S_x, Sp, expB1) + *(-1, B1_y, B1p, S ^ 2, sinh2G) + *(-8, phi_x, phip, S ^ 2, expB1) + *(-6, B2_x, B2p, S ^ 2, expB1) + *(-6, B2p, S, S_x, expB1) + *(-2, G_x, Gp, S ^ 2, expB1) + *(-2, S, coshGsq, *(B1p_x, S) + *(3, B1p, S_x) + *(B1_x, B1p, S), expB1) + *(2, B1p, G_y, S ^ 2, cosh2G) + *(3, B1p, S, S_y, sinh2G) + *(-2, B1p, G_x, S ^ 2, expB1, sinh2G)

    SS[2] = *(-8, S, Sp_y) + *(-2, B2p_y, S ^ 2) + *(8, S_y, Sp) + *(-8, phi_y, phip, S ^ 2) + *(-6, B2_y, B2p, S ^ 2) + *(-6, B2p, S, S_y) + *(-2, G_y, Gp, S ^ 2) + *(2, Gp_x, S ^ 2, expB1) + *(2, S, coshGsq, *(B1p_y, S) + *(3, B1p, S_y) + *(-1, B1_y, B1p, S)) + *(-1, B1p_x, S ^ 2, expB1, sinh2G) + *(2, B1_x, Gp, S ^ 2, expB1) + *(2, B1p, G_y, S ^ 2, sinh2G) + *(6, Gp, S, S_x, expB1) + *(-1, B1_x, B1p, S ^ 2, expB1, sinh2G) + *(-3, B1p, S, S_x, expB1, sinh2G) + *(-2, B1p, G_x, S ^ 2, cosh2G, expB1)

    nothing
end


function Sd_outer_eq_coeff!(ABCS::Vector, vars::AllVars)
    u   = vars.u

    B1     = vars.B1
    B1p    = vars.B1p
    B1t    = vars.B1t
    B1h    = vars.B1h
    B1tp   = vars.B1tp
    B1hp   = vars.B1hp

    B2     = vars.B2
    B2p    = vars.B2p
    B2t    = vars.B2t
    B2h    = vars.B2h
    B2tp   = vars.B2tp
    B2hp   = vars.B2hp

    G      = vars.G
    Gp     = vars.Gp
    Gt     = vars.Gt
    Gh     = vars.Gh
    Gtp    = vars.Gtp
    Ghp    = vars.Ghp

    phi    = vars.phi
    phip   = vars.phip
    phit   = vars.phit
    phih   = vars.phih
    phitp  = vars.phitp
    phihp  = vars.phihp

    S      = vars.S
    Sp     = vars.Sp
    St     = vars.St
    Sh     = vars.Sh
    Stp    = vars.Stp
    Shp    = vars.Shp

    Fx     = vars.Fx
    Fxp    = vars.Fxp
    Fxt    = vars.Fxt
    Fxh    = vars.Fxh
    Fxtp   = vars.Fxtp
    Fxhp   = vars.Fxhp

    Fy     = vars.Fy
    Fyp    = vars.Fyp
    Fyt    = vars.Fyt
    Fyh    = vars.Fyh
    Fytp   = vars.Fytp
    Fyhp   = vars.Fyhp


    expB1   = exp(B1)
    sinh2G  = sinh(*(2, G))
    cosh2G  = cosh(*(2, G))
    coshGsq = cosh(G)^2
    coshG   = cosh(G)
    sinhG   = sinh(G)


    ABCS[1] = 0.0

    ABCS[2] = *(-12, S ^ 3, u ^ 2, expB1)

    ABCS[3] = *(24, Sp, S ^ 2, expB1)

    ABCS[4] = *(*(S, *(*(-8, Gh, St) + *(-8, Gt, Sh), coshG) + *(*(-16, Sth) + *(2, Sh, Fxp + *(-4, B2t)) + *(2, Sp, *(4, Fxh) + *(4, Fyt)) + *(2, St, Fyp + *(-4, B2h)), sinhG)) + *(S ^ 2, *(*(-4, Gth) + *(-2, Gh, B1t + B2t + *(-1, Fxp)) + *(2, Gp, Fxh + Fyt) + *(2, Gt, B1h + Fyp + *(-1, B2h)), coshG) + *(*(-4, B2th) + *(2, Fxhp) + *(2, Fytp) + *(-8, phih, phit) + *(-4, Gh, Gt) + *(2, B2h, Fxp + *(-4, B2t)) + *(2, B2p, Fxh + Fyt) + *(2, Fyp, B2t + *(-1, Fxp)), sinhG)) + *(8, Sh, St, sinhG), expB1 *  expB2) + *(*(S, *(*(8, Shh) + *(-8, Fyh, Sp) + *(-2, Sh, Fyp + *(-4, B2h) + *(4, B1h)), coshG) + *(8, Gh, Sh, sinhG)) + *(S ^ 2, *(*(2, Ghh) + *(-2, Fyh, Gp) + *(2, Gh, B2h + *(-1, Fyp) + *(-2, B1h)), sinhG) + *(Fyp ^ 2 + *(-2, B1hh) + *(-2, Fyhp) + *(2, B2hh) + *(2, B1h ^ 2) + *(2, Gh ^ 2) + *(4, B2h ^ 2) + *(4, phih ^ 2) + *(Fyp, *(-2, B2h) + *(2, B1h)) + *(-2, B1h, B2h) + *(-2, B2p, Fyh) + *(2, B1p, Fyh), coshG)) + *(-4, Sh ^ 2, coshG), expB2) + *(*(S, *(*(8, Stt) + *(-8, Fxt, Sp) + *(-2, Fxp, St) + *(8, B1t, St) + *(8, B2t, St), coshG) + *(8, Gt, St, sinhG)) + *(S ^ 2, *(*(2, Gtt) + *(-2, Fxt, Gp) + *(2, Gt, B2t + *(-1, Fxp) + *(2, B1t)), sinhG) + *(Fxp ^ 2 + *(-2, Fxtp) + *(2, B1tt) + *(2, B2tt) + *(2, B1t ^ 2) + *(2, Gt ^ 2) + *(4, B2t ^ 2) + *(4, phit ^ 2) + *(-1, Fxp, *(2, B1t) + *(2, B2t)) + *(-2, Fxt, B1p + B2p) + *(2, B1t, B2t), coshG)) + *(-4, St ^ 2, coshG), exp(B2 + *(2, B1))) + *(8, V(phi), S ^ 4, expB1)

    nothing
end

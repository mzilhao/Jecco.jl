
# Vf(phi)  = VV(phi)
# Vfp(phi) = ∂(VV)(phi)

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

# this is a coupled equation for Fx and Fy. the notation used is
#
# ( A11 d_uu Fx + A12 d_uu Fy + B11 d_u Fx + B12 d_u Fy + C11 Fx + C12 Fy ) = -S1
# ( A21 d_uu Fx + A22 d_uu Fy + B21 d_u Fx + B22 d_u Fy + C21 Fx + C22 Fy ) = -S2

function Fxy_outer_eq_coeff!(ABCS::Vector, vars::AllVars)
    u    = vars.u

    B1   = vars.B1
    B1p  = vars.B1p
    B1h  = vars.B1h
    B1hp = vars.B1hp
    B1t  = vars.B1t
    B1tp = vars.B1tp

    B2p  = vars.B2p
    B2t  = vars.B2t
    B2h  = vars.B2h
    B2tp = vars.B2tp
    B2hp = vars.B2hp

    G    = vars.G
    Gp   = vars.Gp
    Gh   = vars.Gh
    Gt   = vars.Gt
    Ghp  = vars.Ghp
    Gtp  = vars.Gtp

    phih = vars.phih
    phip = vars.phip
    phit = vars.phit

    S    = vars.S
    Sp   = vars.Sp
    Stp  = vars.Stp
    St   = vars.St
    Sh   = vars.Sh
    Shp  = vars.Shp


    AA   = ABCS[1]
    BB   = ABCS[2]
    CC   = ABCS[3]
    SS   = ABCS[4]


    AA[1,1] = *(2, S ^ 2, u ^ 4, exp(B1))
    AA[1,2] = 0
    AA[2,1] = 0
    AA[2,2] = *(2, S ^ 2, u ^ 4)

    BB[1,1] = *(2, S, u ^ 2, *(-1, Sp) + *(-1, B2p, S) + *(2, S, u) + *(-1, B1p, S, cosh(G) ^ 2), exp(B1))
    BB[1,2] = *(2, S ^ 2, u ^ 2, Gp + *(1 / 2, B1p, sinh(*(2, G))))
    BB[2,1] = *(2, S ^ 2, u ^ 2, Gp + *(-1/2, B1p, sinh(*(2, G))), exp(B1))
    BB[2,2] = *(S, u ^ 2, *(-2, Sp) + *(S, B1p + *(-2, B2p) + *(4, u)) + *(B1p, S, cosh(*(2, G))))

    CC[1,1] = 0
    CC[1,2] = 0
    CC[2,1] = 0
    CC[2,2] = 0

    SS[1] = *(S ^ 2, *(2, Ghp) + *(B1hp + *(-1, B1h, B1p), sinh(*(2, G))) + *(-2, B1h, Gp) + *(2, B1p, Gh, cosh(*(2, G)))) + *(*(S, *(-8, Stp) + *(-6, St, B2p + *(B1p, cosh(G) ^ 2))) + *(-2, S ^ 2, B2tp + *(Gp, Gt) + *(cosh(G) ^ 2, B1tp + *(B1p, B1t)) + *(3, B2p, B2t) + *(4, phip, phit) + *(B1p, Gt, sinh(*(2, G)))) + *(8, Sp, St), exp(B1)) + *(6, S, Sh, Gp + *(1 / 2, B1p, sinh(*(2, G))))

    SS[2] = *(S, *(-8, Shp) + *(-6, B2p, Sh) + *(6, B1p, Sh, cosh(G) ^ 2)) + *(-2, S ^ 2, B2hp + *(Gh, Gp) + *(cosh(G) ^ 2, *(-1, B1hp) + *(B1h, B1p)) + *(3, B2h, B2p) + *(4, phih, phip) + *(-1, B1p, Gh, sinh(*(2, G)))) + *(8, Sh, Sp) + *(S, *(S, *(2, Gtp) + *(-1, B1tp + *(B1p, B1t), sinh(*(2, G))) + *(2, B1t, Gp) + *(-2, B1p, Gt, cosh(*(2, G)))) + *(6, Gp, St) + *(-3, B1p, St, sinh(*(2, G))), exp(B1))

    nothing
end

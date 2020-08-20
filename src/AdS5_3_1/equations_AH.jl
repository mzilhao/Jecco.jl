
function AH_eq_coeff(vars::Tuple, ::Outer)
    (
        sigma0, sigma0_x, sigma0_y, sigma0_xx, sigma0_yy, sigma0_xy,
        xi    , xi_x    , xi_y    , xi_xx    , xi_yy, xi_xy,
        B1   , B2   , G   ,  S    , Fx    , Fy    , Sd ,
        B1p  , B2p  , Gp  ,  Sp   , Fxp   , Fyp   , Sdp,
        B1pp , B2pp , Gpp ,  Spp  , Fxpp  , Fypp  ,
        B1_x , B2_x , G_x ,  S_x  , Fx_x  , Fy_x  , Sd_x,
	B1_y , B2_y , G_y ,  S_y  , Fx_y  , Fy_y  , Sd_y,
        B1p_x, B2p_x, Gp_x,  Sp_x , Fxp_x , Fyp_x ,
        B1p_y, B2p_y, Gp_y,  Sp_y , Fxp_y , Fyp_y ,
    ) = vars

    @tilde_outer("B1")
    @tilde_outer("B2")
    @tilde_outer("G")
    @tilde_outer("S")
    @tilde_outer("Fx")
    @tilde_outer("Fy")
    @tilde_outer("Sd")

    @hat_outer("B1")
    @hat_outer("B2")
    @hat_outer("G")
    @hat_outer("S")
    @hat_outer("Fx")
    @hat_outer("Fy")
    @hat_outer("Sd")

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


    x0 = exp(B1 + B2)
    x1 = cosh(G)
    x2 = S*x1
    x3 = sinh(G)
    x4 = exp(B2)
    x5 = x3*x4
    x6 = 2*x5
    x7 = -B1
    x8 = exp(B2 + x7)
    x9 = 3*Sp
    x10 = sigma0_y*x9
    x11 = exp(B1)
    x12 = x1*x11
    x13 = Fy + xi_y
    x14 = Sp*x13
    x15 = sigma0_x*x9
    x16 = Fx + xi_x
    x17 = Sp*x16
    x18 = B2p*x13
    x19 = B2h + x18
    x20 = x19*x3
    x21 = Gp*x13
    x22 = Gh + x21
    x23 = S*x22
    x24 = Fx*Gp
    x25 = Gp*xi_x
    x26 = B1p + B2p
    x27 = x11*(x1*(B1t + B2t + x16*x26) + x3*(Gt + x24 + x25))
    x28 = 2*Fx
    x29 = 2*xi_x
    x30 = B2p*S
    x31 = -Sh
    x32 = 2*Fy
    x33 = 2*xi_y
    x34 = B1p*S
    x35 = S*x11
    x36 = sigma0_y + x13
    x37 = 6*Sp
    x38 = x1*x4
    x39 = 3*Spp
    x40 = x36^2
    x41 = x38*x40
    x42 = 3*B2p
    x43 = sigma0_x + x16
    x44 = exp(2*B1 + B2)
    x45 = x1*x44
    x46 = x43^2
    x47 = x45*x46
    x48 = Gp*x3
    x49 = Sh + x14
    x50 = x3*x49
    x51 = St + x17
    x52 = x12*x51
    x53 = x1*(-Gh - x21) + x27 + x3*(-B2h - x18)
    x54 = S*x53 - x50 + x52
    x55 = x0*x54
    x56 = -x43*x9 + x51
    x57 = x11*x3
    x58 = x56*x57
    x59 = B1h + B1p*x13
    x60 = x1*x19
    x61 = x22*x3
    x62 = B2p*x16 + B2t
    x63 = x57*x62
    x64 = Gp*x16 + Gt
    x65 = x12*x64
    x66 = x1*x59 - x60 - x61 + x63 + x65
    x67 = S*x66 + x1*(-x14 + x31) + x58
    x68 = x4*x67
    x69 = 2*B2p
    x70 = 2*B1p + x69
    x71 = Fyp*x13
    x72 = Fxp*x13
    x73 = 2*sigma0_xy + 2*xi_xy
    x74 = x0*x3
    x75 = 6*S
    x76 = Fyp*x16
    x77 = Fxp*x16
    x78 = x0*x1
    x79 = x11*(Sd*x75 + x5*(2*Fyt + x73 + 2*x76) - x78*(2*Fxt + 2*sigma0_xx + 2*x77 + 2*xi_xx))
    x80 = -x38*(2*Fyh + 2*sigma0_yy + 2*x71 + 2*xi_yy) + x74*(2*Fxh + 2*x72 + x73) + x79
    x81 = 2*x13
    x82 = Fyh + sigma0_yy + x71 + xi_yy
    x83 = sigma0_xy + xi_xy
    x84 = Fxh + x72 + x83
    x85 = 2*x78
    x86 = Fyt + x76 + x83
    x87 = 2*Gp
    x88 = Fxt + sigma0_xx + x77 + xi_xx
    x89 = Spp*x13
    x90 = Gp*x1
    x91 = Spp*x16 + Spt
    x92 = Gp*x57
    x93 = B2ph + B2pp*x13
    x94 = Gph + Gpp*x13
    x95 = B2pp*x16 + B2pt
    x96 = Gpp*x16 + Gpt
    x97 = x0*(2*sigma0_x + x28 + x29)
    x98 = Gp*x12
    x99 = x4*(2*sigma0_y + x32 + x33)
    x100 = S*x80 + x41*x9 + x47*x9 - x54*x97 + x67*x99
    x101 = exp(x7)/2

    axx = -x0*x2

    axy = S*x6

    ayy = -x2*x8

    bx = x4*(S*x20 - S*x27 + Sh*x3 - St*x12 + x1*x23 - x10*x3 + x12*x15 + 2*x12*x17 - 2*x14*x3)

    by = x8*(x1*(B1h*S - B2h*S - Fy*x30 + Fy*x34 + Gt*x35 + Sp*x32 + Sp*x33 + x10 + x24*x35 + x25*x35 - x30*xi_y + x31 + x34*xi_y) + x3*(-x11*(-B2t*S - Fx*x30 + Sp*x28 + Sp*x29 - St + x15 - x30*xi_x) - x23))

    cc = x101*(-B1p*x100 + Fxp*x37*x43*x45 - 2*Fxp*x55 + Fyp*x36*x37*x38 + 2*Fyp*x68 + Gp*x40*x5*x9 + S*(B1p*x79 - Gp*x6*x82 + Gp*x84*x85 + x11*(Sd*x37 + Sdp*x75 - x26*x85*x88 + x38*x86*x87 + x5*x69*x86 + x6*(Fypp*x16 + Fypt) - x74*x87*x88 - x85*(Fxpp*x16 + Fxpt)) - x38*x69*x82 - x38*(2*Fyph + Fypp*x81) + x70*x74*x84 + x74*(2*Fxph + Fxpp*x81)) + Sp*x41*x42 + Sp*x47*(6*B1p + x42) + Sp*x80 + x36*x68*x69 + x39*x41 + x39*x47 - x43*x55*x70 + x44*x46*x48*x9 - x97*(B1p*x52 + S*(B1p*x27 - Gp*x60 - Gp*x61 - x1*x94 + x11*(x1*x95 + x1*(B1pp*x16 + B1pt) + x3*x96 + x48*x62 + x48*(B1p*x16 + B1t) + x64*x90) - x3*x93) + Sp*x53 + x12*x91 - x3*(Sph + x89) - x49*x90 + x51*x92) + x99*(B1p*x58 - Gp*x50 + S*(B1p*x63 + B1p*x65 - Gp*x20 - x1*x93 + x1*(B1ph + B1pp*x13) + x12*x96 - x22*x90 - x3*x94 + x48*x59 + x57*x95 + x62*x98 + x64*x92) + Sp*x66 + x1*(-Sph - x89) + x56*x98 + x57*(-Fxp*x9 - x39*x43 + x91)))

    SS = x100*x101

    return axx, ayy, axy, bx, by, cc, SS
end

function AH_eq_res(vars::Tuple, ::Outer)
    (
        sigma0, sigma0_x, sigma0_y, sigma0_xx, sigma0_yy, sigma0_xy,
        xi    , xi_x    , xi_y    , xi_xx    , xi_yy, xi_xy,
        B1   , B2   , G   ,  S    , Fx    , Fy    , Sd ,
        B1p  , B2p  , Gp  ,  Sp   , Fxp   , Fyp   , Sdp,
        B1pp , B2pp , Gpp ,  Spp  , Fxpp  , Fypp  ,
        B1_x , B2_x , G_x ,  S_x  , Fx_x  , Fy_x  , Sd_x,
        B1_y , B2_y , G_y ,  S_y  , Fx_y  , Fy_y  , Sd_y,
    ) = vars

    @tilde_outer("B1")
    @tilde_outer("B2")
    @tilde_outer("G")
    @tilde_outer("S")
    @tilde_outer("Fx")
    @tilde_outer("Fy")
    @tilde_outer("Sd")

    @hat_outer("B1")
    @hat_outer("B2")
    @hat_outer("G")
    @hat_outer("S")
    @hat_outer("Fx")
    @hat_outer("Fy")
    @hat_outer("Sd")

    x0 = Fy + xi_y
    x1 = cosh(G)
    x2 = exp(B2)
    x3 = x1*x2
    x4 = 3*Sp
    x5 = Fx + xi_x
    x6 = sigma0_x + x5
    x7 = exp(B1 + B2)
    x8 = sinh(G)
    x9 = Sp*x0
    x10 = Sp*x5 + St
    x11 = exp(B1)
    x12 = x1*x11
    x13 = B2p*x0
    x14 = Gp*x0
    x15 = 2*x0
    x16 = 2*sigma0_xy + 2*xi_xy
    x17 = 2*x5
    x18 = x11*x8

    (S*(x11*(6*S*Sd - x1*x7*(Fxp*x17 + 2*Fxt + 2*sigma0_xx + 2*xi_xx) + x2*x8*(Fyp*x17 + 2*Fyt + x16)) - x3*(2*Fyh + Fyp*x15 + 2*sigma0_yy + 2*xi_yy) + x7*x8*(2*Fxh + Fxp*x15 + x16)) + x1*x4*x6^2*exp(2*B1 + B2) + x2*(2*Fy + 2*sigma0_y + 2*xi_y)*(S*(x1*(B1h + B1p*x0) - x1*(B2h + x13) + x12*(Gp*x5 + Gt) + x18*(B2p*x5 + B2t) - x8*(Gh + x14)) + x1*(-Sh - x9) + x18*(x10 - x4*x6)) + x3*x4*(sigma0_y + x0)^2 - x7*(2*Fx + 2*sigma0_x + 2*xi_x)*(S*(x1*(-Gh - x14) + x11*(x1*(B1t + B2t + x5*(B1p + B2p)) + x8*(Fx*Gp + Gp*xi_x + Gt)) + x8*(-B2h - x13)) + x10*x12 - x8*(Sh + x9)))*exp(-B1)/2


end

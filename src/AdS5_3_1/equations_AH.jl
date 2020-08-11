
function AH_eq_coeff(vars::Tuple, ::Outer)
    (
        sigma0, sigma0_x, sigma0_y, sigma0_xx, sigma0_yy, sigma0_xy,
        xi    , xi_x    , xi_y    , xi_xx    , xi_yy,
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
    x9 = B2h*S
    x10 = Fyp*S
    x11 = Sp*sigma0_y
    x12 = 2*xi_y
    x13 = Sp*x12
    x14 = S*xi_y
    x15 = B2p*x14
    x16 = 2*Fy
    x17 = B2p*S
    x18 = 2*sigma0_y
    x19 = x17*x18
    x20 = exp(B1)
    x21 = S*x20
    x22 = Gt*x21
    x23 = Fx*Gp
    x24 = 2*x21
    x25 = x23*x24
    x26 = Gp*sigma0_x
    x27 = x24*x26
    x28 = S*xi_x
    x29 = Gp*x20
    x30 = x28*x29
    x31 = x16 + x18
    x32 = S*(Gh + Gp*(x31 + xi_y))
    x33 = 2*B1p
    x34 = 2*B2p
    x35 = x33 + x34
    x36 = Fx*S
    x37 = S*x33
    x38 = Fx*Sp
    x39 = Sp*sigma0_x
    x40 = 2*xi_x
    x41 = S*x34
    x42 = B2p*x28 + B2t*S + Fxp*S - Sp*x40 + St + sigma0_x*x41 - x38 - x39
    x43 = -Sh
    x44 = Fy + sigma0_y
    x45 = x44 + xi_y
    x46 = 6*Sp
    x47 = x1*x4
    x48 = 3*Spp
    x49 = x45^2
    x50 = x47*x49
    x51 = 3*B2p
    x52 = 3*Sp
    x53 = Fx + sigma0_x
    x54 = x53 + xi_x
    x55 = exp(2*B1 + B2)
    x56 = x1*x55
    x57 = x54^2
    x58 = x56*x57
    x59 = Gp*x3
    x60 = Sp*x44
    x61 = Sh + x60
    x62 = x3*x61
    x63 = St + x38 + x39
    x64 = x1*x20
    x65 = x63*x64
    x66 = B2p*x44
    x67 = Gp*x44
    x68 = Gt + x23 + x26
    x69 = x3*x68
    x70 = B1p + B2p
    x71 = x20*(x1*(B1t + B2t + x53*x70) + x69)
    x72 = x1*(-Gh - x67) + x3*(-B2h - x66) + x71
    x73 = S*x72 - x62 + x65
    x74 = x0*x73
    x75 = 2*sigma0_xy
    x76 = Fxp*x44
    x77 = x0*x3
    x78 = 6*S
    x79 = Fx*Fyp
    x80 = Fyp*sigma0_x
    x81 = Fx*Fxp
    x82 = Fxp*sigma0_x
    x83 = x0*x1
    x84 = x20*(Sd*x78 + x5*(2*Fyt + x75 + 2*x79 + 2*x80) - x83*(2*Fxt + 2*sigma0_xx + 2*x81 + 2*x82))
    x85 = -x47*(2*Fyh + Fyp*x16 + Fyp*x18 + 2*sigma0_yy) + x77*(2*Fxh + x75 + 2*x76) + x84
    x86 = St - 2*x38 - 2*x39 - x52*xi_x
    x87 = x20*x3
    x88 = x86*x87
    x89 = B1h + B1p*x44
    x90 = B2h + x66
    x91 = x1*x90
    x92 = Gh + x67
    x93 = x3*x92
    x94 = B2p*x53 + B2t
    x95 = x87*x94
    x96 = x64*x68
    x97 = x1*x89 - x91 - x93 + x95 + x96
    x98 = S*x97 + x1*(x43 - x60) + x88
    x99 = x4*x98
    x100 = 2*x44
    x101 = Fy*Fyp + Fyh + Fyp*sigma0_y + sigma0_yy
    x102 = Fxh + sigma0_xy + x76
    x103 = 2*x83
    x104 = Fyt + sigma0_xy + x79 + x80
    x105 = 2*Gp
    x106 = Fxt + sigma0_xx + x81 + x82
    x107 = Sph + Spp*x44
    x108 = Gp*x1
    x109 = B2ph + B2pp*x44
    x110 = Gph + Gpp*x44
    x111 = Gp^2
    x112 = Fx*Gpp + Gpp*sigma0_x + Gpt
    x113 = 2*Fx
    x114 = 2*sigma0_x
    x115 = x0*(x113 + x114 + x40)
    x116 = Gp*x64
    x117 = x4*(x12 + x31)
    x118 = S*x85 - x115*x73 + x117*x98 + x50*x52 + x52*x58
    x119 = exp(x7)/2

    axx = -x0*x2

    axy = S*x6

    ayy = -x2*x8

    bx = x4*(x1*(-x20*(B1p*x28 + B1t*S + sigma0_x*x37 + x35*x36 + x42) + x32) + x3*(-Fy*Sp + Sh + x10 - x11 - x13 + x15 + x16*x17 + x19 - x22 - x25 - x27 - x30 + x9))

    by = x8*(x1*(B1h*S + B1p*S*x18 + B1p*x14 + Fy*(Sp + x37 - x41) - x10 + x11 + x13 - x15 - x19 + x22 + x25 + x27 + x30 + x43 - x9) + x3*(x20*(x34*x36 + x42) - x32))

    cc = x119*(-B1p*x118 + Fxp*x46*x54*x56 - 2*Fxp*x74 + Fyp*x45*x46*x47 + 2*Fyp*x99 + Gp*x49*x5*x52 + S*(B1p*x84 - Gp*x101*x6 + Gp*x102*x103 - x101*x34*x47 + x102*x35*x77 + x20*(Sd*x46 + Sdp*x78 - x103*x106*x70 - x103*(Fx*Fxpp + Fxpp*sigma0_x + Fxpt) + x104*x105*x47 + x104*x34*x5 - x105*x106*x77 + x6*(Fx*Fypp + Fypp*sigma0_x + Fypt)) - x47*(2*Fyph + Fypp*x100) + x77*(2*Fxph + Fxpp*x100)) + Sp*x50*x51 + Sp*x58*(6*B1p + x51) + Sp*x85 - x115*(B1p*x65 + Gp*x63*x87 + S*(B1p*x71 - Gp*x91 - Gp*x93 - x1*x110 - x109*x3 + x20*(x1*(B1pp*Fx + B1pp*sigma0_x + B1pt + B2pp*Fx + B2pp*sigma0_x + B2pt + Fx*x111 + Gp*Gt + sigma0_x*x111) + x3*(B1p*x23 + B1p*x26 + B1t*Gp + B2p*x23 + B2p*x26 + B2t*Gp + x112))) + Sp*x72 - x107*x3 - x108*x61 + x64*(Fx*Spp + Spp*sigma0_x + Spt)) + x117*(B1p*x88 - Gp*x62 + S*(B1p*x95 + B1p*x96 - x1*x109 + x1*(B1ph + B1pp*x44) - x108*x92 - x110*x3 + x112*x64 + x116*x94 + x29*x69 + x59*x89 - x59*x90 + x87*(B2pp*x53 + B2pt)) + Sp*x97 - x1*x107 + x116*x86 - x87*(Fxp*x52 + Spp*x113 + Spp*x114 - Spt + x48*xi_x)) + x34*x45*x99 - x35*x54*x74 + x48*x50 + x48*x58 + x52*x55*x57*x59)

    # SS = x118*x119

    return axx, ayy, axy, bx, by, cc
end

function AH_eq_res(vars::Tuple, ::Outer)
    (
        sigma0, sigma0_x, sigma0_y, sigma0_xx, sigma0_yy, sigma0_xy,
        xi    , xi_x    , xi_y    , xi_xx    , xi_yy,
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

    x0 = Fy + sigma0_y
    x1 = cosh(G)
    x2 = exp(B2)
    x3 = x1*x2
    x4 = 3*Sp
    x5 = Fx + sigma0_x
    x6 = sinh(G)
    x7 = exp(B1 + B2)
    x8 = 2*sigma0_xy
    x9 = 2*Fy
    x10 = 2*sigma0_y
    x11 = exp(B1)
    x12 = 2*Fx
    x13 = 2*sigma0_x
    x14 = Sp*x0
    x15 = x1*x11
    x16 = B2p*x0
    x17 = Gp*x0
    x18 = Fx*Gp + Gp*sigma0_x + Gt
    x19 = x11*x6

    (S*(x11*(6*S*Sd - x1*x7*(Fxp*x12 + Fxp*x13 + 2*Fxt + 2*sigma0_xx) + x2*x6*(Fyp*x12 + Fyp*x13 + 2*Fyt + x8)) - x3*(2*Fyh + Fyp*x10 + Fyp*x9 + 2*sigma0_yy) + x6*x7*(2*Fxh + 2*Fxp*x0 + x8)) + x1*x4*(x5 + xi_x)^2*exp(2*B1 + B2) + x2*(x10 + x9 + 2*xi_y)*(S*(x1*(B1h + B1p*x0) - x1*(B2h + x16) + x15*x18 + x19*(B2p*x5 + B2t) - x6*(Gh + x17)) + x1*(-Sh - x14) + x19*(-Sp*x12 - Sp*x13 + St - x4*xi_x)) + x3*x4*(x0 + xi_y)^2 - x7*(x12 + x13 + 2*xi_x)*(S*(x1*(-Gh - x17) + x11*(x1*(B1t + B2t + x5*(B1p + B2p)) + x18*x6) + x6*(-B2h - x16)) + x15*(Fx*Sp + Sp*sigma0_x + St) - x6*(Sh + x14)))*exp(-B1)/2

end

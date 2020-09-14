
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

    x0 = cosh(G)
    x1 = exp(B1 + B2)
    x2 = x0*x1
    x3 = sigma0^2
    x4 = 1/x3
    x5 = S*x4
    x6 = sinh(G)
    x7 = exp(B2)
    x8 = -B1
    x9 = exp(B2 + x8)
    x10 = sigma0^(-4)
    x11 = 2*B2p
    x12 = S*x11
    x13 = Sp - x12
    x14 = Fy*x3
    x15 = Sh*x3
    x16 = Sp*sigma0_y
    x17 = B2h*x3
    x18 = S*x3
    x19 = x3*xi_y
    x20 = S*sigma0
    x21 = 4*x20
    x22 = Gt*x3
    x23 = exp(B1)
    x24 = S*x23
    x25 = 2*sigma0_x
    x26 = Gp*x25
    x27 = Fx*x3
    x28 = 2*Gp
    x29 = x24*x28
    x30 = x3*xi_x
    x31 = Fyp*x18 + S*x17 - Sp*x19 - sigma0_y*x12 - sigma0_y*x21 + x12*x19 + x15 + x16 - x22*x24 + x24*x26 - x27*x29 - x29*x30
    x32 = Gh*x3
    x33 = -sigma0_y + x14 + x19
    x34 = S*(x28*x33 + x32)
    x35 = Sp*sigma0_x
    x36 = St*x3
    x37 = 2*B1p
    x38 = S*x37
    x39 = sigma0_x*x12
    x40 = sigma0_x*x21
    x41 = B1t*x3
    x42 = B2t*x3
    x43 = S*x42
    x44 = Fxp*x18
    x45 = Sp*x30
    x46 = x12*x30
    x47 = x11 + x37
    x48 = -x35
    x49 = -x36
    x50 = B1h*x3
    x51 = x33^2
    x52 = x0*x7
    x53 = x51*x52
    x54 = 3*Spp
    x55 = 3*B2p
    x56 = Gp*x6
    x57 = 3*Sp
    x58 = x56*x57
    x59 = -sigma0_x + x27 + x30
    x60 = x59^2
    x61 = 2*B1
    x62 = exp(B2 + x61)
    x63 = x0*x62
    x64 = x60*x63
    x65 = 2*sigma0_y
    x66 = -x65
    x67 = Fyp*sigma0
    x68 = x66 + x67
    x69 = sigma0*x52
    x70 = 6*Sp
    x71 = -x25
    x72 = Fxp*sigma0
    x73 = x71 + x72
    x74 = sigma0^3
    x75 = 3*x74
    x76 = Sd*x24
    x77 = sigma0*sigma0_xy
    x78 = -2*x77
    x79 = x74*xi_xy
    x80 = sigma0_x*sigma0_y
    x81 = Fyp*x74
    x82 = Fxp*x74
    x83 = Fyt*x74 - sigma0_x*x67
    x84 = Fxh*x74 - sigma0_y*x72
    x85 = x23*x6
    x86 = x85*(Fx*x81 + Fy*x82 + x78 + 2*x79 + 4*x80 + x81*xi_x + x82*xi_y + x83 + x84)
    x87 = sigma0_x^2
    x88 = sigma0*sigma0_xx
    x89 = Fxt*x74 - sigma0_x*x72 + x74*xi_xx + 2*x87 - x88
    x90 = sigma0_y^2
    x91 = sigma0*sigma0_yy
    x92 = Fyh*x74 - sigma0_y*x67 + x74*xi_yy + 2*x90 - x91
    x93 = x0*(Fy*x81 + x81*xi_y + x92 + (Fx*x82 + x82*xi_x + x89)*exp(x61))
    x94 = x7*(-x86 + x93) - x75*x76
    x95 = 2*Sp
    x96 = 2*sigma0
    x97 = Fy + xi_y
    x98 = x3*x97
    x99 = Sp*x98 + x15 - x16
    x100 = x6*x99
    x101 = Fx + xi_x
    x102 = Sp*x101
    x103 = x102*x3 + x36 + x48
    x104 = x0*x23
    x105 = x103*x104
    x106 = -B2p*sigma0_y + B2p*x98 + x17
    x107 = x106*x6
    x108 = Gp*sigma0_y
    x109 = Gp*x98 - x108 + x32
    x110 = x0*x109
    x111 = x101*x3
    x112 = -B1p*sigma0_x + B1p*x111 + x41
    x113 = -B2p*sigma0_x + B2p*x111 + x42
    x114 = x0*x113
    x115 = Gp*sigma0_x
    x116 = Gp*x101
    x117 = -x115 + x116*x3 + x22
    x118 = x117*x6
    x119 = x23*(x0*x112 + x114 + x118)
    x120 = x107 + x110 - x119
    x121 = -S*x120 - x100 + x105
    x122 = x1*x121
    x123 = x0*x99
    x124 = 2*x35
    x125 = 2*x27
    x126 = Sp*x125
    x127 = -B1p*sigma0_y + B1p*x98 + x50
    x128 = x0*x106
    x129 = x109*x6
    x130 = x113*x85
    x131 = x104*x117
    x132 = x0*x127 - x128 - x129 + x130 + x131
    x133 = -S*x132 + x123 + x85*(-x124 + x126 + 2*x45 + x49)
    x134 = x133*x7
    x135 = -Spp*sigma0_y - x16*x96 + x3*(Sph + Spp*x97)
    x136 = sigma0*x11
    x137 = -B2pp*sigma0_y - sigma0_y*x136 + x3*(B2ph + B2pp*x97)
    x138 = -Gpp*sigma0_y - x108*x96 + x3*(Gph + Gpp*x97)
    x139 = Gp*x0
    x140 = sigma0*x37
    x141 = -B2pp*sigma0_x - sigma0_x*x136 + x3*(B2pp*x101 + B2pt)
    x142 = -Gpp*sigma0_x - sigma0*x26 + x3*(Gpp*x101 + Gpt)
    x143 = 2*x30
    x144 = x1*(x125 + x143 + x71)
    x145 = -x126 + x36 + x95*(sigma0_x - x30)
    x146 = Gp*x23
    x147 = x7*(2*x14 + 2*x19 + x66)
    x148 = x78 + 6*x80
    x149 = sigma0*x6
    x150 = x7*(x81*x97 + x92)
    x151 = sigma0_y*x25 - x77 + x79
    x152 = x1*(x151 + x82*x97 + x84)
    x153 = B1p + B2p
    x154 = S*x75
    x155 = x7*(x101*x81 + x151 + x83)
    x156 = x155*x6
    x157 = x101*x82 + x89
    x158 = x157*x2
    x159 = 2*x20
    x160 = exp(x8)/2
    x161 = -sigma0_y*x4 + x97
    x162 = -sigma0_x*x4 + x101
    x163 = Sp*x161
    x164 = 2*Fx
    x165 = 2*xi_x
    x166 = B2p*x161
    x167 = Gp*x161
    x168 = B2p*x162 + B2t
    x169 = Gt - x115*x4 + x116

    axx = x2*x5

    axy = -2*x5*x6*x7

    ayy = x0*x5*x9

    bx = x10*x7*(x0*(x23*(S*x41 - sigma0_x*x38 + x27*(S*x47 - Sp) + x30*x38 + x35 + x36 - x39 - x40 + x43 + x44 - x45 + x46) - x34) - x6*(-x13*x14 + x31))

    by = x10*x9*(x0*(-S*x50 + sigma0_y*x38 - x14*(x13 + x38) - x19*x38 + x31) + x6*(x23*(x13*x27 + x39 + x40 - x43 - x44 + x45 - x46 + x48 + x49) + x34))

    cc = x160*(B1p*(-x121*x144 - x133*x147 - x159*x94 + x53*x57 + x57*x64) - Sp*x53*x55 - Sp*x64*(6*B1p + x55) - sigma0*x59*x63*x70*x73 + sigma0*x94*x95 + x11*x134*x33 + x122*x47*x59 + x122*x73*x96 + x134*x68*x96 + x144*(B1p*x105 + Gp*x103*x85 - Gp*x123 - S*(-B1p*x119 + Gp*x128 + Gp*x129 + x0*x138 + x137*x6 - x23*(x0*x141 + x0*(-B1pp*sigma0_x - sigma0_x*x140 + x3*(B1pp*x101 + B1pt)) + x112*x56 + x113*x56 + x117*x139 + x142*x6)) - Sp*x120 + x104*(-Spp*sigma0_x - sigma0*x124 + x3*(Spp*x101 + Spt)) - x135*x6) - x147*(B1p*x145*x85 - Gp*x100 + Gp*x104*x145 + S*(B1p*x130 + B1p*x131 - Gp*x107 - Gp*x110 - x0*x137 + x0*(-B1pp*sigma0_y - sigma0_y*x140 + x3*(B1ph + B1pp*x97)) + x104*x142 + x114*x146 + x118*x146 + x127*x56 - x138*x6 + x141*x85) + Sp*x132 - x0*x135 + x85*(-Fxp*x3*x57 - Spp*x125 - Spp*x143 + Spp*x25 + Spt*x3 + 4*sigma0*x35)) - x159*(-B1p*x23*(-Sd*x154 - x156 + x158) - B2p*x0*x150 + x1*x149*(-Fxpp*sigma0_y + x148 + x3*(Fxph + Fxpp*x97) - x65*x72) + x139*x152 - x150*x56 + x152*x153*x6 + x23*(B2p*x156 + Sd*x57*x74 + Sdp*x154 - sigma0*x2*(-Fxpp*sigma0_x - x25*x72 + x3*(Fxpp*x101 + Fxpt) + 6*x87 - 2*x88) - x1*x157*x56 + x139*x155 + x149*x7*(-Fypp*sigma0_x + x148 - x25*x67 + x3*(Fypp*x101 + Fypt)) - x153*x158) - x69*(-Fypp*sigma0_y + x3*(Fyph + Fypp*x97) - x65*x67 + 6*x90 - 2*x91)) - x33*x68*x69*x70 - x51*x58*x7 - x53*x54 - x54*x64 - x58*x60*x62)/sigma0^6

    SS = x160*(S*(-x7*(-2*x86 + 2*x93)/x74 + 6*x76) - x1*(x164 + x165 - x25*x4)*(S*(x0*(-Gh - x167) + x23*(x0*x168 + x0*(B1p*x162 + B1t) + x169*x6) + x6*(-B2h - x166)) + x104*(St + x102 - x35*x4) - x6*(Sh + x163)) + x161^2*x52*x57 + x162^2*x57*x63 + x7*(2*Fy - x4*x65 + 2*xi_y)*(S*(x0*(B1h + B1p*x161) - x0*(B2h + x166) + x104*x169 + x168*x85 - x6*(Gh + x167)) + x0*(-Sh - x163) + x85*(-Sp*x164 - Sp*x165 + St + x124*x4)))

    return axx, ayy, axy, bx, by, cc, SS
end

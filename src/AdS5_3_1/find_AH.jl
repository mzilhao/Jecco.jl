
## FIXME

#= Finding the Apparent Horizon

this is a 2D PDE of the type

  (axx Dxx + ayy Dyy + Axy Dxy + bx Dx + by Dy + c Id) sigma = -S

we first build each operator (Dxx, Dyy, etc) through a Kronecker product and
then overwrite each of them with the corresponding coefficient. we then sum them
all up and solve the linear system.

=#

## FIXME
function find_AH!(sigma::Array, bulkconstrain::BulkConstrained,
                  bulkevol::BulkEvolved, deriv::BulkDeriv, gauge::Gauge,
                  cache::HorizonCache, sys::System{Outer},
                  gaugecondition::ConstantAH)
    _, Nx, Ny = size(sys)
    bulk = Bulk(bulkevol, bulkconstrain)

    Du  = sys.Du
    Dx  = sys.Dx
    Dxx = sys.Dxx
    Dy  = sys.Dy
    Dyy = sys.Dyy

    interp = sys.uinterp

    B1_uAH      = cache.bulkhorizon.B1_uAH
    B2_uAH      = cache.bulkhorizon.B2_uAH
    G_uAH       = cache.bulkhorizon.G_uAH
    S_uAH       = cache.bulkhorizon.S_uAH
    Fx_uAH      = cache.bulkhorizon.Fx_uAH
    Fy_uAH      = cache.bulkhorizon.Fy_uAH
    Sd_uAH      = cache.bulkhorizon.Sd_uAH

    Du_B1_uAH   = cache.bulkhorizon.Du_B1_uAH
    Du_B2_uAH   = cache.bulkhorizon.Du_B2_uAH
    Du_G_uAH    = cache.bulkhorizon.Du_G_uAH
    Du_S_uAH    = cache.bulkhorizon.Du_S_uAH
    Du_Fx_uAH   = cache.bulkhorizon.Du_Fx_uAH
    Du_Fy_uAH   = cache.bulkhorizon.Du_Fy_uAH
    Du_Sd_uAH   = cache.bulkhorizon.Du_Sd_uAH

    Duu_B1_uAH  = cache.bulkhorizon.Duu_B1_uAH
    Duu_B2_uAH  = cache.bulkhorizon.Duu_B2_uAH
    Duu_G_uAH   = cache.bulkhorizon.Duu_G_uAH
    Duu_S_uAH   = cache.bulkhorizon.Duu_S_uAH
    Duu_Fx_uAH  = cache.bulkhorizon.Duu_Fx_uAH
    Duu_Fy_uAH  = cache.bulkhorizon.Duu_Fy_uAH

    axx         = cache.axx
    ayy         = cache.ayy
    axy         = cache.axy
    bx          = cache.bx
    by          = cache.by
    cc          = cache.cc
    b_vec       = cache.b_vec

    Dx_2D       = cache.Dx_2D
    Dxx_2D      = cache.Dxx_2D
    Dy_2D       = cache.Dy_2D
    Dyy_2D      = cache.Dyy_2D
    Dxy_2D      = cache.Dxy_2D
    _Dx_2D      = cache._Dx_2D
    _Dxx_2D     = cache._Dxx_2D
    _Dy_2D      = cache._Dy_2D
    _Dyy_2D     = cache._Dyy_2D
    _Dxy_2D     = cache._Dxy_2D


    uAH = gaugecondition.u_AH
    u2  = uAH * uAH
    u3  = uAH * uAH * uAH
    u4  = uAH * uAH * uAH * uAH

    # take all u-derivatives. we assume that they've not been computed beforehand
    @sync begin
        @spawn mul!(deriv.Du_B1,  Du,  bulk.B1)
        @spawn mul!(deriv.Du_B2,  Du,  bulk.B2)
        @spawn mul!(deriv.Du_G,   Du,  bulk.G)
        @spawn mul!(deriv.Du_S,   Du,  bulk.S)
        @spawn mul!(deriv.Du_Fx,  Du,  bulk.Fx)
        @spawn mul!(deriv.Du_Fy,  Du,  bulk.Fy)
        @spawn mul!(deriv.Du_Sd,  Du,  bulk.Sd)
    end

    # interpolate bulk functions (and u-derivatives) to the u=uAH surface
    @inbounds Threads.@threads for j in 1:Ny
        @inbounds for i in 1:Nx
            B1_uAH[1,i,j]       = interp(view(bulk.B1,  :,i,j))(uAH)
            B2_uAH[1,i,j]       = interp(view(bulk.B2,  :,i,j))(uAH)
            G_uAH[1,i,j]        = interp(view(bulk.G,   :,i,j))(uAH)
            S_uAH[1,i,j]        = interp(view(bulk.S,   :,i,j))(uAH)
            Fx_uAH[1,i,j]       = interp(view(bulk.Fx,  :,i,j))(uAH)
            Fy_uAH[1,i,j]       = interp(view(bulk.Fy,  :,i,j))(uAH)
            Sd_uAH[1,i,j]       = interp(view(bulk.Sd,  :,i,j))(uAH)

            Du_B1_uAH[1,i,j]    = interp(view(deriv.Du_B1,  :,i,j))(uAH)
            Du_B2_uAH[1,i,j]    = interp(view(deriv.Du_B2,  :,i,j))(uAH)
            Du_G_uAH[1,i,j]     = interp(view(deriv.Du_G,   :,i,j))(uAH)
            Du_S_uAH[1,i,j]     = interp(view(deriv.Du_S,   :,i,j))(uAH)
            Du_Fx_uAH[1,i,j]    = interp(view(deriv.Du_Fx,  :,i,j))(uAH)
            Du_Fy_uAH[1,i,j]    = interp(view(deriv.Du_Fy,  :,i,j))(uAH)
            Du_Sd_uAH[1,i,j]    = interp(view(deriv.Du_Sd,  :,i,j))(uAH)

            Duu_B1_uAH[1,i,j]   = interp(view(deriv.Duu_B1,  :,i,j))(uAH)
            Duu_B2_uAH[1,i,j]   = interp(view(deriv.Duu_B2,  :,i,j))(uAH)
            Duu_G_uAH[1,i,j]    = interp(view(deriv.Duu_G,   :,i,j))(uAH)
            Duu_S_uAH[1,i,j]    = interp(view(deriv.Duu_S,   :,i,j))(uAH)
            Duu_Fx_uAH[1,i,j]   = interp(view(deriv.Duu_Fx,  :,i,j))(uAH)
            Duu_Fy_uAH[1,i,j]   = interp(view(deriv.Duu_Fy,  :,i,j))(uAH)
        end
    end

    ind2D  = LinearIndices(B1_uAH[1,:,:])

    # coefficients of the derivative operators
    @fastmath @inbounds Threads.@threads for j in 1:Ny
        @inbounds for i in 1:Nx
            idx   = ind2D[i,j]

            xi    = gauge.xi[1,i,j]
            xi_x  = Dx(gauge.xi, 1,i,j)
            xi_y  = Dy(gauge.xi, 1,i,j)
            xi_xx = Dxx(gauge.xi, 1,i,j)
            xi_yy = Dyy(gauge.xi, 1,i,j)

            sigma0    = sigma[1,i,j]
            sigma0_x  = Dx(sigma, 1,i,j)
            sigma0_y  = Dy(sigma, 1,i,j)
            sigma0_xx = Dxx(sigma, 1,i,j)
            sigma0_yy = Dyy(sigma, 1,i,j)
            sigma0_xy = Dx(Dy, sigma, 1,i,j)

            B1    = B1_uAH[1,i,j]
            B2    = B2_uAH[1,i,j]
            G     = G_uAH[1,i,j]
            S     = S_uAH[1,i,j]
            Fx    = Fx_uAH[1,i,j]
            Fy    = Fy_uAH[1,i,j]
            Sd    = Sd_uAH[1,i,j]

            # r derivatives

            B1p        = -u2 * Du_B1_uAH[1,i,j]
            B2p        = -u2 * Du_B2_uAH[1,i,j]
            Gp         = -u2 * Du_G_uAH[1,i,j]
            Sp         = -u2 * Du_S_uAH[1,i,j]
            Fxp        = -u2 * Du_Fx_uAH[1,i,j]
            Fyp        = -u2 * Du_Fy_uAH[1,i,j]
            Sdp        = -u2 * Du_Sd_uAH[1,i,j]

            B1pp       = 2*u3 * Du_B1_uAH[1,i,j]  + u4 * Duu_B1_uAH[1,i,j]
            B2pp       = 2*u3 * Du_B2_uAH[1,i,j]  + u4 * Duu_B2_uAH[1,i,j]
            Gpp        = 2*u3 * Du_G_uAH[1,i,j]   + u4 * Duu_G_uAH[1,i,j]
            Spp        = 2*u3 * Du_S_uAH[1,i,j]   + u4 * Duu_S_uAH[1,i,j]
            Fxpp       = 2*u3 * Du_Fx_uAH[1,i,j]  + u4 * Duu_Fx_uAH[1,i,j]
            Fypp       = 2*u3 * Du_Fy_uAH[1,i,j]  + u4 * Duu_Fy_uAH[1,i,j]

            # x and y derivatives

            B1_x    = Dx(B1_uAH,   1,i,j)
            B2_x    = Dx(B2_uAH,   1,i,j)
            G_x     = Dx(G_uAH,    1,i,j)
            S_x     = Dx(S_uAH,    1,i,j)
            Fx_x    = Dx(Fx_uAH,   1,i,j)
            Fy_x    = Dx(Fy_uAH,   1,i,j)
            Sd_x    = Dx(Sd_uAH,   1,i,j)

            B1_y    = Dy(B1_uAH,   1,i,j)
            B2_y    = Dy(B2_uAH,   1,i,j)
            G_y     = Dy(G_uAH,    1,i,j)
            S_y     = Dy(S_uAH,    1,i,j)
            Fx_y    = Dy(Fx_uAH,   1,i,j)
            Fy_y    = Dy(Fy_uAH,   1,i,j)
            Sd_y    = Dy(Sd_uAH,   1,i,j)

            B1p_x   = -u2 * Dx(Du_B1_uAH, 1,i,j)
            B2p_x   = -u2 * Dx(Du_B2_uAH, 1,i,j)
            Gp_x    = -u2 * Dx(Du_G_uAH,  1,i,j)
            Sp_x    = -u2 * Dx(Du_S_uAH,  1,i,j)
            Fxp_x   = -u2 * Dx(Du_Fx_uAH, 1,i,j)
            Fyp_x   = -u2 * Dx(Du_Fy_uAH, 1,i,j)

            B1p_y   = -u2 * Dy(Du_B1_uAH, 1,i,j)
            B2p_y   = -u2 * Dy(Du_B2_uAH, 1,i,j)
            Gp_y    = -u2 * Dy(Du_G_uAH,  1,i,j)
            Sp_y    = -u2 * Dy(Du_S_uAH,  1,i,j)
            Fxp_y   = -u2 * Dy(Du_Fx_uAH, 1,i,j)
            Fyp_y   = -u2 * Dy(Du_Fy_uAH, 1,i,j)

            vars = (
                sigma0, sigma0_x, sigma0_y, sigma0_xx, sigma0_yy, sigma0_xy,
                xi    , xi_x    , xi_y    , xi_xx    , xi_yy,
                B1   , B2   , G   ,  S    , Fx    , Fy    , Sd ,
                B1p  , B2p  , Gp  ,  Sp   , Fxp   , Fyp   , Sdp,
                B1pp , B2pp , Gpp ,  Spp  , Fxpp  , Fypp  ,
                B1_x , B2_x , G_x ,  S_x  , Fx_x  , Fy_x  , Sd_x,
	        B1_y , B2_y , G_y ,  S_y  , Fx_y  , Fy_y  , Sd_y,
                B1p_x, B2p_x, Gp_x,  Sp_x , Fxp_x , Fyp_x ,
                B1p_y, B2p_y, Gp_y,  Sp_y , Fxp_y , Fyp_y ,
            )

            a11, a22, a12, b1, b2, c = AH_eq_coeff(vars, sys.gridtype)

            axx[idx]   = a11
            ayy[idx]   = a22
            axy[idx]   = a12
            bx[idx]    = b1
            by[idx]    = b2
            cc[idx]    = c

            # TODO: res
            # b_vec[idx] = -SS
            b_vec[idx] = -0
        end
    end

    # each time this routine is called, the operators Dx_2D, Dxx_2D, etc, are
    # overwritten. so restore them here from _Dx_2D, _Dxx_2D, etc. these are
    # never overwritten.
    copyto!(Dx_2D,  _Dx_2D)
    copyto!(Dxx_2D, _Dxx_2D)
    copyto!(Dy_2D,  _Dy_2D)
    copyto!(Dyy_2D, _Dyy_2D)
    copyto!(Dxy_2D, _Dxy_2D)

    # overwrite the operators with the coefficients computed in the loop above
    mul_col!(axx, Dxx_2D)
    mul_col!(ayy, Dyy_2D)
    mul_col!(axy, Dxy_2D)
    mul_col!(bx,  Dx_2D)
    mul_col!(by,  Dy_2D)
    ccId = Diagonal(cc)

    # build actual operator to be inverted
    A_mat = Dxx_2D + Dyy_2D + Dxy_2D + Dx_2D + Dy_2D + ccId

    # since we're using periodic boundary conditions, the operator A_mat (just
    # like the Dx and Dxx operators) is strictly speaking not invertible (it has
    # zero determinant) since the solution is not unique. indeed, its LU
    # decomposition shouldn't even be defined. for some reason, however, the
    # call to "lu" does in fact factorize the matrix. in any case, to be safer,
    # let's instead use the left division operator. this calls "factorize",
    # which uses fancy algorithms to determine which is the best way to
    # factorize (and which performs a QR decomposition if the LU fails). the
    # inverse that is performed probably returns the minimum norm least squares
    # solution, or something similar. in any case, for our purposes here we
    # mostly care about getting a solution (not necessarily the minimum norm
    # least squares one).
    sol = A_mat \ b_vec

    # FIXME
    # copyto!(xi_t, sol)

    nothing
end

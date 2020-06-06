
function compute_xi_t!(gauge_t::Gauge, bulkconstrain::BulkConstrained, bulkevol::BulkEvolved,
                       deriv::BulkDeriv, gauge::Gauge, sys::System{Outer}, ::ConstantGauge)
    xi_t = getxi(gauge_t)
    fill!(xi_t, 0)
    nothing
end

# TODO
function compute_xi_t!(gauge_t::Gauge, bulkconstrain::BulkConstrained, bulkevol::BulkEvolved,
                       deriv::BulkDeriv, gauge::Gauge, sys::System{Outer}, gaugecondition::ConstantAH)
    _, Nx, Ny = size(sys)
    bulk = Bulk(bulkevol, bulkconstrain)

    Du  = sys.Du
    Dx  = sys.Dx
    Dxx = sys.Dxx
    Dy  = sys.Dy
    Dyy = sys.Dyy

    interp = sys.uinterp

    # take u-derivatives of Sd, B1d, B2d and Gd. these are not needed in the
    # nested system, so they haven't been computed before. since we need them
    # here, compute them now
    @sync begin
        @spawn mul!(deriv.Du_Sd,  Du,  bulkconstrain.Sd)
        @spawn mul!(deriv.Du_B1d, Du,  bulkconstrain.B1d)
        @spawn mul!(deriv.Du_B2d, Du,  bulkconstrain.B2d)
        @spawn mul!(deriv.Du_Gd,  Du,  bulkconstrain.Gd)
    end

    xi_t = getxi(gauge_t)

    uAH   = gaugecondition.u_AH
    kappa = gaugecondition.kappa

    u2 = uAH * uAH
    u3 = uAH * uAH * uAH
    u4 = uAH * uAH * uAH * uAH

    T = typeof(uAH)
    B1_uAH      = zeros(T,1,Nx,Ny)
    B2_uAH      = zeros(T,1,Nx,Ny)
    G_uAH       = zeros(T,1,Nx,Ny)
    phi_uAH     = zeros(T,1,Nx,Ny)
    S_uAH       = zeros(T,1,Nx,Ny)
    Fx_uAH      = zeros(T,1,Nx,Ny)
    Fy_uAH      = zeros(T,1,Nx,Ny)
    Sd_uAH      = zeros(T,1,Nx,Ny)
    B1d_uAH     = zeros(T,1,Nx,Ny)
    B2d_uAH     = zeros(T,1,Nx,Ny)
    Gd_uAH      = zeros(T,1,Nx,Ny)
    phid_uAH    = zeros(T,1,Nx,Ny)
    A_uAH       = zeros(T,1,Nx,Ny)

    Du_B1_uAH   = zeros(T,1,Nx,Ny)
    Du_B2_uAH   = zeros(T,1,Nx,Ny)
    Du_G_uAH    = zeros(T,1,Nx,Ny)
    Du_phi_uAH  = zeros(T,1,Nx,Ny)
    Du_S_uAH    = zeros(T,1,Nx,Ny)
    Du_Fx_uAH   = zeros(T,1,Nx,Ny)
    Du_Fy_uAH   = zeros(T,1,Nx,Ny)
    Du_Sd_uAH   = zeros(T,1,Nx,Ny)
    Du_B1d_uAH  = zeros(T,1,Nx,Ny)
    Du_B2d_uAH  = zeros(T,1,Nx,Ny)
    Du_Gd_uAH   = zeros(T,1,Nx,Ny)
    Du_A_uAH    = zeros(T,1,Nx,Ny)

    Duu_B1_uAH  = zeros(T,1,Nx,Ny)
    Duu_B2_uAH  = zeros(T,1,Nx,Ny)
    Duu_G_uAH   = zeros(T,1,Nx,Ny)
    Duu_S_uAH   = zeros(T,1,Nx,Ny)
    Duu_Fx_uAH  = zeros(T,1,Nx,Ny)
    Duu_Fy_uAH  = zeros(T,1,Nx,Ny)
    Duu_A_uAH   = zeros(T,1,Nx,Ny)


    M = Nx * Ny
    axx     = zeros(T,M)
    ayy     = zeros(T,M)
    axy     = zeros(T,M)
    bx      = zeros(T,M)
    by      = zeros(T,M)
    cc      = zeros(T,M)
    b_vec   = zeros(T,M)


    # interpolate bulk functions (and u-derivatives) to the u=uAH surface
    @inbounds Threads.@threads for j in 1:Ny
        @inbounds for i in 1:Nx
            B1_uAH[1,i,j]       = interp(view(bulk.B1,  :,i,j))(uAH)
            B2_uAH[1,i,j]       = interp(view(bulk.B2,  :,i,j))(uAH)
            G_uAH[1,i,j]        = interp(view(bulk.G,   :,i,j))(uAH)
            phi_uAH[1,i,j]      = interp(view(bulk.phi, :,i,j))(uAH)
            S_uAH[1,i,j]        = interp(view(bulk.S,   :,i,j))(uAH)
            Fx_uAH[1,i,j]       = interp(view(bulk.Fx,  :,i,j))(uAH)
            Fy_uAH[1,i,j]       = interp(view(bulk.Fy,  :,i,j))(uAH)
            Sd_uAH[1,i,j]       = interp(view(bulk.Sd,  :,i,j))(uAH)
            B1d_uAH[1,i,j]      = interp(view(bulk.B1d, :,i,j))(uAH)
            B2d_uAH[1,i,j]      = interp(view(bulk.B2d, :,i,j))(uAH)
            Gd_uAH[1,i,j]       = interp(view(bulk.Gd,  :,i,j))(uAH)
            phid_uAH[1,i,j]     = interp(view(bulk.phid,:,i,j))(uAH)
            A_uAH[1,i,j]        = interp(view(bulk.A,   :,i,j))(uAH)

            Du_B1_uAH[1,i,j]    = interp(view(deriv.Du_B1,  :,i,j))(uAH)
            Du_B2_uAH[1,i,j]    = interp(view(deriv.Du_B2,  :,i,j))(uAH)
            Du_G_uAH[1,i,j]     = interp(view(deriv.Du_G,   :,i,j))(uAH)
            Du_phi_uAH[1,i,j]   = interp(view(deriv.Du_phi, :,i,j))(uAH)
            Du_S_uAH[1,i,j]     = interp(view(deriv.Du_S,   :,i,j))(uAH)
            Du_Fx_uAH[1,i,j]    = interp(view(deriv.Du_Fx,  :,i,j))(uAH)
            Du_Fy_uAH[1,i,j]    = interp(view(deriv.Du_Fy,  :,i,j))(uAH)
            Du_Sd_uAH[1,i,j]    = interp(view(deriv.Du_Sd,  :,i,j))(uAH)
            Du_B1d_uAH[1,i,j]   = interp(view(deriv.Du_B1d, :,i,j))(uAH)
            Du_B2d_uAH[1,i,j]   = interp(view(deriv.Du_B2d, :,i,j))(uAH)
            Du_Gd_uAH[1,i,j]    = interp(view(deriv.Du_Gd,  :,i,j))(uAH)
            Du_A_uAH[1,i,j]     = interp(view(deriv.Du_A,   :,i,j))(uAH)

            Duu_B1_uAH[1,i,j]   = interp(view(deriv.Duu_B1,  :,i,j))(uAH)
            Duu_B2_uAH[1,i,j]   = interp(view(deriv.Duu_B2,  :,i,j))(uAH)
            Duu_G_uAH[1,i,j]    = interp(view(deriv.Duu_G,   :,i,j))(uAH)
            Duu_S_uAH[1,i,j]    = interp(view(deriv.Duu_S,   :,i,j))(uAH)
            Duu_Fx_uAH[1,i,j]   = interp(view(deriv.Duu_Fx,  :,i,j))(uAH)
            Duu_Fy_uAH[1,i,j]   = interp(view(deriv.Duu_Fy,  :,i,j))(uAH)
            Duu_A_uAH[1,i,j]    = interp(view(deriv.Duu_A,   :,i,j))(uAH)
        end
    end


    ind2D   = LinearIndices(B1_uAH)

    # coefficients of the derivative operators
    @inbounds Threads.@threads for j in 1:Ny
        @inbounds for i in 1:Nx
            idx   = ind2D[i,j]

            xi    = gauge.xi[1,i,j]
            xi_x  = Dx(gauge.xi, 1,i,j)
            xi_y  = Dy(gauge.xi, 1,i,j)
            xi_xx = Dxx(gauge.xi, 1,i,j)
            xi_yy = Dyy(gauge.xi, 1,i,j)
            xi_xy = Dx(Dy, gauge.xi, 1,i,j)

            B1    = B1_uAH[1,i,j]
            B2    = B2_uAH[1,i,j]
            G     = G_uAH[1,i,j]
            phi   = phi_uAH[1,i,j]
            S     = S_uAH[1,i,j]
            Fx    = Fx_uAH[1,i,j]
            Fy    = Fy_uAH[1,i,j]
            Sd    = Sd_uAH[1,i,j]
            B1d   = B1d_uAH[1,i,j]
            B2d   = B2d_uAH[1,i,j]
            Gd    = Gd_uAH[1,i,j]
            phid  = phid_uAH[1,i,j]
            A     = A_uAH[1,i,j]

            # r derivatives

            B1p        = -u2 * Du_B1_uAH[1,i,j]
            B2p        = -u2 * Du_B2_uAH[1,i,j]
            Gp         = -u2 * Du_G_uAH[1,i,j]
            phip       = -u2 * Du_phi_uAH[1,i,j]
            Sp         = -u2 * Du_S_uAH[1,i,j]
            Fxp        = -u2 * Du_Fx_uAH[1,i,j]
            Fyp        = -u2 * Du_Fy_uAH[1,i,j]
            Sdp        = -u2 * Du_Sd_uAH[1,i,j]
            B1dp       = -u2 * Du_B1d_uAH[1,i,j]
            B2dp       = -u2 * Du_B2d_uAH[1,i,j]
            Gdp        = -u2 * Du_Gd_uAH[1,i,j]
            Ap         = -u2 * Du_A_uAH[1,i,j]

            B1pp       = 2*u3 * Du_B1_uAH[1,i,j]  + u4 * Duu_B1_uAH[1,i,j]
            B2pp       = 2*u3 * Du_B2_uAH[1,i,j]  + u4 * Duu_B2_uAH[1,i,j]
            Gpp        = 2*u3 * Du_G_uAH[1,i,j]   + u4 * Duu_G_uAH[1,i,j]
            Spp        = 2*u3 * Du_S_uAH[1,i,j]   + u4 * Duu_S_uAH[1,i,j]
            Fxpp       = 2*u3 * Du_Fx_uAH[1,i,j]  + u4 * Duu_Fx_uAH[1,i,j]
            Fypp       = 2*u3 * Du_Fy_uAH[1,i,j]  + u4 * Duu_Fy_uAH[1,i,j]
            App        = 2*u3 * Du_A_uAH[1,i,j]   + u4 * Duu_A_uAH[1,i,j]

            # x and y derivatives

            B1_x    = Dx(B1_uAH,   1,i,j)
            B2_x    = Dx(B2_uAH,   1,i,j)
            G_x     = Dx(G_uAH,    1,i,j)
            phi_x   = Dx(phi_uAH,  1,i,j)
            S_x     = Dx(S_uAH,    1,i,j)
            Fx_x    = Dx(Fx_uAH,   1,i,j)
            Fy_x    = Dx(Fy_uAH,   1,i,j)
            Sd_x    = Dx(Sd_uAH,   1,i,j)
            B1d_x   = Dx(B1d_uAH,  1,i,j)
            B2d_x   = Dx(B2d_uAH,  1,i,j)
            Gd_x    = Dx(Gd_uAH,   1,i,j)
            A_x     = Dx(A_uAH,    1,i,j)

            B1_y    = Dy(B1_uAH,   1,i,j)
            B2_y    = Dy(B2_uAH,   1,i,j)
            G_y     = Dy(G_uAH,    1,i,j)
            phi_y   = Dy(phi_uAH,  1,i,j)
            S_y     = Dy(S_uAH,    1,i,j)
            Fx_y    = Dy(Fx_uAH,   1,i,j)
            Fy_y    = Dy(Fy_uAH,   1,i,j)
            Sd_y    = Dy(Sd_uAH,   1,i,j)
            B1d_y   = Dy(B1d_uAH,  1,i,j)
            B2d_y   = Dy(B2d_uAH,  1,i,j)
            Gd_y    = Dy(Gd_uAH,   1,i,j)
            A_y     = Dy(A_uAH,    1,i,j)

            B1p_x   = -u2 * Dx(Du_B1_uAH, 1,i,j)
            B2p_x   = -u2 * Dx(Du_B2_uAH, 1,i,j)
            Gp_x    = -u2 * Dx(Du_G_uAH,  1,i,j)
            Sp_x    = -u2 * Dx(Du_S_uAH,  1,i,j)
            Fxp_x   = -u2 * Dx(Du_Fx_uAH, 1,i,j)
            Fyp_x   = -u2 * Dx(Du_Fy_uAH, 1,i,j)
            Ap_x    = -u2 * Dx(Du_A_uAH,  1,i,j)

            B1p_y   = -u2 * Dy(Du_B1_uAH, 1,i,j)
            B2p_y   = -u2 * Dy(Du_B2_uAH, 1,i,j)
            Gp_y    = -u2 * Dy(Du_G_uAH,  1,i,j)
            Sp_y    = -u2 * Dy(Du_S_uAH,  1,i,j)
            Fxp_y   = -u2 * Dy(Du_Fx_uAH, 1,i,j)
            Fyp_y   = -u2 * Dy(Du_Fy_uAH, 1,i,j)
            Ap_y    = -u2 * Dy(Du_A_uAH,  1,i,j)

            Fy_xx   = Dxx(Fy_uAH,  1,i,j)
            A_xx    = Dxx(A_uAH,   1,i,j)

            Fx_yy   = Dyy(Fx_uAH,  1,i,j)
            A_yy    = Dyy(A_uAH,   1,i,j)

            Fx_xy   = Dx(Dy, Fx_uAH, 1,i,j)
            Fy_xy   = Dx(Dy, Fy_uAH, 1,i,j)
            A_xy    = Dx(Dy, A_uAH,  1,i,j)

            vars =  (
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
                                                  Fx_xy , Fy_xy ,                                A_xy
            )

            a11, a22, a12, b1, b2, c, S = xi_t_eq_coeff(vars, sys.gridtype)

            axx[idx]   = a11
            ayy[idx]   = a22
            axy[idx]   = a12
            bx[idx]    = b1
            by[idx]    = b2
            cc[idx]    = c

            b_vec[idx] = -S
        end
    end

    nothing
end

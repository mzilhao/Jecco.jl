
function compute_xi_t!(gauge_t::Gauge, bulkconstrain::BulkConstrained, bulkevol::BulkEvolved,
                       deriv::BulkDeriv, gauge::Gauge, sys::System{Outer}, ::ConstantGauge)
    xi_t = getxi(gauge_t)
    fill!(xi_t, 0)
    nothing
end

struct TestGauge <: GaugeCondition end


# TODO
function compute_xi_t!(gauge_t::Gauge, bulkconstrain::BulkConstrained, bulkevol::BulkEvolved,
                       deriv::BulkDeriv, gauge::Gauge, sys::System{Outer}, ::TestGauge)
    bulk = Bulk(bulkevol, bulkconstrain)
    _, Nx, Ny = size(sys)

    Dx  = sys.Dx
    Dxx = sys.Dxx
    Dy  = sys.Dy
    Dyy = sys.Dyy

    interp = sys.uinterp

    xi_t = getxi(gauge_t)

    T = Float64
    uAH = 1.0

    M = Nx * Ny

    B1_uAH   = zeros(T,1,Nx,Ny)
    B2_uAH   = zeros(T,1,Nx,Ny)
    G_uAH    = zeros(T,1,Nx,Ny)
    S_uAH    = zeros(T,1,Nx,Ny)
    Fx_uAH   = zeros(T,1,Nx,Ny)
    Fy_uAH   = zeros(T,1,Nx,Ny)
    Sd_uAH   = zeros(T,1,Nx,Ny)
    B1d_uAH  = zeros(T,1,Nx,Ny)
    B2d_uAH  = zeros(T,1,Nx,Ny)
    Gd_uAH   = zeros(T,1,Nx,Ny)
    phid_uAH = zeros(T,1,Nx,Ny)
    A_uAH    = zeros(T,1,Nx,Ny)


    axx     = zeros(T,M)
    ayy     = zeros(T,M)
    axy     = zeros(T,M)
    bx      = zeros(T,M)
    by      = zeros(T,M)
    cc      = zeros(T,M)
    b_vec   = zeros(T,M)


    # intepolate all bulk functions to the u=uAH surface

    @inbounds Threads.@threads for j in 1:Ny
        @inbounds for i in 1:Nx
            B1_uAH[1,i,j]   = interp(view(bulk.B1,  :,i,j))(uAH)
            B2_uAH[1,i,j]   = interp(view(bulk.B2,  :,i,j))(uAH)
            G_uAH[1,i,j]    = interp(view(bulk.G,   :,i,j))(uAH)
            S_uAH[1,i,j]    = interp(view(bulk.S,   :,i,j))(uAH)
            Fx_uAH[1,i,j]   = interp(view(bulk.Fx,  :,i,j))(uAH)
            Fy_uAH[1,i,j]   = interp(view(bulk.Fy,  :,i,j))(uAH)
            Sd_uAH[1,i,j]   = interp(view(bulk.Sd,  :,i,j))(uAH)
            B1d_uAH[1,i,j]  = interp(view(bulk.B1d, :,i,j))(uAH)
            B2d_uAH[1,i,j]  = interp(view(bulk.B2d, :,i,j))(uAH)
            Gd_uAH[1,i,j]   = interp(view(bulk.Gd,  :,i,j))(uAH)
            phid_uAH[1,i,j] = interp(view(bulk.phid,:,i,j))(uAH)
            A_uAH[1,i,j]    = interp(view(bulk.A,   :,i,j))(uAH)

            # TODO: need the u derivatives here as well...

        end
    end




    ind2D   = LinearIndices(B1_uAH)

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
            S     = S_uAH[1,i,j]
            Fx    = Fx_uAH[1,i,j]
            Fy    = Fy_uAH[1,i,j]
            Sd    = Sd_uAH[1,i,j]
            B1d   = B1d_uAH[1,i,j]
            B2d   = B2d_uAH[1,i,j]
            Gd    = Gd_uAH[1,i,j]
            phid  = phid_uAH[1,i,j]
            A     = A_uAH[1,i,j]

            # TODO: r derivatives


            # x and y derivatives

            B1_x    = Dx(B1_uAH,   1,i,j)
            B2_x    = Dx(B2_uAH,   1,i,j)
            G_x     = Dx(G_uAH,    1,i,j)
            S_x     = Dx(S_uAH,    1,i,j)
            Fx_x    = Dx(Fx_uAH,   1,i,j)
            Fy_x    = Dx(Fy_uAH,   1,i,j)
            Sd_x    = Dx(Sd_uAH,   1,i,j)
            B1d_x   = Dx(B1d_uAH,  1,i,j)
            B2d_x   = Dx(B2d_uAH,  1,i,j)
            Gd_x    = Dx(Gd_uAH,   1,i,j)
            phid_x  = Dx(phid_uAH, 1,i,j)
            A_x     = Dx(A_uAH,    1,i,j)


            # xi_t_eq_coeff

            # FIXME
            axx[idx] = 0
            ayy[idx] = 0
            axy[idx] = 0
            bx[idx]  = 0
            by[idx]  = 0
            cc[idx]  = 0

            b_vec[idx] = 0
        end
    end



    nothing
end


# TODO
function compute_xi_t!(gauge_t::Gauge, bulkconstrain::BulkConstrained, bulkevol::BulkEvolved,
                       deriv::BulkDeriv, gauge::Gauge, sys::System{Outer}, gaugecondition::ConstantAH)
    _, Nx, Ny = size(sys)

    xi_t = getxi(gauge_t)

    fill!(xi_t, 0)

    nothing
end

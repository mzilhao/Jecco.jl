
function compute_bulkevolved_t!(bulkevol_t::BulkEvolved,
                                bulkconstrain::BulkConstrained, gauge_t::Gauge,
                                bulkevol::BulkEvolved, boundary::Boundary,
                                gauge::Gauge, sys::System, ::EvolTest0)

    B1_t, B2_t, G_t, phi_t = unpack(bulkevol_t)
    # B1  , B2  , G  , phi   = unpack(bulkevol)

    fill!(B1_t,  0)
    fill!(B2_t,  0)
    fill!(G_t,   0)
    fill!(phi_t, 0)

    nothing
end

# TODO
function compute_bulkevolved_t!(bulkevol_t::BulkEvolved,
                                bulkconstrain::BulkConstrained, gauge_t::Gauge,
                                bulkevol::BulkEvolved, boundary::Boundary,
                                gauge::Gauge, sys::System{Inner}, evoleq::AffineNull)
    Du  = sys.Du
    Dx  = sys.Dx
    Dy  = sys.Dy

    phi0  = evoleq.phi0
    phi03 = phi0 * phi0 * phi0

    Nu, Nx, Ny = size(sys)

    B1_t, B2_t, G_t, phi_t = unpack(bulkevol_t)

    @fastmath @inbounds for j in 1:Ny
        @inbounds @simd for i in 1:Nx
            xi     = gauge.xi[1,i,j]

            B1     = bulkevol.B1[1,i,j]
            B2     = bulkevol.B2[1,i,j]
            G      = bulkevol.G[1,i,j]
            phi    = bulkevol.phi[1,i,j]

            B1_u   = Du(bulkevol.B1, 1,i,j)
            B2_u   = Du(bulkevol.B2, 1,i,j)
            G_u    = Du(bulkevol.G,  1,i,j)
            phi_u  = Du(bulkevol.phi,1,i,j)

            B1d_u  = Du(bulkconstrain.B1d, 1,i,j)
            B2d_u  = Du(bulkconstrain.B2d, 1,i,j)
            Gd_u   = Du(bulkconstrain.Gd,  1,i,j)
            phid_u = Du(bulkconstrain.phid,1,i,j)

            # u = 0. CHECK!!
            B1_t[1,i,j]  = B1d_u  + 2.5 * B1_u  + 4 * B1  * xi
            B2_t[1,i,j]  = B2d_u  + 2.5 * B2_u  + 4 * B2  * xi
            G_t[1,i,j]   = Gd_u   + 2.5 * G_u   + 4 * G   * xi

            # TODO
            # phi_t[1,i,j] = 0
        end
    end


    # TODO
    @fastmath @inbounds for j in 1:Ny
        @inbounds for i in 1:Nx
            @inbounds @simd for a in 2:Nu
                xi     = gauge.xi[1,i,j]

                B1     = bulkevol.B1[a,i,j]
                B2     = bulkevol.B2[a,i,j]
                G      = bulkevol.G[a,i,j]
                phi    = bulkevol.phi[a,i,j]

                B1_u   = Du(bulkevol.B1, a,i,j)
                B2_u   = Du(bulkevol.B2, a,i,j)
                G_u    = Du(bulkevol.G,  a,i,j)
                phi_u  = Du(bulkevol.phi,a,i,j)

                B1d_u  = Du(bulkconstrain.B1d, a,i,j)
                B2d_u  = Du(bulkconstrain.B2d, a,i,j)
                Gd_u   = Du(bulkconstrain.Gd,  a,i,j)
                phid_u = Du(bulkconstrain.phid,a,i,j)

            end
        end
    end


    fill!(B1_t,  0)
    fill!(B2_t,  0)
    fill!(G_t,   0)
    fill!(phi_t, 0)

    nothing
end


# TODO
function compute_bulkevolved_t!(bulkevol_t::BulkEvolved,
                                bulkconstrain::BulkConstrained, gauge_t::Gauge,
                                bulkevol::BulkEvolved, boundary::Boundary,
                                gauge::Gauge, sys::System{Outer}, evoleq::AffineNull)
    Du  = sys.Du
    Dx  = sys.Dx
    Dy  = sys.Dy

    phi0  = evoleq.phi0
    phi03 = phi0 * phi0 * phi0

    Nu, Nx, Ny = size(sys)


    B1_t, B2_t, G_t, phi_t = unpack(bulkevol_t)
    # B1  , B2  , G  , phi   = unpack(bulkevol)

    fill!(B1_t,  0)
    fill!(B2_t,  0)
    fill!(G_t,   0)
    fill!(phi_t, 0)

    nothing
end

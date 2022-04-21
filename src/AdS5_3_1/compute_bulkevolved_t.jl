
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


function compute_bulkevolved_t!(bulkevol_t::BulkEvolved,
                                bulkconstrain::BulkConstrained, gauge_t::Gauge,
                                bulkevol::BulkEvolved, boundary::Boundary,
                                gauge::Gauge, sys::System{Inner}, evoleq::AffineNull)
    uu  = sys.ucoord
    Du  = sys.Du
    Dx  = sys.Dx
    Dy  = sys.Dy

    phi0  = evoleq.phi0
    phi02 = phi0 * phi0
    phi03 = phi0 * phi02

    Nu, Nx, Ny = size(sys)

    B1_t, B2_t, G_t, phi_t = unpack(bulkevol_t)

    # u = 0
    @fastmath @inbounds @threads for j in 1:Ny
        @inbounds @simd for i in 1:Nx
            xi     = gauge.xi[1,i,j]

            B1     = bulkevol.B1[1,i,j]
            B2     = bulkevol.B2[1,i,j]
            G      = bulkevol.G[1,i,j]

            B1_u   = Du(bulkevol.B1, 1,i,j)
            B2_u   = Du(bulkevol.B2, 1,i,j)
            G_u    = Du(bulkevol.G,  1,i,j)

            B1d_u  = Du(bulkconstrain.B1d, 1,i,j)
            B2d_u  = Du(bulkconstrain.B2d, 1,i,j)
            Gd_u   = Du(bulkconstrain.Gd,  1,i,j)

            B1_t[1,i,j]  = B1d_u  + 2.5 * B1_u  + 4 * B1  * xi
            B2_t[1,i,j]  = B2d_u  + 2.5 * B2_u  + 4 * B2  * xi
            G_t[1,i,j]   = Gd_u   + 2.5 * G_u   + 4 * G   * xi
        end
    end

    # remaining inner grid points
    @fastmath @inbounds @threads for j in 1:Ny
        @inbounds for i in 1:Nx
            xi    = gauge.xi[1,i,j]
            xi_t  = gauge_t.xi[1,i,j]
            @inbounds @simd for a in 2:Nu
                u      = uu[a]
                u4     = u * u * u * u

                B1     = bulkevol.B1[a,i,j]
                B2     = bulkevol.B2[a,i,j]
                G      = bulkevol.G[a,i,j]

                B1d    = bulkconstrain.B1d[a,i,j]
                B2d    = bulkconstrain.B2d[a,i,j]
                Gd     = bulkconstrain.Gd[a,i,j]
                A      = bulkconstrain.A[a,i,j]

                B1_u   = Du(bulkevol.B1, a,i,j)
                B2_u   = Du(bulkevol.B2, a,i,j)
                G_u    = Du(bulkevol.G,  a,i,j)

		B1_t[a,i,j] = ((u * B1_u + 4 * B1) *
                               (-2 * u * u * xi_t + A * u4 +
                                (xi * u + 1) * (xi * u + 1)) +
                               2 * B1d) / (2 * u) -
                               phi02 * u * (4 * B1 + B1_u * u) / 3

		B2_t[a,i,j] = ((u * B2_u + 4 * B2) *
                               (-2 * u * u * xi_t + A * u4 +
                                (xi * u + 1) * (xi * u + 1)) +
                               2 * B2d) / (2 * u) -
                               phi02 * u * (4 * B2 + B2_u * u) / 3

		G_t[a,i,j] = ((u * G_u + 4 * G) *
                               (-2 * u * u * xi_t + A * u4 +
                                (xi * u + 1) * (xi * u + 1)) +
                               2 * Gd) / (2 * u) -
                               phi02 * u * (4 * G + G_u * u) / 3
            end
        end
    end

    # if phi0 = 0 set phi_t to zero and return
    if abs(phi0) < 1e-9
        fill!(phi_t, 0)
        return
    end

    # otherwise, compute phi_t

    # u = 0
    @fastmath @inbounds @threads for j in 1:Ny
        @inbounds @simd for i in 1:Nx
            xi     = gauge.xi[1,i,j]
            xi3    = xi * xi * xi
            xi_t   = gauge_t.xi[1,i,j]

            phi    = bulkevol.phi[1,i,j]

            phi_u  = Du(bulkevol.phi,1,i,j)

            phid_u = Du(bulkconstrain.phid,1,i,j)

            phi_t[1,i,j] = -xi3 / phi02 + 3 * xi * phi +  2//3 * xi +
                2 * xi * xi_t / phi02 + phid_u + 2 * phi_u
        end
    end

    # remaining inner grid points
    @fastmath @inbounds @threads for j in 1:Ny
        @inbounds for i in 1:Nx
            xi     = gauge.xi[1,i,j]
            xi2    = xi * xi
            xi3    = xi * xi * xi
            xi_t   = gauge_t.xi[1,i,j]
            @inbounds @simd for a in 2:Nu
                u      = uu[a]
                u2     = u * u
                u4     = u2 * u2

                phi    = bulkevol.phi[a,i,j]
                phid   = bulkconstrain.phid[a,i,j]
                A      = bulkconstrain.A[a,i,j]

                phi_u  = Du(bulkevol.phi,a,i,j)

                phi_t[a,i,j] = (
                    - 6 * xi3 * u
                    + 3 * xi2 * (-3 + phi02 * u2 * (3 * phi + phi_u * u))
		    + 3 * A * u2 * (1 - 2 * xi * u
				    + phi02 * u2 * (3 * phi + phi_u * u))
		    + 2 * xi * u * (6 * xi_t + 9 * phi * phi02
				    + phi02 * (2 + 3 * phi_u * u))
		    + phi02 * (-2 + 6 * phid - (3 * phi + phi_u * u) *
                               (-3 + 6 * xi_t * u2 + 2 * phi02 * u2))
                ) / (6 * phi02 * u)
            end
        end
    end

    nothing
end


function compute_bulkevolved_t!(bulkevol_t::BulkEvolved,
                                bulkconstrain::BulkConstrained, gauge_t::Gauge,
                                bulkevol::BulkEvolved, boundary::Boundary,
                                gauge::Gauge, sys::System{Outer}, evoleq::AffineNull)
    uu  = sys.ucoord
    Du  = sys.Du
    Dx  = sys.Dx
    Dy  = sys.Dy

    phi0  = evoleq.phi0
    phi02 = phi0 * phi0
    phi03 = phi0 * phi02

    Nu, Nx, Ny = size(sys)

    B1_t, B2_t, G_t, phi_t = unpack(bulkevol_t)

    @fastmath @inbounds @threads for j in 1:Ny
        @inbounds for i in 1:Nx
            xi_t  = gauge_t.xi[1,i,j]
            @inbounds @simd for a in 1:Nu
                u      = uu[a]
                u2     = u * u

                B1d    = bulkconstrain.B1d[a,i,j]
                B2d    = bulkconstrain.B2d[a,i,j]
                Gd     = bulkconstrain.Gd[a,i,j]
                A      = bulkconstrain.A[a,i,j]

                B1_u   = Du(bulkevol.B1, a,i,j)
                B2_u   = Du(bulkevol.B2, a,i,j)
                G_u    = Du(bulkevol.G,  a,i,j)

		B1_t[a,i,j] = B1d + u2 * (A/2 - xi_t) * B1_u
		B2_t[a,i,j] = B2d + u2 * (A/2 - xi_t) * B2_u
		G_t[a,i,j]  = Gd  + u2 * (A/2 - xi_t) * G_u
            end
        end
    end

    # if phi0 = 0 set phi_t to zero and return
    if abs(phi0) < 1e-9
        fill!(phi_t, 0)
        return
    end

    # otherwise, compute phi_t

    @fastmath @inbounds @threads for j in 1:Ny
        @inbounds for i in 1:Nx
            xi_t   = gauge_t.xi[1,i,j]
            @inbounds @simd for a in 1:Nu
                u      = uu[a]
                u2     = u * u

                phid   = bulkconstrain.phid[a,i,j]
                A      = bulkconstrain.A[a,i,j]

                phi_u  = Du(bulkevol.phi, a,i,j)

		phi_t[a,i,j] = phid + u2 * (A/2 - xi_t) * phi_u
            end
        end
    end

    nothing
end

function sync_bulkevolved!(bulkevols_t, bulkconstrains, gauge_t::Gauge,
                           systems::SystemPartition, evoleq::AffineNull)
    Nsys = length(systems)

    @inbounds for i in 1:Nsys-1
        sync_bulkevolved!(bulkevols_t[i], bulkevols_t[i+1], bulkconstrains[i+1],
                          gauge_t, systems[i], systems[i+1], evoleq)
    end

    nothing
end

function sync_bulkevolved!(bulkevol1_t::BulkEvolved, bulkevol2_t::BulkEvolved,
                           bulkconstrain2::BulkConstrained, gauge_t::Gauge,
                           sys1::System{Outer}, sys2::System{Outer}, evoleq::AffineNull)
    u0 = sys2.ucoord[1]

    _, Nx, Ny = size(sys2)

    B11_t, B21_t, G1_t, phi1_t = unpack(bulkevol1_t)
    B12_t, B22_t, G2_t, phi2_t = unpack(bulkevol2_t)

    @fastmath @inbounds @threads for j in 1:Ny
        @inbounds for i in 1:Nx
            xi_t   = gauge_t.xi[1,i,j]
            A      = bulkconstrain2.A[1,i,j]

            # characteristic speed
            c = u0 * u0 * (A/2 - xi_t)

            # if c > 0, mode is entering grid1 from grid2;
            # if c < 0, mode is entering grid2 from grid1.
            # we assume here that grids merely touch at the interface
            if c > 0
                B11_t[end,i,j]  = B12_t[1,i,j]
                B21_t[end,i,j]  = B22_t[1,i,j]
                G1_t[end,i,j]   = G2_t[1,i,j]
                phi1_t[end,i,j] = phi2_t[1,i,j]
            elseif c < 0
                B12_t[1,i,j]  = B11_t[end,i,j]
                B22_t[1,i,j]  = B21_t[end,i,j]
                G2_t[1,i,j]   = G1_t[end,i,j]
                phi2_t[1,i,j] = phi1_t[end,i,j]
            end
        end
    end

    nothing
end

function sync_bulkevolved!(bulkevol1_t::BulkEvolved, bulkevol2_t::BulkEvolved,
                           bulkconstrain2::BulkConstrained, gauge_t::Gauge,
                           sys1::System{Inner}, sys2::System{Outer}, evoleq::AffineNull)
    u0  = sys2.ucoord[1]
    u02 = u0 * u0
    u03 = u0 * u02
    u04 = u02 * u02

    phi0  = evoleq.phi0
    phi02 = phi0 * phi0
    phi03 = phi0 * phi02

    _, Nx, Ny = size(sys2)

    B11_t, B21_t, G1_t, phi1_t = unpack(bulkevol1_t)
    B12_t, B22_t, G2_t, phi2_t = unpack(bulkevol2_t)

    @fastmath @inbounds @threads for j in 1:Ny
        @inbounds for i in 1:Nx
            xi_t   = gauge_t.xi[1,i,j]
            A      = bulkconstrain2.A[1,i,j]

            # characteristic speed
            c = u0 * u0 * (A/2 - xi_t)

            # if c > 0, mode is entering grid1 from grid2;
            # if c < 0, mode is entering grid2 from grid1.
            # we assume here that grids merely touch at the interface
            if c > 0
                B11_t[end,i,j]  = B12_t[1,i,j] / u04
                B21_t[end,i,j]  = B22_t[1,i,j] / u04
                G1_t[end,i,j]   = G2_t[1,i,j] / u04
            elseif c < 0
                B12_t[1,i,j]  = u04 * B11_t[end,i,j]
                B22_t[1,i,j]  = u04 * B21_t[end,i,j]
                G2_t[1,i,j]   = u04 * G1_t[end,i,j]
            end
        end
    end

    # if phi0 = 0 return
    if abs(phi0) < 1e-9
        return
    end

    @fastmath @inbounds @threads for j in 1:Ny
        @inbounds for i in 1:Nx
            xi_t   = gauge_t.xi[1,i,j]
            A      = bulkconstrain2.A[1,i,j]

            # characteristic speed
            c = u0 * u0 * (A/2 - xi_t)

            # if c > 0, mode is entering grid1 from grid2;
            # if c < 0, mode is entering grid2 from grid1.
            # we assume here that grids merely touch at the interface
            if c > 0
                phi1_t[end,i,j] = phi2_t[1,i,j] / (u03 * phi03) + xi_t / (u0 * phi02)
            elseif c < 0
                phi2_t[1,i,j] = -phi0 * u02 * xi_t + u03 * phi03 * phi1_t[end,i,j]
            end
        end
    end

    nothing
end

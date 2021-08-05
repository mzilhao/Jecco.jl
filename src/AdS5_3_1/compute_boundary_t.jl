
function compute_boundary_t!(boundary_t::Boundary, bulkevol::BulkEvolved,
                             boundary::Boundary, gauge::Gauge, sys::System, ::EvolTest0)

    a4_t, fx2_t, fy2_t = unpack(boundary_t)
    # a4  , fx2  , fy2   = unpack(boundary)

    fill!(a4_t,  0)
    fill!(fx2_t, 0)
    fill!(fy2_t, 0)

    nothing
end

function compute_boundary_t!(boundary_t::Boundary, bulk::BulkEvolved,
                             boundary::Boundary, gauge::Gauge, sys::System{Inner},
                             evoleq::AffineNull)
    Du  = sys.Du
    Dx  = sys.Dx
    Dy  = sys.Dy

    phi0  = evoleq.phi0
    phi03 = phi0 * phi0 * phi0

    _, Nx, Ny = size(sys)

    a4_t, fx2_t, fy2_t = unpack(boundary_t)

    @fastmath @inbounds for j in 1:Ny
        @inbounds for i in 1:Nx
            xi      = gauge.xi[1,i,j]
            xi3     = xi*xi*xi

            phi     = bulk.phi[1,i,j]

            xi_x    = Dx(gauge.xi, 1,i,j)
            xi_y    = Dy(gauge.xi, 1,i,j)

            phi_u   = Du(bulk.phi, 1,i,j)
            phi_x   = Dx(bulk.phi, 1,i,j)
            phi_y   = Dy(bulk.phi, 1,i,j)

            b14_x   = Dx(bulk.B1, 1,i,j)
            b14_y   = Dy(bulk.B1, 1,i,j)

            b24_x   = Dx(bulk.B2, 1,i,j)
            b24_y   = Dy(bulk.B2, 1,i,j)

            g4_x    = Dx(bulk.G, 1,i,j)
            g4_y    = Dy(bulk.G, 1,i,j)

            a4_x    = Dx(boundary.a4, 1,i,j)
            a4_y    = Dy(boundary.a4, 1,i,j)

            fx2_x   = Dx(boundary.fx2, 1,i,j)
            fy2_y   = Dy(boundary.fy2, 1,i,j)

            # phi2 = phi0^3 phi(u=0) - phi0 xi^2
            phi2    = phi03 * phi - phi0 * xi * xi
            phi2_x  = phi03 * phi_x - 2 * phi0 * xi * xi_x
            phi2_y  = phi03 * phi_y - 2 * phi0 * xi * xi_y

            # phi2_t = 3 xi phi2 + phi0 xi^3 + phi0^3 phi(u=0)_u
            phi2_t  = 3 * xi * phi2 + phi0 * xi3 + phi03 * phi_u

            a4_t[1,i,j]  = -4//3 * (fx2_x + fy2_y + phi0 * phi2_t)

            fx2_t[1,i,j] = g4_y - b14_x - b24_x - 0.25 * a4_x + 1//3 * phi0 * phi2_x
            fy2_t[1,i,j] = g4_x + b14_y - b24_y - 0.25 * a4_y + 1//3 * phi0 * phi2_y
        end
    end

    nothing
end

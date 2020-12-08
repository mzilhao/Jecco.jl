
function compute_boundary_t!(boundary_t::Boundary, bulkevol::BulkEvolved,
                             boundary::Boundary, gauge::Gauge, sys::System, ::EvolTest0)

    a4_t  = geta4(boundary_t)
    fx2_t = getfx2(boundary_t)
    fy2_t = getfy2(boundary_t)

    fill!(a4_t,  0)
    fill!(fx2_t, 0)
    fill!(fy2_t, 0)

    nothing
end

function compute_boundary_t!(boundary_t::Boundary, bulk::BulkEvolved,
                             boundary::Boundary, gauge::Gauge, sys::System{Inner},
                             evoleq::AffineNull)
    B1GF    = getB1(bulk)
    B2GF    = getB2(bulk)
    GGF     = getG(bulk)
    phiGF   = getphi(bulk)
    a4GF    = geta4(boundary)
    fx2GF   = getfx2(boundary)
    fy2GF   = getfy2(boundary)
    xiGF    = getxi(gauge)

    Du  = sys.Du
    Dx  = sys.Dx
    Dy  = sys.Dy

    phi0  = evoleq.phi0
    phi03 = phi0 * phi0 * phi0

    _, Nx, Ny = size(sys)

    a4_t  = geta4(boundary_t)
    fx2_t = getfx2(boundary_t)
    fy2_t = getfy2(boundary_t)

    @fastmath @inbounds for j in 1:Ny
        @inbounds for i in 1:Nx
            xi      = xiGF[1,i,j]
            xi3     = xi*xi*xi

            phi     = phiGF[1,i,j]

            xi_x    = Dx(xiGF, 1,i,j)
            xi_y    = Dy(xiGF, 1,i,j)

            phi_u   = Du(phiGF, 1,i,j)
            phi_x   = Dx(phiGF, 1,i,j)
            phi_y   = Dy(phiGF, 1,i,j)

            b14_x   = Dx(B1GF, 1,i,j)
            b14_y   = Dy(B1GF, 1,i,j)

            b24_x   = Dx(B2GF, 1,i,j)
            b24_y   = Dy(B2GF, 1,i,j)

            g4_x    = Dx(GGF, 1,i,j)
            g4_y    = Dy(GGF, 1,i,j)

            a4_x    = Dx(a4GF, 1,i,j)
            a4_y    = Dy(a4GF, 1,i,j)

            fx2_x   = Dx(fx2GF, 1,i,j)
            fy2_y   = Dy(fy2GF, 1,i,j)

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

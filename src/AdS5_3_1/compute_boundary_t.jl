
function compute_boundary_t!(boundary_t::Boundary, bulkevol::BulkEvolved,
                             boundary::Boundary, gauge::Gauge, sys::System, ::EvolTest0)

    a4_t, fx2_t, fy2_t = unpack(boundary_t)
    # a4  , fx2  , fy2   = unpack(boundary)

    fill!(a4_t,  0)
    fill!(fx2_t, 0)
    fill!(fy2_t, 0)

    nothing
end

function compute_boundary_t!(boundary_t::Boundary, bulkevol::BulkEvolved,
                             boundary::Boundary, gauge::Gauge, sys::System{Inner},
                             ::AffineNull)
    Dx  = sys.Dx
    Dy  = sys.Dy

    _, Nx, Ny = size(sys)

    a4_t, fx2_t, fy2_t = unpack(boundary_t)
    a4  , fx2  , fy2   = unpack(boundary)

    xi = getxi(gauge)

    B1 = getB1(bulkevol)
    B2 = getB2(bulkevol)

    @fastmath @inbounds for j in 1:Ny
        @inbounds for i in 1:Nx
            b14_x = Dx(B1, 1,i,j)
            b14_y = Dy(B1, 1,i,j)

            # TODO
            a4_t[1,i,j]  = 0
            fx2_t[1,i,j] = 0
            fy2_t[1,i,j] = 0
        end
    end

    nothing
end

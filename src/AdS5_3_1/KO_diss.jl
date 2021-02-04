
function apply_dissipation_3D!(var_t, var, sys::System)
    Nu, Nx, Ny = size(sys)

    DKOx = sys.DKOx
    DKOy = sys.DKOy

    @fastmath @inbounds @threads for j in 1:Ny
        @inbounds for i in 1:Nx
            @inbounds for a in 1:Nu
                var_t[a,i,j] += DKOx(var, a,i,j) + DKOy(var, a,i,j)
            end
        end
    end

    nothing
end

function apply_dissipation_2D!(var_t, var, sys::System)
    _, Nx, Ny = size(sys)

    DKOx = sys.DKOx
    DKOy = sys.DKOy

    @fastmath @inbounds @threads for j in 1:Ny
        @inbounds for i in 1:Nx
            var_t[1,i,j] += DKOx(var, 1,i,j) + DKOy(var, 1,i,j)
        end
    end

    nothing
end


function apply_dissipation!(boundary_t::Boundary, boundary::Boundary, sys::System)
    a4  = geta4(boundary)
    fx2 = getfx2(boundary)
    fy2 = getfy2(boundary)

    a4_t  = geta4(boundary_t)
    fx2_t = getfx2(boundary_t)
    fy2_t = getfy2(boundary_t)

    @sync begin
        @spawn apply_dissipation_2D!( a4_t,  a4, sys)
        @spawn apply_dissipation_2D!(fx2_t, fx2, sys)
        @spawn apply_dissipation_2D!(fy2_t, fy2, sys)
    end

    nothing
end

function apply_dissipation!(gauge_t::Gauge, gauge::Gauge, sys::System)
    xi   = getxi(gauge)
    xi_t = getxi(gauge_t)

    apply_dissipation_2D!(xi_t, xi, sys)

    nothing
end

function apply_dissipation!(bulkevol_t::BulkEvolved, bulkevol::BulkEvolved,
                            sys::System)
    B1  = getB1(bulkevol)
    B2  = getB2(bulkevol)
    G   = getG(bulkevol)
    phi = getphi(bulkevol)

    B1_t  = getB1(bulkevol_t)
    B2_t  = getB2(bulkevol_t)
    G_t   = getG(bulkevol_t)
    phi_t = getphi(bulkevol_t)

    @sync begin
        @spawn apply_dissipation_3D!( B1_t,  B1, sys)
        @spawn apply_dissipation_3D!( B2_t,  B2, sys)
        @spawn apply_dissipation_3D!(  G_t,   G, sys)
        @spawn apply_dissipation_3D!(phi_t, phi, sys)
    end

    nothing
end

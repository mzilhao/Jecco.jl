
function apply_dissipation_3D!(var_out, var, sys::System)
    Nu, Nx, Ny = size(sys)

    DKOx = sys.DKOx
    DKOy = sys.DKOy

    copyto!(var_out, var)

    @fastmath @inbounds @threads for j in 1:Ny
        @inbounds for i in 1:Nx
            @inbounds for a in 1:Nu
                var_out[a,i,j] += DKOx(var, a,i,j) + DKOy(var, a,i,j)
            end
        end
    end

    nothing
end

function apply_dissipation_2D!(var_out, var, sys::System)
    _, Nx, Ny = size(sys)

    DKOx = sys.DKOx
    DKOy = sys.DKOy

    copyto!(var_out, var)

    @fastmath @inbounds @threads for j in 1:Ny
        @inbounds for i in 1:Nx
            var_out[1,i,j] += DKOx(var, 1,i,j) + DKOy(var, 1,i,j)
        end
    end

    nothing
end

function apply_dissipation!(boundary::Boundary, cache::Boundary, sys::System)
    a4  = geta4(boundary)
    fx2 = getfx2(boundary)
    fy2 = getfy2(boundary)

    a4_cache  = geta4(cache)
    fx2_cache = getfx2(cache)
    fy2_cache = getfy2(cache)

    # the loops in these functions are threaded, so it's probably not worth it
    # to @spawn here
    apply_dissipation_2D!( a4_cache,  a4, sys)
    apply_dissipation_2D!(fx2_cache, fx2, sys)
    apply_dissipation_2D!(fy2_cache, fy2, sys)

    copyto!( a4,  a4_cache)
    copyto!(fx2, fx2_cache)
    copyto!(fy2, fy2_cache)

    nothing
end

function apply_dissipation!(gauge::Gauge, cache::Gauge, sys::System)
    xi = getxi(gauge)
    xi_cache = getxi(cache)

    apply_dissipation_2D!(xi_cache, xi, sys)
    copyto!(xi, xi_cache)

    nothing
end

function apply_dissipation!(bulkevol::BulkEvolved, cache::BulkEvolved,
                            sys::System)
    B1  = getB1(bulkevol)
    B2  = getB2(bulkevol)
    G   = getG(bulkevol)
    phi = getphi(bulkevol)

    B1_cache  = getB1(cache)
    B2_cache  = getB2(cache)
    G_cache   = getG(cache)
    phi_cache = getphi(cache)

    # the loops in these functions are threaded, so it's probably not worth it
    # to @spawn here
    apply_dissipation_3D!( B1_cache,  B1, sys)
    apply_dissipation_3D!( B2_cache,  B2, sys)
    apply_dissipation_3D!(  G_cache,   G, sys)
    apply_dissipation_3D!(phi_cache, phi, sys)

    copyto!( B1,  B1_cache)
    copyto!( B2,  B2_cache)
    copyto!(  G,   G_cache)
    copyto!(phi, phi_cache)

    nothing
end

# exponential filtering
function (filters::Filters)(bulkevol::BulkEvolved)
    @sync begin
        @spawn filters.exp_filter(bulkevol.B1)
        @spawn filters.exp_filter(bulkevol.B2)
        @spawn filters.exp_filter(bulkevol.G)
        @spawn filters.exp_filter(bulkevol.phi)
    end
    nothing
end

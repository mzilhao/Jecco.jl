
using LinearAlgebra

# FIXME: need to set Sd and also A

function solve_nested_g1!(bulk::BulkVars, boundary::BoundaryVars, sys)
    coords = sys.coords
    derivs = sys.derivs
    uderiv = derivs[1]

    uu, xx, yy = Vivi.xx(coords)
    Nu = length(uu)

    Du_phi    = Vivi.D(bulk.phi, derivs[1], 1)
    Dxx_phi   = Vivi.D2(bulk.phi, derivs[2], 2)
    Dyy_phi   = Vivi.D2(bulk.phi, derivs[3], 3)


    # TODO: parallelize here
    @fastmath @inbounds for j in eachindex(yy)
        @fastmath @inbounds for i in eachindex(xx)
            # TODO: pre-allocate
            A_mat = zeros(Nu, Nu)
            b_vec = zeros(Nu)
            ABCS  = zeros(4)
            vars  = AllVars{Float64}()

            @fastmath @inbounds @simd for a in eachindex(uu)
                vars.u       = uu[a]
                vars.Sd_d0   = bulk.Sd[a,i,j]
                vars.phi_d0  = bulk.phi[a,i,j]
                vars.phi_du  = Du_phi[a,i,j]
                vars.phi_dxx = Dxx_phi[a,i,j]
                vars.phi_dyy = Dyy_phi[a,i,j]

                phig1_eq_coeff!(ABCS, vars)

                b_vec[a]     = -ABCS[4]

                @inbounds @simd for aa in eachindex(uu)
                    A_mat[a,aa] = ABCS[1] * uderiv.D2[a,aa] + ABCS[2] * uderiv.D[a,aa]
                end
                A_mat[a,a] += ABCS[3]
            end

            # boundary condition
            b_vec[1]    = bulk.phi[1,i,j]  # phi2[i,j]
            A_mat[1,:] .= 0.0
            A_mat[1,1]  = 1.0

            sol = view(bulk.phid, :, i, j)
            solve_lin_system!(sol, A_mat, b_vec)
        end
    end

    nothing
end

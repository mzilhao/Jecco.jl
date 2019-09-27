
using LinearAlgebra

## TODO
# Need signature Vf(phi::Real)
function Vf end
function Vfp end

function solve_lin_system!(sol, A_mat, b_vec)
    A_fact = lu!(A_mat)
    ldiv!(A_fact, b_vec)
    sol .= b_vec
    nothing
end


# assuming
# (A d_uu + B d_u + C Id) f = -S
function phig1_eq_coeff!(ABCS::Vector, vars::AllVars)
    u  = vars.u
    u2 = u * u
    u3 = u2 * u
    u4 = u2 * u2
    u5 = u3 * u2
    u6 = u4 * u2
    u9 = u6 * u3

    phi  = vars.phi_d0
    phi3 = phi * phi * phi
    phi4 = phi * phi3

    ABCS[1] = 0
    ABCS[2] = 9.0 * u
    ABCS[3] = 4.5

    ABCS[4] = -4.5 * phi - 27 * u4 * phi * vars.Sd_d0 -
        4 * u6 * phi3 * Vf(phi) - u9 * phi4 * Vfp(phi) +
        3 * u2 * (vars.phi_dxx + vars.phi_dyy) -
        4.5 * u * vars.phi_du - 9 * u5 * vars.Sd_d0 * vars.phi_du

    nothing
end

# FIXME: need to set Sd (and maybe make sure such cases are never
# uninitialized?). and then also A

function solve_phidg1!(bulk::BulkVars, sys)
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

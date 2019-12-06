
using LinearAlgebra

function solve_lin_system!(sol, A_mat, b_vec)
    A_fact = lu!(A_mat)
    ldiv!(A_fact, b_vec)
    sol .= b_vec
    nothing
end

struct Nested{S,D,T<:Real}
    sys     :: S
    uu      :: Vector{T}
    xx      :: Vector{T}
    yy      :: Vector{T}
    Du_phi  :: D
    Dxx_phi :: D
    Dyy_phi :: D
    Du_phid :: D
    A_mat   :: Matrix{T}
    b_vec   :: Vector{T}
    vars    :: AllVars{T}
end
function Nested(sys::System)
    coords = sys.coords

    uu, xx, yy = Vivi.xx(coords)
    Nu = length(uu)
    Nx = length(xx)
    Ny = length(yy)

    Du_phi    = zeros(Nu, Nx, Ny)
    Dxx_phi   = zeros(Nu, Nx, Ny)
    Dyy_phi   = zeros(Nu, Nx, Ny)
    Du_phid   = zeros(Nu, Nx, Ny)

    A_mat = zeros(Nu, Nu)
    b_vec = zeros(Nu)
    vars  = AllVars{eltype(A_mat)}()

    Nested{typeof(sys), typeof(Du_phi),
           eltype(A_mat)}(sys, uu, xx, yy, Du_phi, Dxx_phi, Dyy_phi, Du_phid,
                          A_mat, b_vec, vars)
end

Nested(systems::Vector) = [Nested(sys) for sys in systems]


function nested_g1!(nested::Nested, bulk::BulkVars, BC::BulkVars)
    sys  = nested.sys
    uu   = nested.uu
    xx   = nested.xx
    yy   = nested.yy

    Du_phi  = nested.Du_phi
    Dxx_phi = nested.Dxx_phi
    Dyy_phi = nested.Dyy_phi
    Du_phid = nested.Du_phid

    A_mat   = nested.A_mat
    b_vec   = nested.b_vec
    vars    = nested.vars

    uderiv = sys.uderiv
    xderiv = sys.xderiv
    yderiv = sys.yderiv

    Vivi.D!(Du_phi, bulk.phi, uderiv, 1)
    Vivi.D2!(Dxx_phi, bulk.phi, xderiv, 2)
    Vivi.D2!(Dyy_phi, bulk.phi, yderiv, 3)

    ABCS  = zeros(4)

    # set Sd
    @fastmath @inbounds for j in eachindex(yy)
        @fastmath @inbounds for i in eachindex(xx)
            @fastmath @inbounds @simd for a in eachindex(uu)
                bulk.Sd[a,i,j] = BC.Sd[i,j]
            end
        end
    end


    # solve for phidg1

    # TODO: parallelize here
    @fastmath @inbounds for j in eachindex(yy)
        @fastmath @inbounds for i in eachindex(xx)

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
            b_vec[1]    = BC.phid[i,j]
            A_mat[1,:] .= 0.0
            A_mat[1,1]  = 1.0

            sol = view(bulk.phid, :, i, j)
            solve_lin_system!(sol, A_mat, b_vec)
        end
    end


    # set Ag1
    @fastmath @inbounds for j in eachindex(yy)
        @fastmath @inbounds for i in eachindex(xx)
            @fastmath @inbounds @simd for a in eachindex(uu)
                bulk.A[a,i,j] = BC.A[i,j]
            end
        end
    end


    # finally compute dphidt_g1

    Vivi.D!(Du_phid, bulk.phid, uderiv, 1)

    # TODO: parallelize here
    @fastmath @inbounds for j in eachindex(yy)
        @inbounds for i in eachindex(xx)
            @inbounds @simd for a in eachindex(uu)
                vars.u       = uu[a]

                vars.phi_d0  = bulk.phi[a,i,j]
                vars.phid_d0 = bulk.phid[a,i,j]
                vars.A_d0    = bulk.A[a,i,j]

                vars.phi_du  = Du_phi[a,i,j]
                vars.phid_du = Du_phid[a,i,j]

                if vars.u > 1.e-9
                    bulk.dphidt[a,i,j]  = dphig1dt(vars)
                else
                    bulk.dphidt[a,i,j]  = dphig1dt_u0(vars)
                end
            end
        end
    end

    nothing
end

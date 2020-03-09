
using Jecco
using Jecco.AdS5_3_1

import Base.Threads.@threads
import Base.Threads.@spawn
using LinearAlgebra

# par_base = ParamBase(
#     which_potential = "square",
# )

function IDtest0(sys::System)
    Nu, Nx, Ny = size(sys.grid)

    B1  = zeros(Nu, Nx, Ny)
    B2  = zeros(Nu, Nx, Ny)
    G   = zeros(Nu, Nx, Ny)
    phi = zeros(Nu, Nx, Ny)

    b14 = 0.01
    b24 = 0.02

    for j in 1:Ny
        for i in 1:Nx
            for a in 1:Nu
                u,x,y = sys.grid[a,i,j]
                B1[a,i,j]  = u^4 * b14
                B2[a,i,j]  = u^4 * b24
                phi[a,i,j] = 0.0
                G[a,i,j]   = 0.0
            end
        end
    end

    BulkVars(B1, B2, G, phi)
end


par_grid = ParamGrid(
    xmin        = -5.0,
    xmax        =  5.0,
    xnodes      =  128,
    ymin        = -5.0,
    ymax        =  5.0,
    ynodes      =  128,
    umin        =  0.01,
    umax        =  1.0,
    udomains    =  1,
    unodes      =  16,
    # unodes      =  32,
)


systems = Jecco.AdS5_3_1.create_systems(par_grid)

sys = systems[1]

# Nu, Nx, Ny = size(sys.grid)
Nu_, Nx, Ny = size(sys.grid)
u, x, y = sys.grid[:]

nested = Jecco.AdS5_3_1.Nested(sys)

bulk = IDtest0(sys)

BC  = BulkVars(Nx, Ny)
dBC = BulkVars(Nx, Ny)

u0 = u[1]

fx2_0 = 0.005
fy2_0 = 0.002

BC.S  .= 1.0/u0
dBC.S .= -1.0/(u0*u0)

BC.Fx .= fx2_0 * u0 * u0
BC.Fy .= fy2_0 * u0 * u0

dBC.Fx .= 2 * fx2_0 * u0
dBC.Fy .= 2 * fy2_0 * u0

Jecco.AdS5_3_1.solve_nested_outer!(bulk, BC, dBC, nested)

Du_B1   = nested.Du_B1
Du_B2   = nested.Du_B2
Du_G    = nested.Du_G
Du_phi  = nested.Du_phi
Du_S    = nested.Du_S
Duu_B1  = nested.Duu_B1
Duu_B2  = nested.Duu_B2
Duu_G   = nested.Duu_G
Duu_phi = nested.Duu_phi
Duu_S   = nested.Duu_S

Du  = nested.sys.Du
Duu = nested.sys.Duu

@sync begin
    @spawn mul!(Du_B1,  Du,  bulk.B1)
    @spawn mul!(Du_B2,  Du,  bulk.B2)
    @spawn mul!(Du_G,   Du,  bulk.G)
    @spawn mul!(Du_phi, Du,  bulk.phi)
    @spawn mul!(Duu_B1, Duu, bulk.B1)
    @spawn mul!(Duu_B2, Duu, bulk.B2)
    @spawn mul!(Duu_G,  Duu, bulk.G)
    @spawn mul!(Duu_phi,Duu, bulk.phi)
end


@sync begin
    @spawn mul!(Du_S,   Du,  bulk.S)
    @spawn mul!(Duu_S,  Duu, bulk.S)
end


# T = Float64


#             AA = zeros(2,2)
#             BB = zeros(2,2)
#             CC = zeros(2,2)
#             SS = zeros(2)

#             A_mat2 = zeros(T, 2*Nu, 2*Nu)
#             b_vec2 = zeros(T, 2*Nu)

#             vars = Jecco.AdS5_3_1.FxyVars{T}()

function Fxy(bulk::BulkVars, BC::BulkVars, dBC::BulkVars, nested::Jecco.AdS5_3_1.Nested)
    sys  = nested.sys
    uu   = nested.uu
    xx   = nested.xx
    yy   = nested.yy

    Nu = length(uu)

    Du_B1   = nested.Du_B1
    Du_B2   = nested.Du_B2
    Du_G    = nested.Du_G
    Du_phi  = nested.Du_phi
    Du_S    = nested.Du_S
    Duu_B1  = nested.Duu_B1
    Duu_B2  = nested.Duu_B2
    Duu_G   = nested.Duu_G
    Duu_phi = nested.Duu_phi
    Duu_S   = nested.Duu_S

    aux_acc = nested.aux_acc

    Du  = sys.Du
    Duu = sys.Duu
    Dx  = sys.Dx
    Dxx = sys.Dxx
    Dy  = sys.Dy
    Dyy = sys.Dyy

    @fastmath @inbounds @threads for j in eachindex(yy)
        @inbounds for i in eachindex(xx)

            id  = Threads.threadid()
            aux = aux_acc[id]

            @inbounds @simd for a in eachindex(uu)
                u          = uu[a]
                u2         = u * u
                u3         = u * u2
                u4         = u2 * u2

                aux.varsFxy.u     = u

                # maybe it's worth to make some structs (or macros), here...

                aux.varsFxy.B1    = bulk.B1[a,i,j]
                aux.varsFxy.B1p   = -u2 * Du_B1[a,i,j]
                aux.varsFxy.B1_x  = Dx(bulk.B1, a,i,j)
                aux.varsFxy.B1_y  = Dy(bulk.B1, a,i,j)
                aux.varsFxy.B1pp  = 2*u3 * Du_B1[a,i,j] + u4 * Duu_B1[a,i,j]
                aux.varsFxy.B1p_x = -u2 * Dx(Du_B1, a,i,j)
                aux.varsFxy.B1p_y = -u2 * Dy(Du_B1, a,i,j)

                aux.varsFxy.B2    = bulk.B2[a,i,j]
                aux.varsFxy.B2p   = -u2 * Du_B2[a,i,j]
                aux.varsFxy.B2_x  = Dx(bulk.B2, a,i,j)
                aux.varsFxy.B2_y  = Dy(bulk.B2, a,i,j)
                aux.varsFxy.B2pp  = 2*u3 * Du_B2[a,i,j] + u4 * Duu_B2[a,i,j]
                aux.varsFxy.B2p_x = -u2 * Dx(Du_B2, a,i,j)
                aux.varsFxy.B2p_y = -u2 * Dy(Du_B2, a,i,j)

                aux.varsFxy.G     = bulk.G[a,i,j]
                aux.varsFxy.Gp    = -u2 * Du_G[a,i,j]
                aux.varsFxy.G_x   = Dx(bulk.G, a,i,j)
                aux.varsFxy.G_y   = Dy(bulk.G, a,i,j)
                aux.varsFxy.Gpp   = 2*u3 * Du_G[a,i,j] + u4 * Duu_G[a,i,j]
                aux.varsFxy.Gp_x  = -u2 * Dx(Du_G, a,i,j)
                aux.varsFxy.Gp_y  = -u2 * Dy(Du_G, a,i,j)

                aux.varsFxy.phi   = bulk.phi[a,i,j]
                aux.varsFxy.phip  = -u2 * Du_phi[a,i,j]
                aux.varsFxy.phi_x = Dx(bulk.phi, a,i,j)
                aux.varsFxy.phi_y = Dy(bulk.phi, a,i,j)
                # aux.varsFxy.phipp   = 2*u3 * Du_phi[a,i,j] + u4 * Duu_phi[a,i,j]
                # aux.varsFxy.phip_x  = -u2 * Dx(Du_phi, a,i,j)
                # aux.varsFxy.phip_y  = -u2 * Dy(Du_phi, a,i,j)

                aux.varsFxy.S     = bulk.S[a,i,j]
                aux.varsFxy.Sp    = -u2 * Du_S[a,i,j]
                aux.varsFxy.S_x   = Dx(bulk.S, a,i,j)
                aux.varsFxy.S_y   = Dy(bulk.S, a,i,j)
                aux.varsFxy.Spp   = 2*u3 * Du_S[a,i,j] + u4 * Duu_S[a,i,j]
                aux.varsFxy.Sp_x  = -u2 * Dx(Du_S, a,i,j)
                aux.varsFxy.Sp_y  = -u2 * Dy(Du_S, a,i,j)

                Jecco.AdS5_3_1.Fxy_outer_eq_coeff!(aux.AA, aux.BB, aux.CC, aux.SS, aux.varsFxy)

                @inbounds @simd for aa in eachindex(uu)
                    aux.A_mat2[a,aa]         = aux.AA[1,1] .* Duu[a,aa] .+ aux.BB[1,1] .* Du[a,aa] .+ aux.CC[1,1]
                    aux.A_mat2[a,aa+Nu]      = aux.AA[1,2] .* Duu[a,aa] .+ aux.BB[1,2] .* Du[a,aa] .+ aux.CC[1,2]
                    aux.A_mat2[a+Nu,aa]      = aux.AA[2,1] .* Duu[a,aa] .+ aux.BB[2,1] .* Du[a,aa] .+ aux.CC[2,1]
                    aux.A_mat2[a+Nu,aa+Nu]   = aux.AA[2,2] .* Duu[a,aa] .+ aux.BB[2,2] .* Du[a,aa] .+ aux.CC[2,2]
                end

                aux.b_vec2[a]    = -aux.SS[1]
                aux.b_vec2[a+Nu] = -aux.SS[2]

            end

            # BC

            aux.b_vec2[1]    = BC.Fx[i,j]
            aux.b_vec2[1+Nu] = BC.Fy[i,j]

            aux.A_mat2[1,:]   .= 0.0
            aux.A_mat2[1,1]    = 1.0
            aux.A_mat2[1,1+Nu] = 1.0

            aux.b_vec2[Nu]   = dBC.Fx[i,j]
            aux.b_vec2[2*Nu] = dBC.Fy[i,j]

            @inbounds @simd for aa in eachindex(uu)
                aux.A_mat2[end,aa]    = Du[1,aa]
                aux.A_mat2[end,aa+Nu] = Du[1,aa]
            end

            Jecco.AdS5_3_1.solve_lin_system!(aux.sol2, aux.A_mat2, aux.b_vec2)

            @inbounds @simd for aa in eachindex(uu)
                bulk.Fx[aa,i,j] = aux.sol2[aa]
                bulk.Fy[aa,i,j] = aux.sol2[aa+Nu]
            end

        end
    end

    nothing
end

# Dx = nested.sys.Dx
# Dy = nested.sys.Dy


Fxy(bulk, BC, dBC, nested)

# @inbounds @simd for a in eachindex(uu)
#     u              = uu[a]
#     aux.vars.u     = u

#     aux.vars.B1p   = -u*u * Du_B1[a,i,j]
#     aux.vars.B2p   = -u*u * Du_B2[a,i,j]

#     aux.vars.G     = bulk.G[a,i,j]
#     aux.vars.Gp    = -u*u * Du_G[a,i,j]

#     aux.vars.phip  = -u*u * Du_phi[a,i,j]

#     S_outer_eq_coeff!(aux.ABCS, aux.vars)

#     aux.b_vec[a]   = -aux.ABCS[4]
#     @inbounds @simd for aa in eachindex(uu)
#         aux.A_mat[a,aa] = aux.ABCS[1] * Duu.D[a,aa] + aux.ABCS[2] * Du.D[a,aa]
#     end
#     aux.A_mat[a,a] += aux.ABCS[3]
# end

# aux.b_vec[1]    = BC.S[i,j]
# aux.A_mat[1,:] .= 0.0
# aux.A_mat[1,1]  = 1.0

# aux.b_vec[end]    = dBC.S[i,j]
# aux.A_mat[end,:]  = Du.D[1,:]

# sol = view(bulk.S, :, i, j)
# solve_lin_system!(sol, aux.A_mat, aux.b_vec)

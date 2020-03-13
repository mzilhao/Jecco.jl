
import Base.Threads.@threads
import Base.Threads.@spawn
using LinearAlgebra

function solve_lin_system!(sol, A_mat, b_vec)
    A_fact = lu!(A_mat)
    ldiv!(A_fact, b_vec)
    sol .= b_vec
    nothing
end

struct Aux{T<:Real}
    A_mat   :: Matrix{T}
    b_vec   :: Vector{T}
    ABCS    :: Vector{T}
    vars    :: AllVars{T}

    A_mat2  :: Matrix{T}
    b_vec2  :: Vector{T}
    sol2    :: Vector{T}
    AA      :: Matrix{T}
    BB      :: Matrix{T}
    CC      :: Matrix{T}
    SS      :: Vector{T}
    varsFxy :: FxyVars{T}

    function Aux{T}(N::Int) where {T<:Real}
        A_mat  = zeros(T, N, N)
        b_vec  = zeros(T, N)
        ABCS   = zeros(T, 4)
        vars   = AllVars{T}()

        A_mat2 = zeros(T, 2*N, 2*N)
        b_vec2 = zeros(T, 2*N)
        sol2   = zeros(T, 2*N)

        AA = zeros(2,2)
        BB = zeros(2,2)
        CC = zeros(2,2)
        SS = zeros(2)

        varsFxy = FxyVars{T}()

        new(A_mat, b_vec, ABCS, vars, A_mat2, b_vec2, sol2, AA, BB, CC, SS, varsFxy)
    end
end

struct Nested{S,D,T<:Real}
    sys     :: S
    uu      :: Vector{T}
    xx      :: Vector{T}
    yy      :: Vector{T}
    Du_B1   :: D
    Du_B2   :: D
    Du_G    :: D
    Du_phi  :: D
    Du_S    :: D
    Duu_B1  :: D
    Duu_B2  :: D
    Duu_G   :: D
    Duu_phi :: D
    Duu_S   :: D
    aux_acc :: Vector{Aux{T}}
end
function Nested(sys::System)
    Nu, Nx, Ny = size(sys.grid)
    uu, xx, yy = sys.grid[:]

    Du_B1    = zeros(Nu, Nx, Ny)
    Du_B2    = zeros(Nu, Nx, Ny)
    Du_G     = zeros(Nu, Nx, Ny)
    Du_phi   = zeros(Nu, Nx, Ny)
    Du_S     = zeros(Nu, Nx, Ny)
    Duu_B1   = zeros(Nu, Nx, Ny)
    Duu_B2   = zeros(Nu, Nx, Ny)
    Duu_G    = zeros(Nu, Nx, Ny)
    Duu_phi  = zeros(Nu, Nx, Ny)
    Duu_S    = zeros(Nu, Nx, Ny)

    nt = Threads.nthreads()
    # pre-allocate thread-local aux quantities
    aux_acc = [Aux{eltype(uu)}(Nu) for _ in 1:nt]

    Nested{typeof(sys),typeof(Du_B1),
           eltype(uu)}(sys, uu, xx, yy, Du_B1, Du_B2, Du_G, Du_phi, Du_S,
                       Duu_B1, Duu_B2, Duu_G, Duu_phi, Duu_S, aux_acc)
end

Nested(systems::Vector) = [Nested(sys) for sys in systems]


#= Notes

for each metric function there are radial ODEs at each z point. since they are
x,y-independent, they can all be solved independently and simultaneously. this is
achieved via the trivial Threads.@thread parallelisation below.

the matrix A_mat is obtained, for each equation, through

  A_mat = A D_uu + B D_u + C Id

which builds the differential operator, and then the top and bottom lines are
replaced to enforce the boundary conditions, as shown below. the vector b_vec is
simply built through

  b_vec = -S

the first and last entries are replaced to enforce the boundary conditions,
leading to the schematic form

  (  1   0   0  ...  0  ) ( x0 )   (  u0  )
  ( a10 a11 a12 ... a1N ) ( x1 )   (  b1  )
  ( ................... ) ( x2 ) = (  b2  )                               (*)
  ( ................... ) ( .. )   (  ..  )
  ( d00 d01 d02 ... d0N ) ( xN )   (  u'0 )

where N = Nu and d00, d01, d02, ..., d0N are the coefficients of the first
line in the first derivative operator D_u:

         ( d00 d01 d02 d03 ... d0N )
         ( d10 d11 d12 d13 ... d1N )
  D_u =  ( d20 d21 d22 d23 ... d2N )
         ( ....................... )
         ( dN0 dN1 dN2 dN3 ... dNN )

the first equation from the system (*) then enforces x0 = u0 while the last
equation enforces dx/du_{u=u0} = u'0, where u0 and u'0 need to be given for each
equation. the remaining equations from (*) enforce the differential equations
themselves.

note that this example is valid only for the second order ODEs. for the first
order ones, obviously, just one BC is needed. thus, we accordingly skip the step
of replacing the last line of the A_mat matrix and last entry of b_vec vector.

=#

function solve_nested_outer!(bulk::BulkVars, BC::BulkVars, dBC::BulkVars, nested::Nested)
    sys  = nested.sys
    uu   = nested.uu
    xx   = nested.xx
    yy   = nested.yy

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

    Nu = length(uu)

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

    # solve for S

    @fastmath @inbounds @threads for j in eachindex(yy)
        @inbounds for i in eachindex(xx)
            id  = Threads.threadid()
            aux = aux_acc[id]

            @inbounds @simd for a in eachindex(uu)
                u              = uu[a]
                aux.vars.u     = u

                aux.vars.B1p   = -u*u * Du_B1[a,i,j]
                aux.vars.B2p   = -u*u * Du_B2[a,i,j]

                aux.vars.G     = bulk.G[a,i,j]
                aux.vars.Gp    = -u*u * Du_G[a,i,j]

                aux.vars.phip  = -u*u * Du_phi[a,i,j]

                S_outer_eq_coeff!(aux.ABCS, aux.vars)

                aux.b_vec[a]   = -aux.ABCS[4]
                @inbounds @simd for aa in eachindex(uu)
                    aux.A_mat[a,aa] = aux.ABCS[1] * Duu[a,aa] + aux.ABCS[2] * Du[a,aa]
                end
                aux.A_mat[a,a] += aux.ABCS[3]
            end

            aux.b_vec[1]    = BC.S[i,j]
            aux.A_mat[1,:] .= 0.0
            aux.A_mat[1,1]  = 1.0

            aux.b_vec[end]    = dBC.S[i,j]
            @inbounds @simd for aa in eachindex(uu)
                aux.A_mat[end,aa]  = Du[1,aa]
            end

            sol = view(bulk.S, :, i, j)
            solve_lin_system!(sol, aux.A_mat, aux.b_vec)
        end
    end

    # take u-derivatives of S
    @sync begin
        @spawn mul!(Du_S,   Du,  bulk.S)
        @spawn mul!(Duu_S,  Duu, bulk.S)
    end


    # solve for Fx and Fy

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

                Fxy_outer_eq_coeff!(aux.AA, aux.BB, aux.CC, aux.SS, aux.varsFxy)

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

            aux.b_vec2[1]     = BC.Fx[i,j]
            aux.A_mat2[1,:]  .= 0.0
            aux.A_mat2[1,1]   = 1.0

            aux.b_vec2[1+Nu]      = BC.Fy[i,j]
            aux.A_mat2[1+Nu,:]   .= 0.0
            aux.A_mat2[1+Nu,1+Nu] = 1.0

            aux.b_vec2[Nu]   = dBC.Fx[i,j]
            aux.b_vec2[2*Nu] = dBC.Fy[i,j]

            aux.A_mat2[Nu,:]   .= 0.0
            aux.A_mat2[2*Nu,:] .= 0.0
            @inbounds @simd for aa in eachindex(uu)
                aux.A_mat2[Nu,aa]      = Du[1,aa]
                aux.A_mat2[2*Nu,aa+Nu] = Du[1,aa]
            end

            solve_lin_system!(aux.sol2, aux.A_mat2, aux.b_vec2)

            @inbounds @simd for aa in eachindex(uu)
                bulk.Fx[aa,i,j] = aux.sol2[aa]
                bulk.Fy[aa,i,j] = aux.sol2[aa+Nu]
            end

        end
    end



    # # finally compute dphidt_g1

    # mul!(Du_phid, Du, bulk.phid)

    # @fastmath @inbounds @threads for j in eachindex(yy)
    #     @inbounds for i in eachindex(xx)
    #         id  = Threads.threadid()
    #         aux = aux_acc[id]
    #         @inbounds @simd for a in eachindex(uu)
    #             aux.vars.u       = uu[a]

    #             aux.vars.phi_d0  = bulk.phi[a,i,j]
    #             aux.vars.phid_d0 = bulk.phid[a,i,j]
    #             aux.vars.A_d0    = bulk.A[a,i,j]

    #             aux.vars.phi_du  = Du_phi[a,i,j]
    #             aux.vars.phid_du = Du_phid[a,i,j]

    #             if aux.vars.u > 1.e-9
    #                 bulk.dphidt[a,i,j]  = dphig1dt(aux.vars)
    #             else
    #                 bulk.dphidt[a,i,j]  = dphig1dt_u0(aux.vars)
    #             end
    #         end
    #     end
    # end

    nothing
end


# function solve_nested_g1!(bulk::BulkVars, BC::BulkVars, boundary::BoundaryVars,
#                           nested::Nested)
#     # u=0 boundary
#     BC.Sd   .= 0.5 * boundary.a4
#     BC.phid .= bulk.phi[1,:,:] # phi2
#     BC.A    .= boundary.a4

#     solve_nested_g1!(bulk, BC, nested)

#     nothing
# end

# function solve_nested_g1!(bulks::Vector, BCs::Vector, boundary::BoundaryVars,
#                           nesteds::Vector)
#     Nsys = length(nesteds)

#     # u=0 boundary
#     BCs[1].Sd   .= 0.5 * boundary.a4
#     BCs[1].phid .= bulks[1].phi[1,:,:] # phi2
#     BCs[1].A    .= boundary.a4

#     for i in 1:Nsys-1
#         solve_nested_g1!(bulks[i], BCs[i], nesteds[i])
#         BCs[i+1] = bulks[i][end,:,:]
#     end
#     solve_nested_g1!(bulks[Nsys], BCs[Nsys], nesteds[Nsys])

#     # sync boundary points. note: in a more general situation we may need to
#     # check the characteristic speeds (in this case we just know where the
#     # horizon is)
#     for i in 1:Nsys-1
#         bulks[i].dphidt[end,:,:] .= bulks[i+1].dphidt[1,:,:]
#     end

#     nothing
# end


# function solve_nested_g1(phi::Array{<:Number,N}, sys::System) where {N}
#     a4 = -ones2D(sys)
#     boundary = BoundaryVars(a4)

#     bulk = BulkVars(phi)
#     BC = bulk[1,:,:]

#     nested = Nested(sys)

#     solve_nested_g1!(bulk, BC, boundary, nested)
#     bulk
# end

# function solve_nested_g1(phis::Vector, systems::Vector)
#     a4 = -ones2D(systems[1])
#     boundary = BoundaryVars(a4)

#     bulks = BulkVars(phis)
#     phis_slice  = [phi[1,:,:] for phi in phis]
#     BCs  = BulkVars(phis_slice)

#     Nsys    = length(systems)
#     nesteds = Nested(systems)

#     solve_nested_g1!(bulks, BCs, boundary, nesteds)
#     bulks
# end

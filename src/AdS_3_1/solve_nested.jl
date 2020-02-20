
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

    function Aux{T}(N::Int) where {T<:Real}
        A_mat  = zeros(T, N, N)
        b_vec  = zeros(T, N)
        ABCS   = zeros(T, 4)
        vars   = AllVars{T}()

        new(A_mat, b_vec, ABCS, vars)
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
    aux_acc :: Vector{Aux{T}}
end
function Nested(sys::System)
    Nu, Nx, Ny = size(sys.grid)
    uu, xx, yy = sys.grid[:]

    Du_B1    = zeros(Nu, Nx, Ny)
    Du_B2    = zeros(Nu, Nx, Ny)
    Du_G     = zeros(Nu, Nx, Ny)
    Du_phi   = zeros(Nu, Nx, Ny)

    nt = Threads.nthreads()
    # pre-allocate thread-local aux quantities
    aux_acc = [Aux{eltype(uu)}(Nu) for _ in 1:nt]

    Nested{typeof(sys),typeof(Du_phi), eltype(uu)}(sys, uu, xx, yy, Du_B1, Du_B2,
                                                   Du_G, Du_phi, aux_acc)
end

Nested(systems::Vector) = [Nested(sys) for sys in systems]


#= Notes

for each metric function there are radial ODEs at each z point. since they are
z-independent, they can all be solved independently and simultaneously. this is
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

function solve_nested_outer!(bulk::BulkVars, BC::BulkVars, nested::Nested)
    sys  = nested.sys
    uu   = nested.uu
    xx   = nested.xx
    yy   = nested.yy

    Du_B1  = nested.Du_B1
    Du_B2  = nested.Du_B2
    Du_G   = nested.Du_G
    Du_phi = nested.Du_phi

    aux_acc = nested.aux_acc

    Du  = sys.Du
    Duu = sys.Duu
    # Dx  = sys.Dx
    Dxx = sys.Dxx
    # Dy  = sys.Dy
    Dyy = sys.Dyy

    # TODO @spawn here ?
    mul!(Du_B1,  Du, bulk.B1)
    mul!(Du_B2,  Du, bulk.B2)
    mul!(Du_G,   Du, bulk.G)
    mul!(Du_phi, Du, bulk.phi)

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

                aux.b_vec[a]     = -aux.ABCS[4]

                # TODO: replace Duu.D and Du.D with something more general and portable
                @inbounds @simd for aa in eachindex(uu)
                    aux.A_mat[a,aa] = aux.ABCS[1] * Duu.D[a,aa] + aux.ABCS[2] * Du.D[a,aa]
                end
                aux.A_mat[a,a] += aux.ABCS[3]
            end

            # boundary condition. FIXME
            dfunc_du_ui     = 0.0
            func_ui         = 0.01
            aux.b_vec[1]    = func_ui
            aux.A_mat[1,:] .= 0.0
            aux.A_mat[1,1]  = 1.0

            aux.b_vec[end]    = dfunc_du_ui
            aux.A_mat[end,:]  = Du.D[1,:]

            sol = view(bulk.S, :, i, j)
            solve_lin_system!(sol, aux.A_mat, aux.b_vec)
        end
    end


    # # set Ag1
    # @fastmath @inbounds for j in eachindex(yy)
    #     @inbounds for i in eachindex(xx)
    #         @inbounds @simd for a in eachindex(uu)
    #             bulk.A[a,i,j] = BC.A[i,j]
    #         end
    #     end
    # end


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

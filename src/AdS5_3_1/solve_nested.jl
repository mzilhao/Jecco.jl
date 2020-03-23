
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
    vars    :: AllVarsOuter{T}

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
        vars   = AllVarsOuter{T}()

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
    Du_Fx   :: D
    Du_Fy   :: D
    Duu_B1  :: D
    Duu_B2  :: D
    Duu_G   :: D
    Duu_phi :: D
    Duu_S   :: D
    Duu_Fx  :: D
    Duu_Fy  :: D
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
    Du_Fx    = zeros(Nu, Nx, Ny)
    Du_Fy    = zeros(Nu, Nx, Ny)
    Duu_B1   = zeros(Nu, Nx, Ny)
    Duu_B2   = zeros(Nu, Nx, Ny)
    Duu_G    = zeros(Nu, Nx, Ny)
    Duu_phi  = zeros(Nu, Nx, Ny)
    Duu_S    = zeros(Nu, Nx, Ny)
    Duu_Fx   = zeros(Nu, Nx, Ny)
    Duu_Fy   = zeros(Nu, Nx, Ny)

    nt = Threads.nthreads()
    # pre-allocate thread-local aux quantities
    aux_acc = [Aux{eltype(uu)}(Nu) for _ in 1:nt]

    Nested{typeof(sys),typeof(Du_B1),
           eltype(uu)}(sys, uu, xx, yy, Du_B1, Du_B2, Du_G, Du_phi, Du_S, Du_Fx, Du_Fy,
                       Duu_B1, Duu_B2, Duu_G, Duu_phi, Duu_S, Duu_Fx, Duu_Fy, aux_acc)
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

function solve_S_outer!(bulk::BulkVars, BC::BulkVars, dBC::BulkVars, nested::Nested)
    sys  = nested.sys
    uu   = nested.uu
    xx   = nested.xx
    yy   = nested.yy

    Du_B1   = nested.Du_B1
    Du_B2   = nested.Du_B2
    Du_G    = nested.Du_G
    Du_phi  = nested.Du_phi
    # Duu_B1  = nested.Duu_B1
    # Duu_B2  = nested.Duu_B2
    # Duu_G   = nested.Duu_G

    aux_acc = nested.aux_acc

    Du  = sys.Du
    Duu = sys.Duu
    # Dx  = sys.Dx
    # Dxx = sys.Dxx
    # Dy  = sys.Dy
    # Dyy = sys.Dyy

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

            # BC

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

    nothing
end

function solve_Fxy_outer!(bulk::BulkVars, BC::BulkVars, dBC::BulkVars, nested::Nested)
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

                aux.b_vec2[a]    = -aux.SS[1]
                aux.b_vec2[a+Nu] = -aux.SS[2]
                @inbounds @simd for aa in eachindex(uu)
                    aux.A_mat2[a,aa]         = aux.AA[1,1] * Duu[a,aa] + aux.BB[1,1] * Du[a,aa]
                    aux.A_mat2[a,aa+Nu]      = aux.AA[1,2] * Duu[a,aa] + aux.BB[1,2] * Du[a,aa]
                    aux.A_mat2[a+Nu,aa]      = aux.AA[2,1] * Duu[a,aa] + aux.BB[2,1] * Du[a,aa]
                    aux.A_mat2[a+Nu,aa+Nu]   = aux.AA[2,2] * Duu[a,aa] + aux.BB[2,2] * Du[a,aa]
                end
                aux.A_mat2[a,a]       += aux.CC[1,1]
                aux.A_mat2[a,a+Nu]    += aux.CC[1,2]
                aux.A_mat2[a+Nu,a]    += aux.CC[2,1]
                aux.A_mat2[a+Nu,a+Nu] += aux.CC[2,2]
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

    nothing
end

function solve_Sd_outer!(bulk::BulkVars, BC::BulkVars, nested::Nested)
    sys  = nested.sys
    uu   = nested.uu
    xx   = nested.xx
    yy   = nested.yy

    Du_B1   = nested.Du_B1
    Du_B2   = nested.Du_B2
    Du_G    = nested.Du_G
    Du_phi  = nested.Du_phi
    Du_S    = nested.Du_S
    Du_Fx   = nested.Du_Fx
    Du_Fy   = nested.Du_Fy
    Duu_B1  = nested.Duu_B1
    Duu_B2  = nested.Duu_B2
    Duu_G   = nested.Duu_G
    Duu_phi = nested.Duu_phi
    Duu_S   = nested.Duu_S
    Duu_Fx  = nested.Duu_Fx
    Duu_Fy  = nested.Duu_Fy

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

                B1p        = -u2 * Du_B1[a,i,j]
                B2p        = -u2 * Du_B2[a,i,j]
                Gp         = -u2 * Du_G[a,i,j]
                phip       = -u2 * Du_phi[a,i,j]
                Sp         = -u2 * Du_S[a,i,j]
                Fxp        = -u2 * Du_Fx[a,i,j]
                Fyp        = -u2 * Du_Fy[a,i,j]

                B1pp       = 2*u3 * Du_B1[a,i,j]  + u4 * Duu_B1[a,i,j]
                B2pp       = 2*u3 * Du_B2[a,i,j]  + u4 * Duu_B2[a,i,j]
                Gpp        = 2*u3 * Du_G[a,i,j]   + u4 * Duu_G[a,i,j]
                phipp      = 2*u3 * Du_phi[a,i,j] + u4 * Duu_phi[a,i,j]
                Spp        = 2*u3 * Du_S[a,i,j]   + u4 * Duu_S[a,i,j]
                Fxpp       = 2*u3 * Du_Fx[a,i,j]  + u4 * Duu_Fx[a,i,j]
                Fypp       = 2*u3 * Du_Fy[a,i,j]  + u4 * Duu_Fy[a,i,j]

                Fx         = bulk.Fx[a,i,j]
                Fy         = bulk.Fy[a,i,j]

                aux.vars.u     = u

                aux.vars.B1    = bulk.B1[a,i,j]
                aux.vars.B1p   = B1p
                aux.vars.B1t   = Dx(bulk.B1, a,i,j) - Fx * B1p
                aux.vars.B1h   = Dy(bulk.B1, a,i,j) - Fy * B1p

                aux.vars.B1tt  = Dxx(bulk.B1, a,i,j) - Dx(bulk.Fx, a,i,j) * B1p +
                    2.0 * Fx * u2 * Dx(Du_B1, a,i,j) + Fx * Fxp * B1p + Fx * Fx * B1pp
                aux.vars.B1hh  = Dyy(bulk.B1, a,i,j) - Dy(bulk.Fy, a,i,j) * B1p +
                    2.0 * Fy * u2 * Dy(Du_B1, a,i,j) + Fy * Fyp * B1p + Fy * Fy * B1pp

                aux.vars.B1tp  = -u2 * Dx(Du_B1, a,i,j) - Fxp * B1p - Fx * B1pp
                aux.vars.B1hp  = -u2 * Dy(Du_B1, a,i,j) - Fyp * B1p - Fy * B1pp


                aux.vars.B2    = bulk.B2[a,i,j]
                aux.vars.B2p   = B2p
                aux.vars.B2t   = Dx(bulk.B2, a,i,j) - Fx * B2p
                aux.vars.B2h   = Dy(bulk.B2, a,i,j) - Fy * B2p

                aux.vars.B2tt  = Dxx(bulk.B2, a,i,j) - Dx(bulk.Fx, a,i,j) * B2p +
                    2.0 * Fx * u2 * Dx(Du_B2, a,i,j) + Fx * Fxp * B2p + Fx * Fx * B2pp
                aux.vars.B2hh  = Dyy(bulk.B2, a,i,j) - Dy(bulk.Fy, a,i,j) * B2p +
                    2.0 * Fy * u2 * Dy(Du_B2, a,i,j) + Fy * Fyp * B2p + Fy * Fy * B2pp

                aux.vars.B2tp  = -u2 * Dx(Du_B2, a,i,j) - Fxp * B2p - Fx * B2pp
                aux.vars.B2hp  = -u2 * Dy(Du_B2, a,i,j) - Fyp * B2p - Fy * B2pp


                aux.vars.G    = bulk.G[a,i,j]
                aux.vars.Gp   = Gp
                aux.vars.Gt   = Dx(bulk.G, a,i,j) - Fx * Gp
                aux.vars.Gh   = Dy(bulk.G, a,i,j) - Fy * Gp

                aux.vars.Gtt  = Dxx(bulk.G, a,i,j) - Dx(bulk.Fx, a,i,j) * Gp +
                    2.0 * Fx * u2 * Dx(Du_G, a,i,j) + Fx * Fxp * Gp + Fx * Fx * Gpp
                aux.vars.Ghh  = Dyy(bulk.G, a,i,j) - Dy(bulk.Fy, a,i,j) * Gp +
                    2.0 * Fy * u2 * Dy(Du_G, a,i,j) + Fy * Fyp * Gp + Fy * Fy * Gpp

                aux.vars.Gtp  = -u2 * Dx(Du_G, a,i,j) - Fxp * Gp - Fx * Gpp
                aux.vars.Ghp  = -u2 * Dy(Du_G, a,i,j) - Fyp * Gp - Fy * Gpp


                aux.vars.phi    = bulk.phi[a,i,j]
                aux.vars.phip   = phip
                aux.vars.phit   = Dx(bulk.phi, a,i,j) - Fx * phip
                aux.vars.phih   = Dy(bulk.phi, a,i,j) - Fy * phip

                aux.vars.phitt  = Dxx(bulk.phi, a,i,j) - Dx(bulk.Fx, a,i,j) * phip +
                    2.0 * Fx * u2 * Dx(Du_phi, a,i,j) + Fx * Fxp * phip + Fx * Fx * phipp
                aux.vars.phihh  = Dyy(bulk.phi, a,i,j) - Dy(bulk.Fy, a,i,j) * phip +
                    2.0 * Fy * u2 * Dy(Du_phi, a,i,j) + Fy * Fyp * phip + Fy * Fy * phipp

                aux.vars.phitp  = -u2 * Dx(Du_phi, a,i,j) - Fxp * phip - Fx * phipp
                aux.vars.phihp  = -u2 * Dy(Du_phi, a,i,j) - Fyp * phip - Fy * phipp


                aux.vars.S    = bulk.S[a,i,j]
                aux.vars.Sp   = Sp
                aux.vars.St   = Dx(bulk.S, a,i,j) - Fx * Sp
                aux.vars.Sh   = Dy(bulk.S, a,i,j) - Fy * Sp

                aux.vars.Stt  = Dxx(bulk.S, a,i,j) - Dx(bulk.Fx, a,i,j) * Sp +
                    2.0 * Fx * u2 * Dx(Du_S, a,i,j) + Fx * Fxp * Sp + Fx * Fx * Spp
                aux.vars.Shh  = Dyy(bulk.S, a,i,j) - Dy(bulk.Fy, a,i,j) * Sp +
                    2.0 * Fy * u2 * Dy(Du_S, a,i,j) + Fy * Fyp * Sp + Fy * Fy * Spp

                aux.vars.Stp  = -u2 * Dx(Du_S, a,i,j) - Fxp * Sp - Fx * Spp
                aux.vars.Shp  = -u2 * Dy(Du_S, a,i,j) - Fyp * Sp - Fy * Spp


                aux.vars.Fx    = bulk.Fx[a,i,j]
                aux.vars.Fxp   = Fxp
                aux.vars.Fxt   = Dx(bulk.Fx, a,i,j) - Fx * Fxp
                aux.vars.Fxh   = Dy(bulk.Fx, a,i,j) - Fy * Fxp

                aux.vars.Fxtp  = -u2 * Dx(Du_Fx, a,i,j) - Fxp * Fxp - Fx * Fxpp
                aux.vars.Fxhp  = -u2 * Dy(Du_Fx, a,i,j) - Fyp * Fxp - Fy * Fxpp


                aux.vars.Fy    = bulk.Fy[a,i,j]
                aux.vars.Fyp   = Fyp
                aux.vars.Fyt   = Dx(bulk.Fy, a,i,j) - Fy * Fyp
                aux.vars.Fyh   = Dy(bulk.Fy, a,i,j) - Fy * Fyp

                aux.vars.Fytp  = -u2 * Dx(Du_Fy, a,i,j) - Fyp * Fyp - Fy * Fypp
                aux.vars.Fyhp  = -u2 * Dy(Du_Fy, a,i,j) - Fyp * Fyp - Fy * Fypp


                aux.vars.B2th = Dx(Dy, bulk.B2, a,i,j) - Dy(bulk.Fx, a,i,j) * B2p +
                    Fx * u2 * Dy(Du_B2, a,i,j) + Fy * u2 * Dx(Du_B2, a,i,j) +
                    Fy * Fxp * B2p + Fx * Fy * B2pp

                aux.vars.Gth = Dx(Dy, bulk.G, a,i,j) - Dy(bulk.Fx, a,i,j) * Gp +
                    Fx * u2 * Dy(Du_G, a,i,j) + Fy * u2 * Dx(Du_G, a,i,j) +
                    Fy * Fxp * Gp + Fx * Fy * Gpp

                aux.vars.Sth = Dx(Dy, bulk.S, a,i,j) - Dy(bulk.Fx, a,i,j) * Sp +
                    Fx * u2 * Dy(Du_S, a,i,j) + Fy * u2 * Dx(Du_S, a,i,j) +
                    Fy * Fxp * Sp + Fx * Fy * Spp

                Sd_outer_eq_coeff!(aux.ABCS, aux.vars)

                aux.b_vec[a]   = -aux.ABCS[4]
                @inbounds @simd for aa in eachindex(uu)
                    aux.A_mat[a,aa] = aux.ABCS[1] * Duu[a,aa] + aux.ABCS[2] * Du[a,aa]
                end
                aux.A_mat[a,a] += aux.ABCS[3]
            end

            # BC (first order equation)

            aux.b_vec[1]    = BC.Sd[i,j]
            aux.A_mat[1,:] .= 0.0
            aux.A_mat[1,1]  = 1.0

            sol = view(bulk.Sd, :, i, j)
            solve_lin_system!(sol, aux.A_mat, aux.b_vec)
        end
    end

    nothing
end

function solve_B2d_outer!(bulk::BulkVars, BC::BulkVars, nested::Nested)
    sys  = nested.sys
    uu   = nested.uu
    xx   = nested.xx
    yy   = nested.yy

    Du_B1   = nested.Du_B1
    Du_B2   = nested.Du_B2
    Du_G    = nested.Du_G
    Du_phi  = nested.Du_phi
    Du_S    = nested.Du_S
    Du_Fx   = nested.Du_Fx
    Du_Fy   = nested.Du_Fy
    Duu_B1  = nested.Duu_B1
    Duu_B2  = nested.Duu_B2
    Duu_G   = nested.Duu_G
    Duu_phi = nested.Duu_phi
    Duu_S   = nested.Duu_S
    Duu_Fx  = nested.Duu_Fx
    Duu_Fy  = nested.Duu_Fy

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

                B1p        = -u2 * Du_B1[a,i,j]
                B2p        = -u2 * Du_B2[a,i,j]
                Gp         = -u2 * Du_G[a,i,j]
                phip       = -u2 * Du_phi[a,i,j]
                Sp         = -u2 * Du_S[a,i,j]
                Fxp        = -u2 * Du_Fx[a,i,j]
                Fyp        = -u2 * Du_Fy[a,i,j]

                B1pp       = 2*u3 * Du_B1[a,i,j]  + u4 * Duu_B1[a,i,j]
                B2pp       = 2*u3 * Du_B2[a,i,j]  + u4 * Duu_B2[a,i,j]
                Gpp        = 2*u3 * Du_G[a,i,j]   + u4 * Duu_G[a,i,j]
                phipp      = 2*u3 * Du_phi[a,i,j] + u4 * Duu_phi[a,i,j]
                Spp        = 2*u3 * Du_S[a,i,j]   + u4 * Duu_S[a,i,j]
                Fxpp       = 2*u3 * Du_Fx[a,i,j]  + u4 * Duu_Fx[a,i,j]
                Fypp       = 2*u3 * Du_Fy[a,i,j]  + u4 * Duu_Fy[a,i,j]

                Fx         = bulk.Fx[a,i,j]
                Fy         = bulk.Fy[a,i,j]

                aux.vars.u     = u

                aux.vars.B1    = bulk.B1[a,i,j]
                aux.vars.B1p   = B1p
                aux.vars.B1t   = Dx(bulk.B1, a,i,j) - Fx * B1p
                aux.vars.B1h   = Dy(bulk.B1, a,i,j) - Fy * B1p

                aux.vars.B1tt  = Dxx(bulk.B1, a,i,j) - Dx(bulk.Fx, a,i,j) * B1p +
                    2.0 * Fx * u2 * Dx(Du_B1, a,i,j) + Fx * Fxp * B1p + Fx * Fx * B1pp
                aux.vars.B1hh  = Dyy(bulk.B1, a,i,j) - Dy(bulk.Fy, a,i,j) * B1p +
                    2.0 * Fy * u2 * Dy(Du_B1, a,i,j) + Fy * Fyp * B1p + Fy * Fy * B1pp

                aux.vars.B1tp  = -u2 * Dx(Du_B1, a,i,j) - Fxp * B1p - Fx * B1pp
                aux.vars.B1hp  = -u2 * Dy(Du_B1, a,i,j) - Fyp * B1p - Fy * B1pp


                aux.vars.B2    = bulk.B2[a,i,j]
                aux.vars.B2p   = B2p
                aux.vars.B2t   = Dx(bulk.B2, a,i,j) - Fx * B2p
                aux.vars.B2h   = Dy(bulk.B2, a,i,j) - Fy * B2p

                aux.vars.B2tt  = Dxx(bulk.B2, a,i,j) - Dx(bulk.Fx, a,i,j) * B2p +
                    2.0 * Fx * u2 * Dx(Du_B2, a,i,j) + Fx * Fxp * B2p + Fx * Fx * B2pp
                aux.vars.B2hh  = Dyy(bulk.B2, a,i,j) - Dy(bulk.Fy, a,i,j) * B2p +
                    2.0 * Fy * u2 * Dy(Du_B2, a,i,j) + Fy * Fyp * B2p + Fy * Fy * B2pp

                aux.vars.B2tp  = -u2 * Dx(Du_B2, a,i,j) - Fxp * B2p - Fx * B2pp
                aux.vars.B2hp  = -u2 * Dy(Du_B2, a,i,j) - Fyp * B2p - Fy * B2pp


                aux.vars.G    = bulk.G[a,i,j]
                aux.vars.Gp   = Gp
                aux.vars.Gt   = Dx(bulk.G, a,i,j) - Fx * Gp
                aux.vars.Gh   = Dy(bulk.G, a,i,j) - Fy * Gp

                aux.vars.Gtt  = Dxx(bulk.G, a,i,j) - Dx(bulk.Fx, a,i,j) * Gp +
                    2.0 * Fx * u2 * Dx(Du_G, a,i,j) + Fx * Fxp * Gp + Fx * Fx * Gpp
                aux.vars.Ghh  = Dyy(bulk.G, a,i,j) - Dy(bulk.Fy, a,i,j) * Gp +
                    2.0 * Fy * u2 * Dy(Du_G, a,i,j) + Fy * Fyp * Gp + Fy * Fy * Gpp

                aux.vars.Gtp  = -u2 * Dx(Du_G, a,i,j) - Fxp * Gp - Fx * Gpp
                aux.vars.Ghp  = -u2 * Dy(Du_G, a,i,j) - Fyp * Gp - Fy * Gpp


                aux.vars.phi    = bulk.phi[a,i,j]
                aux.vars.phip   = phip
                aux.vars.phit   = Dx(bulk.phi, a,i,j) - Fx * phip
                aux.vars.phih   = Dy(bulk.phi, a,i,j) - Fy * phip

                aux.vars.phitt  = Dxx(bulk.phi, a,i,j) - Dx(bulk.Fx, a,i,j) * phip +
                    2.0 * Fx * u2 * Dx(Du_phi, a,i,j) + Fx * Fxp * phip + Fx * Fx * phipp
                aux.vars.phihh  = Dyy(bulk.phi, a,i,j) - Dy(bulk.Fy, a,i,j) * phip +
                    2.0 * Fy * u2 * Dy(Du_phi, a,i,j) + Fy * Fyp * phip + Fy * Fy * phipp

                aux.vars.phitp  = -u2 * Dx(Du_phi, a,i,j) - Fxp * phip - Fx * phipp
                aux.vars.phihp  = -u2 * Dy(Du_phi, a,i,j) - Fyp * phip - Fy * phipp


                aux.vars.S    = bulk.S[a,i,j]
                aux.vars.Sp   = Sp
                aux.vars.St   = Dx(bulk.S, a,i,j) - Fx * Sp
                aux.vars.Sh   = Dy(bulk.S, a,i,j) - Fy * Sp

                aux.vars.Stt  = Dxx(bulk.S, a,i,j) - Dx(bulk.Fx, a,i,j) * Sp +
                    2.0 * Fx * u2 * Dx(Du_S, a,i,j) + Fx * Fxp * Sp + Fx * Fx * Spp
                aux.vars.Shh  = Dyy(bulk.S, a,i,j) - Dy(bulk.Fy, a,i,j) * Sp +
                    2.0 * Fy * u2 * Dy(Du_S, a,i,j) + Fy * Fyp * Sp + Fy * Fy * Spp

                aux.vars.Stp  = -u2 * Dx(Du_S, a,i,j) - Fxp * Sp - Fx * Spp
                aux.vars.Shp  = -u2 * Dy(Du_S, a,i,j) - Fyp * Sp - Fy * Spp


                aux.vars.Fx    = bulk.Fx[a,i,j]
                aux.vars.Fxp   = Fxp
                aux.vars.Fxt   = Dx(bulk.Fx, a,i,j) - Fx * Fxp
                aux.vars.Fxh   = Dy(bulk.Fx, a,i,j) - Fy * Fxp

                aux.vars.Fxtp  = -u2 * Dx(Du_Fx, a,i,j) - Fxp * Fxp - Fx * Fxpp
                aux.vars.Fxhp  = -u2 * Dy(Du_Fx, a,i,j) - Fyp * Fxp - Fy * Fxpp


                aux.vars.Fy    = bulk.Fy[a,i,j]
                aux.vars.Fyp   = Fyp
                aux.vars.Fyt   = Dx(bulk.Fy, a,i,j) - Fy * Fyp
                aux.vars.Fyh   = Dy(bulk.Fy, a,i,j) - Fy * Fyp

                aux.vars.Fytp  = -u2 * Dx(Du_Fy, a,i,j) - Fyp * Fyp - Fy * Fypp
                aux.vars.Fyhp  = -u2 * Dy(Du_Fy, a,i,j) - Fyp * Fyp - Fy * Fypp

                aux.vars.Sd    = bulk.Sd[a,i,j]

                aux.vars.B2th = Dx(Dy, bulk.B2, a,i,j) - Dy(bulk.Fx, a,i,j) * B2p +
                    Fx * u2 * Dy(Du_B2, a,i,j) + Fy * u2 * Dx(Du_B2, a,i,j) +
                    Fy * Fxp * B2p + Fx * Fy * B2pp

                aux.vars.Gth = Dx(Dy, bulk.G, a,i,j) - Dy(bulk.Fx, a,i,j) * Gp +
                    Fx * u2 * Dy(Du_G, a,i,j) + Fy * u2 * Dx(Du_G, a,i,j) +
                    Fy * Fxp * Gp + Fx * Fy * Gpp

                aux.vars.Sth = Dx(Dy, bulk.S, a,i,j) - Dy(bulk.Fx, a,i,j) * Sp +
                    Fx * u2 * Dy(Du_S, a,i,j) + Fy * u2 * Dx(Du_S, a,i,j) +
                    Fy * Fxp * Sp + Fx * Fy * Spp


                B2d_outer_eq_coeff!(aux.ABCS, aux.vars)

                aux.b_vec[a]   = -aux.ABCS[4]
                @inbounds @simd for aa in eachindex(uu)
                    aux.A_mat[a,aa] = aux.ABCS[1] * Duu[a,aa] + aux.ABCS[2] * Du[a,aa]
                end
                aux.A_mat[a,a] += aux.ABCS[3]
            end

            # BC (first order equation)

            aux.b_vec[1]    = BC.B2d[i,j]
            aux.A_mat[1,:] .= 0.0
            aux.A_mat[1,1]  = 1.0

            sol = view(bulk.B2d, :, i, j)
            solve_lin_system!(sol, aux.A_mat, aux.b_vec)
        end
    end

    nothing
end

function solve_B1dGd_outer!(bulk::BulkVars, BC::BulkVars, nested::Nested)
    sys  = nested.sys
    uu   = nested.uu
    xx   = nested.xx
    yy   = nested.yy

    Du_B1   = nested.Du_B1
    Du_B2   = nested.Du_B2
    Du_G    = nested.Du_G
    Du_phi  = nested.Du_phi
    Du_S    = nested.Du_S
    Du_Fx   = nested.Du_Fx
    Du_Fy   = nested.Du_Fy
    Duu_B1  = nested.Duu_B1
    Duu_B2  = nested.Duu_B2
    Duu_G   = nested.Duu_G
    Duu_phi = nested.Duu_phi
    Duu_S   = nested.Duu_S
    Duu_Fx  = nested.Duu_Fx
    Duu_Fy  = nested.Duu_Fy

    aux_acc = nested.aux_acc

    Du  = sys.Du
    Duu = sys.Duu
    Dx  = sys.Dx
    Dxx = sys.Dxx
    Dy  = sys.Dy
    Dyy = sys.Dyy

    Nu = length(uu)

    @fastmath @inbounds @threads for j in eachindex(yy)
        @inbounds for i in eachindex(xx)
            id  = Threads.threadid()
            aux = aux_acc[id]

            @inbounds @simd for a in eachindex(uu)
                u          = uu[a]
                u2         = u * u
                u3         = u * u2
                u4         = u2 * u2

                B1p        = -u2 * Du_B1[a,i,j]
                B2p        = -u2 * Du_B2[a,i,j]
                Gp         = -u2 * Du_G[a,i,j]
                phip       = -u2 * Du_phi[a,i,j]
                Sp         = -u2 * Du_S[a,i,j]
                Fxp        = -u2 * Du_Fx[a,i,j]
                Fyp        = -u2 * Du_Fy[a,i,j]

                B1pp       = 2*u3 * Du_B1[a,i,j]  + u4 * Duu_B1[a,i,j]
                B2pp       = 2*u3 * Du_B2[a,i,j]  + u4 * Duu_B2[a,i,j]
                Gpp        = 2*u3 * Du_G[a,i,j]   + u4 * Duu_G[a,i,j]
                phipp      = 2*u3 * Du_phi[a,i,j] + u4 * Duu_phi[a,i,j]
                Spp        = 2*u3 * Du_S[a,i,j]   + u4 * Duu_S[a,i,j]
                Fxpp       = 2*u3 * Du_Fx[a,i,j]  + u4 * Duu_Fx[a,i,j]
                Fypp       = 2*u3 * Du_Fy[a,i,j]  + u4 * Duu_Fy[a,i,j]

                Fx         = bulk.Fx[a,i,j]
                Fy         = bulk.Fy[a,i,j]

                aux.vars.u     = u

                aux.vars.B1    = bulk.B1[a,i,j]
                aux.vars.B1p   = B1p
                aux.vars.B1t   = Dx(bulk.B1, a,i,j) - Fx * B1p
                aux.vars.B1h   = Dy(bulk.B1, a,i,j) - Fy * B1p

                aux.vars.B1tt  = Dxx(bulk.B1, a,i,j) - Dx(bulk.Fx, a,i,j) * B1p +
                    2.0 * Fx * u2 * Dx(Du_B1, a,i,j) + Fx * Fxp * B1p + Fx * Fx * B1pp
                aux.vars.B1hh  = Dyy(bulk.B1, a,i,j) - Dy(bulk.Fy, a,i,j) * B1p +
                    2.0 * Fy * u2 * Dy(Du_B1, a,i,j) + Fy * Fyp * B1p + Fy * Fy * B1pp

                aux.vars.B1tp  = -u2 * Dx(Du_B1, a,i,j) - Fxp * B1p - Fx * B1pp
                aux.vars.B1hp  = -u2 * Dy(Du_B1, a,i,j) - Fyp * B1p - Fy * B1pp


                aux.vars.B2    = bulk.B2[a,i,j]
                aux.vars.B2p   = B2p
                aux.vars.B2t   = Dx(bulk.B2, a,i,j) - Fx * B2p
                aux.vars.B2h   = Dy(bulk.B2, a,i,j) - Fy * B2p

                aux.vars.B2tt  = Dxx(bulk.B2, a,i,j) - Dx(bulk.Fx, a,i,j) * B2p +
                    2.0 * Fx * u2 * Dx(Du_B2, a,i,j) + Fx * Fxp * B2p + Fx * Fx * B2pp
                aux.vars.B2hh  = Dyy(bulk.B2, a,i,j) - Dy(bulk.Fy, a,i,j) * B2p +
                    2.0 * Fy * u2 * Dy(Du_B2, a,i,j) + Fy * Fyp * B2p + Fy * Fy * B2pp

                aux.vars.B2tp  = -u2 * Dx(Du_B2, a,i,j) - Fxp * B2p - Fx * B2pp
                aux.vars.B2hp  = -u2 * Dy(Du_B2, a,i,j) - Fyp * B2p - Fy * B2pp


                aux.vars.G    = bulk.G[a,i,j]
                aux.vars.Gp   = Gp
                aux.vars.Gt   = Dx(bulk.G, a,i,j) - Fx * Gp
                aux.vars.Gh   = Dy(bulk.G, a,i,j) - Fy * Gp

                aux.vars.Gtt  = Dxx(bulk.G, a,i,j) - Dx(bulk.Fx, a,i,j) * Gp +
                    2.0 * Fx * u2 * Dx(Du_G, a,i,j) + Fx * Fxp * Gp + Fx * Fx * Gpp
                aux.vars.Ghh  = Dyy(bulk.G, a,i,j) - Dy(bulk.Fy, a,i,j) * Gp +
                    2.0 * Fy * u2 * Dy(Du_G, a,i,j) + Fy * Fyp * Gp + Fy * Fy * Gpp

                aux.vars.Gtp  = -u2 * Dx(Du_G, a,i,j) - Fxp * Gp - Fx * Gpp
                aux.vars.Ghp  = -u2 * Dy(Du_G, a,i,j) - Fyp * Gp - Fy * Gpp


                aux.vars.phi    = bulk.phi[a,i,j]
                aux.vars.phip   = phip
                aux.vars.phit   = Dx(bulk.phi, a,i,j) - Fx * phip
                aux.vars.phih   = Dy(bulk.phi, a,i,j) - Fy * phip

                aux.vars.phitt  = Dxx(bulk.phi, a,i,j) - Dx(bulk.Fx, a,i,j) * phip +
                    2.0 * Fx * u2 * Dx(Du_phi, a,i,j) + Fx * Fxp * phip + Fx * Fx * phipp
                aux.vars.phihh  = Dyy(bulk.phi, a,i,j) - Dy(bulk.Fy, a,i,j) * phip +
                    2.0 * Fy * u2 * Dy(Du_phi, a,i,j) + Fy * Fyp * phip + Fy * Fy * phipp

                aux.vars.phitp  = -u2 * Dx(Du_phi, a,i,j) - Fxp * phip - Fx * phipp
                aux.vars.phihp  = -u2 * Dy(Du_phi, a,i,j) - Fyp * phip - Fy * phipp


                aux.vars.S    = bulk.S[a,i,j]
                aux.vars.Sp   = Sp
                aux.vars.St   = Dx(bulk.S, a,i,j) - Fx * Sp
                aux.vars.Sh   = Dy(bulk.S, a,i,j) - Fy * Sp

                aux.vars.Stt  = Dxx(bulk.S, a,i,j) - Dx(bulk.Fx, a,i,j) * Sp +
                    2.0 * Fx * u2 * Dx(Du_S, a,i,j) + Fx * Fxp * Sp + Fx * Fx * Spp
                aux.vars.Shh  = Dyy(bulk.S, a,i,j) - Dy(bulk.Fy, a,i,j) * Sp +
                    2.0 * Fy * u2 * Dy(Du_S, a,i,j) + Fy * Fyp * Sp + Fy * Fy * Spp

                aux.vars.Stp  = -u2 * Dx(Du_S, a,i,j) - Fxp * Sp - Fx * Spp
                aux.vars.Shp  = -u2 * Dy(Du_S, a,i,j) - Fyp * Sp - Fy * Spp


                aux.vars.Fx    = bulk.Fx[a,i,j]
                aux.vars.Fxp   = Fxp
                aux.vars.Fxt   = Dx(bulk.Fx, a,i,j) - Fx * Fxp
                aux.vars.Fxh   = Dy(bulk.Fx, a,i,j) - Fy * Fxp

                aux.vars.Fxtp  = -u2 * Dx(Du_Fx, a,i,j) - Fxp * Fxp - Fx * Fxpp
                aux.vars.Fxhp  = -u2 * Dy(Du_Fx, a,i,j) - Fyp * Fxp - Fy * Fxpp


                aux.vars.Fy    = bulk.Fy[a,i,j]
                aux.vars.Fyp   = Fyp
                aux.vars.Fyt   = Dx(bulk.Fy, a,i,j) - Fy * Fyp
                aux.vars.Fyh   = Dy(bulk.Fy, a,i,j) - Fy * Fyp

                aux.vars.Fytp  = -u2 * Dx(Du_Fy, a,i,j) - Fyp * Fyp - Fy * Fypp
                aux.vars.Fyhp  = -u2 * Dy(Du_Fy, a,i,j) - Fyp * Fyp - Fy * Fypp

                aux.vars.Sd    = bulk.Sd[a,i,j]

                aux.vars.B2th = Dx(Dy, bulk.B2, a,i,j) - Dy(bulk.Fx, a,i,j) * B2p +
                    Fx * u2 * Dy(Du_B2, a,i,j) + Fy * u2 * Dx(Du_B2, a,i,j) +
                    Fy * Fxp * B2p + Fx * Fy * B2pp

                aux.vars.Gth = Dx(Dy, bulk.G, a,i,j) - Dy(bulk.Fx, a,i,j) * Gp +
                    Fx * u2 * Dy(Du_G, a,i,j) + Fy * u2 * Dx(Du_G, a,i,j) +
                    Fy * Fxp * Gp + Fx * Fy * Gpp

                aux.vars.Sth = Dx(Dy, bulk.S, a,i,j) - Dy(bulk.Fx, a,i,j) * Sp +
                    Fx * u2 * Dy(Du_S, a,i,j) + Fy * u2 * Dx(Du_S, a,i,j) +
                    Fy * Fxp * Sp + Fx * Fy * Spp


                B1dGd_outer_eq_coeff!(aux.AA, aux.BB, aux.CC, aux.SS, aux.vars)

                aux.b_vec2[a]    = -aux.SS[1]
                aux.b_vec2[a+Nu] = -aux.SS[2]
                @inbounds @simd for aa in eachindex(uu)
                    aux.A_mat2[a,aa]         = aux.AA[1,1] * Duu[a,aa] + aux.BB[1,1] * Du[a,aa]
                    aux.A_mat2[a,aa+Nu]      = aux.AA[1,2] * Duu[a,aa] + aux.BB[1,2] * Du[a,aa]
                    aux.A_mat2[a+Nu,aa]      = aux.AA[2,1] * Duu[a,aa] + aux.BB[2,1] * Du[a,aa]
                    aux.A_mat2[a+Nu,aa+Nu]   = aux.AA[2,2] * Duu[a,aa] + aux.BB[2,2] * Du[a,aa]
                end
                aux.A_mat2[a,a]       += aux.CC[1,1]
                aux.A_mat2[a,a+Nu]    += aux.CC[1,2]
                aux.A_mat2[a+Nu,a]    += aux.CC[2,1]
                aux.A_mat2[a+Nu,a+Nu] += aux.CC[2,2]
            end

            # BC (first order system)

            aux.b_vec2[1]     = BC.B1d[i,j]
            aux.A_mat2[1,:]  .= 0.0
            aux.A_mat2[1,1]   = 1.0

            aux.b_vec2[1+Nu]      = BC.Gd[i,j]
            aux.A_mat2[1+Nu,:]   .= 0.0
            aux.A_mat2[1+Nu,1+Nu] = 1.0

            solve_lin_system!(aux.sol2, aux.A_mat2, aux.b_vec2)

            @inbounds @simd for aa in eachindex(uu)
                bulk.B1d[aa,i,j] = aux.sol2[aa]
                bulk.Gd[aa,i,j]  = aux.sol2[aa+Nu]
            end

        end
    end

    nothing
end

function solve_phid_outer!(bulk::BulkVars, BC::BulkVars, nested::Nested)
    sys  = nested.sys
    uu   = nested.uu
    xx   = nested.xx
    yy   = nested.yy

    Du_B1   = nested.Du_B1
    Du_B2   = nested.Du_B2
    Du_G    = nested.Du_G
    Du_phi  = nested.Du_phi
    Du_S    = nested.Du_S
    Du_Fx   = nested.Du_Fx
    Du_Fy   = nested.Du_Fy
    Duu_B1  = nested.Duu_B1
    Duu_B2  = nested.Duu_B2
    Duu_G   = nested.Duu_G
    Duu_phi = nested.Duu_phi
    Duu_S   = nested.Duu_S
    Duu_Fx  = nested.Duu_Fx
    Duu_Fy  = nested.Duu_Fy

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

                B1p        = -u2 * Du_B1[a,i,j]
                B2p        = -u2 * Du_B2[a,i,j]
                Gp         = -u2 * Du_G[a,i,j]
                phip       = -u2 * Du_phi[a,i,j]
                Sp         = -u2 * Du_S[a,i,j]
                Fxp        = -u2 * Du_Fx[a,i,j]
                Fyp        = -u2 * Du_Fy[a,i,j]
                phip       = -u2 * Du_phi[a,i,j]

                B1pp       = 2*u3 * Du_B1[a,i,j]  + u4 * Duu_B1[a,i,j]
                B2pp       = 2*u3 * Du_B2[a,i,j]  + u4 * Duu_B2[a,i,j]
                Gpp        = 2*u3 * Du_G[a,i,j]   + u4 * Duu_G[a,i,j]
                phipp      = 2*u3 * Du_phi[a,i,j] + u4 * Duu_phi[a,i,j]
                Spp        = 2*u3 * Du_S[a,i,j]   + u4 * Duu_S[a,i,j]
                Fxpp       = 2*u3 * Du_Fx[a,i,j]  + u4 * Duu_Fx[a,i,j]
                Fypp       = 2*u3 * Du_Fy[a,i,j]  + u4 * Duu_Fy[a,i,j]
                phipp      = 2*u3 * Du_phi[a,i,j] + u4 * Duu_phi[a,i,j]

                Fx         = bulk.Fx[a,i,j]
                Fy         = bulk.Fy[a,i,j]

                aux.vars.u     = u

                aux.vars.B1    = bulk.B1[a,i,j]
                aux.vars.B1p   = B1p
                aux.vars.B1t   = Dx(bulk.B1, a,i,j) - Fx * B1p
                aux.vars.B1h   = Dy(bulk.B1, a,i,j) - Fy * B1p

                aux.vars.B1tt  = Dxx(bulk.B1, a,i,j) - Dx(bulk.Fx, a,i,j) * B1p +
                    2.0 * Fx * u2 * Dx(Du_B1, a,i,j) + Fx * Fxp * B1p + Fx * Fx * B1pp
                aux.vars.B1hh  = Dyy(bulk.B1, a,i,j) - Dy(bulk.Fy, a,i,j) * B1p +
                    2.0 * Fy * u2 * Dy(Du_B1, a,i,j) + Fy * Fyp * B1p + Fy * Fy * B1pp

                aux.vars.B1tp  = -u2 * Dx(Du_B1, a,i,j) - Fxp * B1p - Fx * B1pp
                aux.vars.B1hp  = -u2 * Dy(Du_B1, a,i,j) - Fyp * B1p - Fy * B1pp


                aux.vars.B2    = bulk.B2[a,i,j]
                aux.vars.B2p   = B2p
                aux.vars.B2t   = Dx(bulk.B2, a,i,j) - Fx * B2p
                aux.vars.B2h   = Dy(bulk.B2, a,i,j) - Fy * B2p

                aux.vars.B2tt  = Dxx(bulk.B2, a,i,j) - Dx(bulk.Fx, a,i,j) * B2p +
                    2.0 * Fx * u2 * Dx(Du_B2, a,i,j) + Fx * Fxp * B2p + Fx * Fx * B2pp
                aux.vars.B2hh  = Dyy(bulk.B2, a,i,j) - Dy(bulk.Fy, a,i,j) * B2p +
                    2.0 * Fy * u2 * Dy(Du_B2, a,i,j) + Fy * Fyp * B2p + Fy * Fy * B2pp

                aux.vars.B2tp  = -u2 * Dx(Du_B2, a,i,j) - Fxp * B2p - Fx * B2pp
                aux.vars.B2hp  = -u2 * Dy(Du_B2, a,i,j) - Fyp * B2p - Fy * B2pp


                aux.vars.G    = bulk.G[a,i,j]
                aux.vars.Gp   = Gp
                aux.vars.Gt   = Dx(bulk.G, a,i,j) - Fx * Gp
                aux.vars.Gh   = Dy(bulk.G, a,i,j) - Fy * Gp

                aux.vars.Gtt  = Dxx(bulk.G, a,i,j) - Dx(bulk.Fx, a,i,j) * Gp +
                    2.0 * Fx * u2 * Dx(Du_G, a,i,j) + Fx * Fxp * Gp + Fx * Fx * Gpp
                aux.vars.Ghh  = Dyy(bulk.G, a,i,j) - Dy(bulk.Fy, a,i,j) * Gp +
                    2.0 * Fy * u2 * Dy(Du_G, a,i,j) + Fy * Fyp * Gp + Fy * Fy * Gpp

                aux.vars.Gtp  = -u2 * Dx(Du_G, a,i,j) - Fxp * Gp - Fx * Gpp
                aux.vars.Ghp  = -u2 * Dy(Du_G, a,i,j) - Fyp * Gp - Fy * Gpp


                aux.vars.phi    = bulk.phi[a,i,j]
                aux.vars.phip   = phip
                aux.vars.phit   = Dx(bulk.phi, a,i,j) - Fx * phip
                aux.vars.phih   = Dy(bulk.phi, a,i,j) - Fy * phip

                aux.vars.phitt  = Dxx(bulk.phi, a,i,j) - Dx(bulk.Fx, a,i,j) * phip +
                    2.0 * Fx * u2 * Dx(Du_phi, a,i,j) + Fx * Fxp * phip + Fx * Fx * phipp
                aux.vars.phihh  = Dyy(bulk.phi, a,i,j) - Dy(bulk.Fy, a,i,j) * phip +
                    2.0 * Fy * u2 * Dy(Du_phi, a,i,j) + Fy * Fyp * phip + Fy * Fy * phipp

                aux.vars.phitp  = -u2 * Dx(Du_phi, a,i,j) - Fxp * phip - Fx * phipp
                aux.vars.phihp  = -u2 * Dy(Du_phi, a,i,j) - Fyp * phip - Fy * phipp


                aux.vars.S    = bulk.S[a,i,j]
                aux.vars.Sp   = Sp
                aux.vars.St   = Dx(bulk.S, a,i,j) - Fx * Sp
                aux.vars.Sh   = Dy(bulk.S, a,i,j) - Fy * Sp

                aux.vars.Stt  = Dxx(bulk.S, a,i,j) - Dx(bulk.Fx, a,i,j) * Sp +
                    2.0 * Fx * u2 * Dx(Du_S, a,i,j) + Fx * Fxp * Sp + Fx * Fx * Spp
                aux.vars.Shh  = Dyy(bulk.S, a,i,j) - Dy(bulk.Fy, a,i,j) * Sp +
                    2.0 * Fy * u2 * Dy(Du_S, a,i,j) + Fy * Fyp * Sp + Fy * Fy * Spp

                aux.vars.Stp  = -u2 * Dx(Du_S, a,i,j) - Fxp * Sp - Fx * Spp
                aux.vars.Shp  = -u2 * Dy(Du_S, a,i,j) - Fyp * Sp - Fy * Spp


                aux.vars.Fx    = bulk.Fx[a,i,j]
                aux.vars.Fxp   = Fxp
                aux.vars.Fxt   = Dx(bulk.Fx, a,i,j) - Fx * Fxp
                aux.vars.Fxh   = Dy(bulk.Fx, a,i,j) - Fy * Fxp

                aux.vars.Fxtp  = -u2 * Dx(Du_Fx, a,i,j) - Fxp * Fxp - Fx * Fxpp
                aux.vars.Fxhp  = -u2 * Dy(Du_Fx, a,i,j) - Fyp * Fxp - Fy * Fxpp


                aux.vars.Fy    = bulk.Fy[a,i,j]
                aux.vars.Fyp   = Fyp
                aux.vars.Fyt   = Dx(bulk.Fy, a,i,j) - Fy * Fyp
                aux.vars.Fyh   = Dy(bulk.Fy, a,i,j) - Fy * Fyp

                aux.vars.Fytp  = -u2 * Dx(Du_Fy, a,i,j) - Fyp * Fyp - Fy * Fypp
                aux.vars.Fyhp  = -u2 * Dy(Du_Fy, a,i,j) - Fyp * Fyp - Fy * Fypp

                aux.vars.Sd    = bulk.Sd[a,i,j]

                aux.vars.B2th = Dx(Dy, bulk.B2, a,i,j) - Dy(bulk.Fx, a,i,j) * B2p +
                    Fx * u2 * Dy(Du_B2, a,i,j) + Fy * u2 * Dx(Du_B2, a,i,j) +
                    Fy * Fxp * B2p + Fx * Fy * B2pp

                aux.vars.Gth = Dx(Dy, bulk.G, a,i,j) - Dy(bulk.Fx, a,i,j) * Gp +
                    Fx * u2 * Dy(Du_G, a,i,j) + Fy * u2 * Dx(Du_G, a,i,j) +
                    Fy * Fxp * Gp + Fx * Fy * Gpp

                aux.vars.Sth = Dx(Dy, bulk.S, a,i,j) - Dy(bulk.Fx, a,i,j) * Sp +
                    Fx * u2 * Dy(Du_S, a,i,j) + Fy * u2 * Dx(Du_S, a,i,j) +
                    Fy * Fxp * Sp + Fx * Fy * Spp

                aux.vars.phith = Dx(Dy, bulk.phi, a,i,j) - Dy(bulk.Fx, a,i,j) * phip +
                    Fx * u2 * Dy(Du_phi, a,i,j) + Fy * u2 * Dx(Du_phi, a,i,j) +
                    Fy * Fxp * phip + Fx * Fy * phipp


                phid_outer_eq_coeff!(aux.ABCS, aux.vars)

                aux.b_vec[a]   = -aux.ABCS[4]
                @inbounds @simd for aa in eachindex(uu)
                    aux.A_mat[a,aa] = aux.ABCS[1] * Duu[a,aa] + aux.ABCS[2] * Du[a,aa]
                end
                aux.A_mat[a,a] += aux.ABCS[3]
            end

            # BC (first order equation)

            aux.b_vec[1]    = BC.phid[i,j]
            aux.A_mat[1,:] .= 0.0
            aux.A_mat[1,1]  = 1.0

            sol = view(bulk.phid, :, i, j)
            solve_lin_system!(sol, aux.A_mat, aux.b_vec)
        end
    end

    nothing
end

function solve_A_outer!(bulk::BulkVars, BC::BulkVars, dBC::BulkVars, nested::Nested)
    sys  = nested.sys
    uu   = nested.uu
    xx   = nested.xx
    yy   = nested.yy

    Du_B1   = nested.Du_B1
    Du_B2   = nested.Du_B2
    Du_G    = nested.Du_G
    Du_phi  = nested.Du_phi
    Du_S    = nested.Du_S
    Du_Fx   = nested.Du_Fx
    Du_Fy   = nested.Du_Fy
    Duu_B1  = nested.Duu_B1
    Duu_B2  = nested.Duu_B2
    Duu_G   = nested.Duu_G
    Duu_phi = nested.Duu_phi
    Duu_S   = nested.Duu_S
    Duu_Fx  = nested.Duu_Fx
    Duu_Fy  = nested.Duu_Fy

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

                B1p        = -u2 * Du_B1[a,i,j]
                B2p        = -u2 * Du_B2[a,i,j]
                Gp         = -u2 * Du_G[a,i,j]
                phip       = -u2 * Du_phi[a,i,j]
                Sp         = -u2 * Du_S[a,i,j]
                Fxp        = -u2 * Du_Fx[a,i,j]
                Fyp        = -u2 * Du_Fy[a,i,j]
                phip       = -u2 * Du_phi[a,i,j]

                B1pp       = 2*u3 * Du_B1[a,i,j]  + u4 * Duu_B1[a,i,j]
                B2pp       = 2*u3 * Du_B2[a,i,j]  + u4 * Duu_B2[a,i,j]
                Gpp        = 2*u3 * Du_G[a,i,j]   + u4 * Duu_G[a,i,j]
                phipp      = 2*u3 * Du_phi[a,i,j] + u4 * Duu_phi[a,i,j]
                Spp        = 2*u3 * Du_S[a,i,j]   + u4 * Duu_S[a,i,j]
                Fxpp       = 2*u3 * Du_Fx[a,i,j]  + u4 * Duu_Fx[a,i,j]
                Fypp       = 2*u3 * Du_Fy[a,i,j]  + u4 * Duu_Fy[a,i,j]
                phipp      = 2*u3 * Du_phi[a,i,j] + u4 * Duu_phi[a,i,j]

                Fx         = bulk.Fx[a,i,j]
                Fy         = bulk.Fy[a,i,j]

                aux.vars.u     = u

                aux.vars.B1    = bulk.B1[a,i,j]
                aux.vars.B1p   = B1p
                aux.vars.B1t   = Dx(bulk.B1, a,i,j) - Fx * B1p
                aux.vars.B1h   = Dy(bulk.B1, a,i,j) - Fy * B1p

                aux.vars.B1tt  = Dxx(bulk.B1, a,i,j) - Dx(bulk.Fx, a,i,j) * B1p +
                    2.0 * Fx * u2 * Dx(Du_B1, a,i,j) + Fx * Fxp * B1p + Fx * Fx * B1pp
                aux.vars.B1hh  = Dyy(bulk.B1, a,i,j) - Dy(bulk.Fy, a,i,j) * B1p +
                    2.0 * Fy * u2 * Dy(Du_B1, a,i,j) + Fy * Fyp * B1p + Fy * Fy * B1pp

                aux.vars.B1tp  = -u2 * Dx(Du_B1, a,i,j) - Fxp * B1p - Fx * B1pp
                aux.vars.B1hp  = -u2 * Dy(Du_B1, a,i,j) - Fyp * B1p - Fy * B1pp


                aux.vars.B2    = bulk.B2[a,i,j]
                aux.vars.B2p   = B2p
                aux.vars.B2t   = Dx(bulk.B2, a,i,j) - Fx * B2p
                aux.vars.B2h   = Dy(bulk.B2, a,i,j) - Fy * B2p

                aux.vars.B2tt  = Dxx(bulk.B2, a,i,j) - Dx(bulk.Fx, a,i,j) * B2p +
                    2.0 * Fx * u2 * Dx(Du_B2, a,i,j) + Fx * Fxp * B2p + Fx * Fx * B2pp
                aux.vars.B2hh  = Dyy(bulk.B2, a,i,j) - Dy(bulk.Fy, a,i,j) * B2p +
                    2.0 * Fy * u2 * Dy(Du_B2, a,i,j) + Fy * Fyp * B2p + Fy * Fy * B2pp

                aux.vars.B2tp  = -u2 * Dx(Du_B2, a,i,j) - Fxp * B2p - Fx * B2pp
                aux.vars.B2hp  = -u2 * Dy(Du_B2, a,i,j) - Fyp * B2p - Fy * B2pp


                aux.vars.G    = bulk.G[a,i,j]
                aux.vars.Gp   = Gp
                aux.vars.Gt   = Dx(bulk.G, a,i,j) - Fx * Gp
                aux.vars.Gh   = Dy(bulk.G, a,i,j) - Fy * Gp

                aux.vars.Gtt  = Dxx(bulk.G, a,i,j) - Dx(bulk.Fx, a,i,j) * Gp +
                    2.0 * Fx * u2 * Dx(Du_G, a,i,j) + Fx * Fxp * Gp + Fx * Fx * Gpp
                aux.vars.Ghh  = Dyy(bulk.G, a,i,j) - Dy(bulk.Fy, a,i,j) * Gp +
                    2.0 * Fy * u2 * Dy(Du_G, a,i,j) + Fy * Fyp * Gp + Fy * Fy * Gpp

                aux.vars.Gtp  = -u2 * Dx(Du_G, a,i,j) - Fxp * Gp - Fx * Gpp
                aux.vars.Ghp  = -u2 * Dy(Du_G, a,i,j) - Fyp * Gp - Fy * Gpp


                aux.vars.phi    = bulk.phi[a,i,j]
                aux.vars.phip   = phip
                aux.vars.phit   = Dx(bulk.phi, a,i,j) - Fx * phip
                aux.vars.phih   = Dy(bulk.phi, a,i,j) - Fy * phip

                aux.vars.phitt  = Dxx(bulk.phi, a,i,j) - Dx(bulk.Fx, a,i,j) * phip +
                    2.0 * Fx * u2 * Dx(Du_phi, a,i,j) + Fx * Fxp * phip + Fx * Fx * phipp
                aux.vars.phihh  = Dyy(bulk.phi, a,i,j) - Dy(bulk.Fy, a,i,j) * phip +
                    2.0 * Fy * u2 * Dy(Du_phi, a,i,j) + Fy * Fyp * phip + Fy * Fy * phipp

                aux.vars.phitp  = -u2 * Dx(Du_phi, a,i,j) - Fxp * phip - Fx * phipp
                aux.vars.phihp  = -u2 * Dy(Du_phi, a,i,j) - Fyp * phip - Fy * phipp


                aux.vars.S    = bulk.S[a,i,j]
                aux.vars.Sp   = Sp
                aux.vars.St   = Dx(bulk.S, a,i,j) - Fx * Sp
                aux.vars.Sh   = Dy(bulk.S, a,i,j) - Fy * Sp

                aux.vars.Stt  = Dxx(bulk.S, a,i,j) - Dx(bulk.Fx, a,i,j) * Sp +
                    2.0 * Fx * u2 * Dx(Du_S, a,i,j) + Fx * Fxp * Sp + Fx * Fx * Spp
                aux.vars.Shh  = Dyy(bulk.S, a,i,j) - Dy(bulk.Fy, a,i,j) * Sp +
                    2.0 * Fy * u2 * Dy(Du_S, a,i,j) + Fy * Fyp * Sp + Fy * Fy * Spp

                aux.vars.Stp  = -u2 * Dx(Du_S, a,i,j) - Fxp * Sp - Fx * Spp
                aux.vars.Shp  = -u2 * Dy(Du_S, a,i,j) - Fyp * Sp - Fy * Spp


                aux.vars.Fx    = bulk.Fx[a,i,j]
                aux.vars.Fxp   = Fxp
                aux.vars.Fxt   = Dx(bulk.Fx, a,i,j) - Fx * Fxp
                aux.vars.Fxh   = Dy(bulk.Fx, a,i,j) - Fy * Fxp

                aux.vars.Fxtp  = -u2 * Dx(Du_Fx, a,i,j) - Fxp * Fxp - Fx * Fxpp
                aux.vars.Fxhp  = -u2 * Dy(Du_Fx, a,i,j) - Fyp * Fxp - Fy * Fxpp


                aux.vars.Fy    = bulk.Fy[a,i,j]
                aux.vars.Fyp   = Fyp
                aux.vars.Fyt   = Dx(bulk.Fy, a,i,j) - Fy * Fyp
                aux.vars.Fyh   = Dy(bulk.Fy, a,i,j) - Fy * Fyp

                aux.vars.Fytp  = -u2 * Dx(Du_Fy, a,i,j) - Fyp * Fyp - Fy * Fypp
                aux.vars.Fyhp  = -u2 * Dy(Du_Fy, a,i,j) - Fyp * Fyp - Fy * Fypp

                aux.vars.Sd    = bulk.Sd[a,i,j]
                aux.vars.B1d   = bulk.B1d[a,i,j]
                aux.vars.B2d   = bulk.B2d[a,i,j]
                aux.vars.Gd    = bulk.Gd[a,i,j]
                aux.vars.phid  = bulk.phid[a,i,j]

                aux.vars.B2th = Dx(Dy, bulk.B2, a,i,j) - Dy(bulk.Fx, a,i,j) * B2p +
                    Fx * u2 * Dy(Du_B2, a,i,j) + Fy * u2 * Dx(Du_B2, a,i,j) +
                    Fy * Fxp * B2p + Fx * Fy * B2pp

                aux.vars.Gth = Dx(Dy, bulk.G, a,i,j) - Dy(bulk.Fx, a,i,j) * Gp +
                    Fx * u2 * Dy(Du_G, a,i,j) + Fy * u2 * Dx(Du_G, a,i,j) +
                    Fy * Fxp * Gp + Fx * Fy * Gpp

                aux.vars.Sth = Dx(Dy, bulk.S, a,i,j) - Dy(bulk.Fx, a,i,j) * Sp +
                    Fx * u2 * Dy(Du_S, a,i,j) + Fy * u2 * Dx(Du_S, a,i,j) +
                    Fy * Fxp * Sp + Fx * Fy * Spp

                aux.vars.phith = Dx(Dy, bulk.phi, a,i,j) - Dy(bulk.Fx, a,i,j) * phip +
                    Fx * u2 * Dy(Du_phi, a,i,j) + Fy * u2 * Dx(Du_phi, a,i,j) +
                    Fy * Fxp * phip + Fx * Fy * phipp


                A_outer_eq_coeff!(aux.ABCS, aux.vars)

                aux.b_vec[a]   = -aux.ABCS[4]
                @inbounds @simd for aa in eachindex(uu)
                    aux.A_mat[a,aa] = aux.ABCS[1] * Duu[a,aa] + aux.ABCS[2] * Du[a,aa]
                end
                aux.A_mat[a,a] += aux.ABCS[3]
            end

            # BC

            aux.b_vec[1]    = BC.A[i,j]
            aux.A_mat[1,:] .= 0.0
            aux.A_mat[1,1]  = 1.0

            aux.b_vec[end]    = dBC.A[i,j]
            @inbounds @simd for aa in eachindex(uu)
                aux.A_mat[end,aa]  = Du[1,aa]
            end

            sol = view(bulk.A, :, i, j)
            solve_lin_system!(sol, aux.A_mat, aux.b_vec)
        end
    end

    nothing
end

function solve_nested_outer!(bulk::BulkVars, BC::BulkVars, dBC::BulkVars, nested::Nested)
    sys  = nested.sys

    Du_B1   = nested.Du_B1
    Du_B2   = nested.Du_B2
    Du_G    = nested.Du_G
    Du_phi  = nested.Du_phi
    Du_S    = nested.Du_S
    Du_Fx   = nested.Du_Fx
    Du_Fy   = nested.Du_Fy
    Duu_B1  = nested.Duu_B1
    Duu_B2  = nested.Duu_B2
    Duu_G   = nested.Duu_G
    Duu_phi = nested.Duu_phi
    Duu_S   = nested.Duu_S
    Duu_Fx  = nested.Duu_Fx
    Duu_Fy  = nested.Duu_Fy

    Du  = sys.Du
    Duu = sys.Duu
    Dx  = sys.Dx
    Dxx = sys.Dxx
    Dy  = sys.Dy
    Dyy = sys.Dyy

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
    solve_S_outer!(bulk, BC, dBC, nested)

    # take u-derivatives of S
    @sync begin
        @spawn mul!(Du_S,   Du,  bulk.S)
        @spawn mul!(Duu_S,  Duu, bulk.S)
    end

    # solve for Fx and Fy
    solve_Fxy_outer!(bulk, BC, dBC, nested)

    # take u-derivatives of Fx and Fy
    @sync begin
        @spawn mul!(Du_Fx,   Du,  bulk.Fx)
        @spawn mul!(Du_Fy,   Du,  bulk.Fy)
        @spawn mul!(Duu_Fx,  Duu, bulk.Fx)
        @spawn mul!(Duu_Fy,  Duu, bulk.Fy)
    end

    # solve for Sd
    solve_Sd_outer!(bulk, BC, nested)

    # solving for B2d, (B1d,Gd) and phid are independent processes. we can
    # therefore @spawn, here
    @sync begin
        @spawn solve_B2d_outer!(bulk, BC, nested)
        @spawn solve_B1dGd_outer!(bulk, BC, nested)
        @spawn solve_phid_outer!(bulk, BC, nested)
    end

    # solve for A
    solve_A_outer!(bulk, BC, dBC, nested)

    nothing
end


function syncBCs!(BC::BulkVars, dBC::BulkVars, bulk::BulkVars, nested::Nested)
    Du = nested.sys.Du

    Nu, Nx, Ny = size(bulk.S)

    @fastmath @inbounds @threads for j in 1:Ny
        @inbounds @simd for i in 1:Nx
            BC.S[i,j]    = bulk.S[end,i,j]
            BC.Fx[i,j]   = bulk.Fx[end,i,j]
            BC.Fy[i,j]   = bulk.Fy[end,i,j]
            BC.Sd[i,j]   = bulk.Sd[end,i,j]
            BC.B1d[i,j]  = bulk.B1d[end,i,j]
            BC.B2d[i,j]  = bulk.B2d[end,i,j]
            BC.Gd[i,j]   = bulk.Gd[end,i,j]
            BC.phid[i,j] = bulk.phid[end,i,j]
            BC.A[i,j]    = bulk.A[end,i,j]

            dBC.S[i,j]   = nested.Du_S[end,i,j]
            dBC.Fx[i,j]  = nested.Du_Fx[end,i,j]
            dBC.Fy[i,j]  = nested.Du_Fy[end,i,j]
            dBC.A[i,j]   = Du(bulk.A, Nu,i,j)
        end
    end

    nothing
end

function solve_nested!(bulks::Vector, BCs::Vector, dBCs::Vector,
                       nesteds::Vector)
    Nsys = length(nesteds)

    # We assume that the first entry on these arrays is the inner grid, and that
    # there is only one domain spanning this grid. If we ever change this
    # construction we must remember to make the appropriate changes here.
    for i in 2:Nsys-1
        solve_nested_outer!(bulks[i], BCs[i], dBCs[i], nesteds[i])
        syncBCs!(BCs[i+1], dBCs[i+1], bulks[i], nesteds[i])
    end
    solve_nested_outer!(bulks[Nsys], BCs[Nsys], dBCs[Nsys], nesteds[Nsys])

    nothing
end

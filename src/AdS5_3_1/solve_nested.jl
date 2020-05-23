
import Base.Threads.@threads
import Base.Threads.@spawn
using LinearAlgebra

function solve_lin_system!(A_mat, b_vec)
    # passing Val(false) to the second argument turns off pivoting. it seems to
    # improve speed for the small matrices that we typically consider. we can
    # revisit this (or make it a parameter) if needed.
    A_fact = lu!(A_mat, Val(false))
    ldiv!(A_fact, b_vec)        # b_vec is overwritten to store the result
    nothing
end

struct Aux{T<:Real}
    A_mat   :: Matrix{T}
    b_vec   :: Vector{T}
    ABCS    :: Vector{T}
    A_mat2  :: Matrix{T}
    b_vec2  :: Vector{T}
    AA      :: Matrix{T}
    BB      :: Matrix{T}
    CC      :: Matrix{T}
    SS      :: Vector{T}
    function Aux{T}(N::Int) where {T<:Real}
        A_mat   = zeros(T, N, N)
        b_vec   = zeros(T, N)
        ABCS    = zeros(T, 4)
        A_mat2  = zeros(T, 2*N, 2*N)
        b_vec2  = zeros(T, 2*N)
        AA      = zeros(T, 2,2)
        BB      = zeros(T, 2,2)
        CC      = zeros(T, 2,2)
        SS      = zeros(T, 2)
        new(A_mat, b_vec, ABCS, A_mat2, b_vec2, AA, BB, CC, SS)
    end
end

struct Nested{S,T<:Real,D}
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
    Nu, Nx, Ny = size(sys)
    uu = sys.ucoord[:]
    xx = sys.xcoord[:]
    yy = sys.ycoord[:]
    T  = Jecco.coord_eltype(sys.ucoord)

    Du_B1    = zeros(T, Nu, Nx, Ny)
    Du_B2    = zeros(T, Nu, Nx, Ny)
    Du_G     = zeros(T, Nu, Nx, Ny)
    Du_phi   = zeros(T, Nu, Nx, Ny)
    Du_S     = zeros(T, Nu, Nx, Ny)
    Du_Fx    = zeros(T, Nu, Nx, Ny)
    Du_Fy    = zeros(T, Nu, Nx, Ny)
    Duu_B1   = zeros(T, Nu, Nx, Ny)
    Duu_B2   = zeros(T, Nu, Nx, Ny)
    Duu_G    = zeros(T, Nu, Nx, Ny)
    Duu_phi  = zeros(T, Nu, Nx, Ny)
    Duu_S    = zeros(T, Nu, Nx, Ny)
    Duu_Fx   = zeros(T, Nu, Nx, Ny)
    Duu_Fy   = zeros(T, Nu, Nx, Ny)

    nt = Threads.nthreads()
    # pre-allocate thread-local aux quantities
    aux_acc = [Aux{T}(Nu) for _ in 1:nt]

    Nested{typeof(sys),T,typeof(Du_B1)}(sys, uu, xx, yy, Du_B1,
                                        Du_B2, Du_G, Du_phi, Du_S, Du_Fx,
                                        Du_Fy, Duu_B1, Duu_B2, Duu_G, Duu_phi,
                                        Duu_S, Duu_Fx, Duu_Fy, aux_acc)
end

Nested(systems::SystemPartition) = [Nested(sys) for sys in systems]



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

function solve_S!(bulk::Bulk, BC::Bulk, dBC::Bulk, gauge::Gauge,
                  nested::Nested, evoleq::EvolEq)
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

    phi0  = evoleq.phi0

    @fastmath @inbounds @threads for j in eachindex(yy)
        @inbounds for i in eachindex(xx)
            id  = Threads.threadid()
            aux = aux_acc[id]

            xi  = gauge.xi[1,i,j]

            @inbounds @simd for a in eachindex(uu)
                u     = uu[a]

                B1    = bulk.B1[a,i,j]
                B2    = bulk.B2[a,i,j]
                G     = bulk.G[a,i,j]
                phi   = bulk.phi[a,i,j]

                B1p   = -u*u * Du_B1[a,i,j]
                B2p   = -u*u * Du_B2[a,i,j]
                Gp    = -u*u * Du_G[a,i,j]
                phip  = -u*u * Du_phi[a,i,j]

                vars = (phi0, u, xi, B1, B1p, B2, B2p, G, Gp, phi, phip)

                S_eq_coeff!(aux.ABCS, vars, sys.gridtype)

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

            solve_lin_system!(aux.A_mat, aux.b_vec)

            @inbounds @simd for aa in eachindex(uu)
                bulk.S[aa,i,j] = aux.b_vec[aa]
            end

        end
    end

    nothing
end

function solve_Fxy!(bulk::Bulk, BC::Bulk, dBC::Bulk, gauge::Gauge,
                    nested::Nested, evoleq::EvolEq)
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

    phi0  = evoleq.phi0

    Nu = length(uu)

    @fastmath @inbounds @threads for j in eachindex(yy)
        @inbounds for i in eachindex(xx)
            id  = Threads.threadid()
            aux = aux_acc[id]

            xi    = gauge.xi[1,i,j]
            xi_x  = Dx(gauge.xi, 1,i,j)
            xi_y  = Dy(gauge.xi, 1,i,j)

            @inbounds @simd for a in eachindex(uu)
                u     = uu[a]
                u2    = u * u
                u3    = u * u2
                u4    = u2 * u2

                B1    = bulk.B1[a,i,j]
                B1p   = -u2 * Du_B1[a,i,j]
                B1_x  = Dx(bulk.B1, a,i,j)
                B1_y  = Dy(bulk.B1, a,i,j)
                B1pp  = 2*u3 * Du_B1[a,i,j] + u4 * Duu_B1[a,i,j]
                B1p_x = -u2 * Dx(Du_B1, a,i,j)
                B1p_y = -u2 * Dy(Du_B1, a,i,j)

                B2    = bulk.B2[a,i,j]
                B2p   = -u2 * Du_B2[a,i,j]
                B2_x  = Dx(bulk.B2, a,i,j)
                B2_y  = Dy(bulk.B2, a,i,j)
                B2pp  = 2*u3 * Du_B2[a,i,j] + u4 * Duu_B2[a,i,j]
                B2p_x = -u2 * Dx(Du_B2, a,i,j)
                B2p_y = -u2 * Dy(Du_B2, a,i,j)

                G     = bulk.G[a,i,j]
                Gp    = -u2 * Du_G[a,i,j]
                G_x   = Dx(bulk.G, a,i,j)
                G_y   = Dy(bulk.G, a,i,j)
                Gpp   = 2*u3 * Du_G[a,i,j] + u4 * Duu_G[a,i,j]
                Gp_x  = -u2 * Dx(Du_G, a,i,j)
                Gp_y  = -u2 * Dy(Du_G, a,i,j)

                phi   = bulk.phi[a,i,j]
                phip  = -u2 * Du_phi[a,i,j]
                phi_x = Dx(bulk.phi, a,i,j)
                phi_y = Dy(bulk.phi, a,i,j)
                # phipp   = 2*u3 * Du_phi[a,i,j] + u4 * Duu_phi[a,i,j]
                # phip_x  = -u2 * Dx(Du_phi, a,i,j)
                # phip_y  = -u2 * Dy(Du_phi, a,i,j)

                S     = bulk.S[a,i,j]
                Sp    = -u2 * Du_S[a,i,j]
                S_x   = Dx(bulk.S, a,i,j)
                S_y   = Dy(bulk.S, a,i,j)
                Spp   = 2*u3 * Du_S[a,i,j] + u4 * Duu_S[a,i,j]
                Sp_x  = -u2 * Dx(Du_S, a,i,j)
                Sp_y  = -u2 * Dy(Du_S, a,i,j)

                vars = (
                    phi0, u, xi, xi_x, xi_y,
                    B1     ,    B2     ,    G      ,    phi    ,    S      ,
                    B1p    ,    B2p    ,    Gp     ,    phip   ,    Sp     ,
                    B1pp   ,    B2pp   ,    Gpp    ,                Spp    ,
                    B1_x   ,    B2_x   ,    G_x    ,    phi_x  ,    S_x    ,
                    B1_y   ,    B2_y   ,    G_y    ,    phi_y  ,    S_y    ,
                    B1p_x  ,    B2p_x  ,    Gp_x   ,                Sp_x   ,
                    B1p_y  ,    B2p_y  ,    Gp_y   ,                Sp_y
                )

                Fxy_eq_coeff!(aux.AA, aux.BB, aux.CC, aux.SS, vars, sys.gridtype)

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

            solve_lin_system!(aux.A_mat2, aux.b_vec2)

            @inbounds @simd for aa in eachindex(uu)
                bulk.Fx[aa,i,j] = aux.b_vec2[aa]
                bulk.Fy[aa,i,j] = aux.b_vec2[aa+Nu]
            end

        end
    end

    nothing
end

function solve_Sd!(bulk::Bulk, BC::Bulk, gauge::Gauge,
                   nested::Nested, evoleq::EvolEq)
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

    potential = evoleq.potential
    phi0      = evoleq.phi0

    @fastmath @inbounds @threads for j in eachindex(yy)
        @inbounds for i in eachindex(xx)
            id   = Threads.threadid()
            aux  = aux_acc[id]

            xi    = gauge.xi[1,i,j]
            xi_x  = Dx(gauge.xi, 1,i,j)
            xi_y  = Dy(gauge.xi, 1,i,j)
            xi_xx = Dxx(gauge.xi, 1,i,j)
            xi_yy = Dyy(gauge.xi, 1,i,j)
            xi_xy = Dx(Dy, gauge.xi, 1,i,j)

            @inbounds @simd for a in eachindex(uu)
                u     = uu[a]
                u2    = u * u
                u3    = u * u2
                u4    = u2 * u2

                B1    = bulk.B1[a,i,j]
                B2    = bulk.B2[a,i,j]
                G     = bulk.G[a,i,j]
                phi   = bulk.phi[a,i,j]
                S     = bulk.S[a,i,j]
                Fx    = bulk.Fx[a,i,j]
                Fy    = bulk.Fy[a,i,j]

                # r derivatives

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

                # x and y derivatives

                B1_x       = Dx(bulk.B1, a,i,j)
                B2_x       = Dx(bulk.B2, a,i,j)
                G_x        = Dx(bulk.G,  a,i,j)
                phi_x      = Dx(bulk.phi,a,i,j)
                S_x        = Dx(bulk.S,  a,i,j)
                Fx_x       = Dx(bulk.Fx, a,i,j)
                Fy_x       = Dx(bulk.Fy, a,i,j)

                B1_y       = Dy(bulk.B1, a,i,j)
                B2_y       = Dy(bulk.B2, a,i,j)
                G_y        = Dy(bulk.G,  a,i,j)
                phi_y      = Dy(bulk.phi,a,i,j)
                S_y        = Dy(bulk.S,  a,i,j)
                Fx_y       = Dy(bulk.Fx, a,i,j)
                Fy_y       = Dy(bulk.Fy, a,i,j)

                B1p_x      = -u2 * Dx(Du_B1, a,i,j)
                B2p_x      = -u2 * Dx(Du_B2, a,i,j)
                Gp_x       = -u2 * Dx(Du_G,  a,i,j)
                phip_x     = -u2 * Dx(Du_phi,a,i,j)
                Sp_x       = -u2 * Dx(Du_S,  a,i,j)
                Fxp_x      = -u2 * Dx(Du_Fx,  a,i,j)
                Fyp_x      = -u2 * Dx(Du_Fy,  a,i,j)

                B1p_y      = -u2 * Dy(Du_B1, a,i,j)
                B2p_y      = -u2 * Dy(Du_B2, a,i,j)
                Gp_y       = -u2 * Dy(Du_G,  a,i,j)
                phip_y     = -u2 * Dy(Du_phi,a,i,j)
                Sp_y       = -u2 * Dy(Du_S,  a,i,j)
                Fxp_y      = -u2 * Dy(Du_Fx, a,i,j)
                Fyp_y      = -u2 * Dy(Du_Fy, a,i,j)

                B1_xx      = Dxx(bulk.B1, a,i,j)
                B2_xx      = Dxx(bulk.B2, a,i,j)
                G_xx       = Dxx(bulk.G,  a,i,j)
                phi_xx     = Dxx(bulk.phi,a,i,j)
                S_xx       = Dxx(bulk.S,  a,i,j)

                B1_yy      = Dyy(bulk.B1, a,i,j)
                B2_yy      = Dyy(bulk.B2, a,i,j)
                G_yy       = Dyy(bulk.G,  a,i,j)
                phi_yy     = Dyy(bulk.phi,a,i,j)
                S_yy       = Dyy(bulk.S,  a,i,j)

                B2_xy      = Dx(Dy, bulk.B2, a,i,j)
                G_xy       = Dx(Dy, bulk.G,  a,i,j)
                S_xy       = Dx(Dy, bulk.S,  a,i,j)

                vars = (
                    potential, phi0, u, xi, xi_x, xi_y, xi_xx, xi_yy, xi_xy,
                    B1     ,    B2     ,    G      ,    phi    ,    S      ,    Fx     ,    Fy     ,
                    B1p    ,    B2p    ,    Gp     ,    phip   ,    Sp     ,    Fxp    ,    Fyp    ,
                    B1pp   ,    B2pp   ,    Gpp    ,    phipp  ,    Spp    ,    Fxpp   ,    Fypp   ,
                    B1_x   ,    B2_x   ,    G_x    ,    phi_x  ,    S_x    ,    Fx_x   ,    Fy_x   ,
                    B1_y   ,    B2_y   ,    G_y    ,    phi_y  ,    S_y    ,    Fx_y   ,    Fy_y   ,
                    B1p_x  ,    B2p_x  ,    Gp_x   ,    phip_x ,    Sp_x   ,    Fxp_x  ,    Fyp_x  ,
                    B1p_y  ,    B2p_y  ,    Gp_y   ,    phip_y ,    Sp_y   ,    Fxp_y  ,    Fyp_y  ,
                    B1_xx  ,    B2_xx  ,    G_xx   ,    phi_xx ,    S_xx   ,
                    B1_yy  ,    B2_yy  ,    G_yy   ,    phi_yy ,    S_yy   ,
                                B2_xy  ,    G_xy   ,                S_xy
                )

                Sd_eq_coeff!(aux.ABCS, vars, sys.gridtype)

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

            solve_lin_system!(aux.A_mat, aux.b_vec)

            @inbounds @simd for aa in eachindex(uu)
                bulk.Sd[aa,i,j] = aux.b_vec[aa]
            end

        end
    end

    nothing
end

function solve_B2d!(bulk::Bulk, BC::Bulk, gauge::Gauge,
                    nested::Nested, evoleq::EvolEq)
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

    phi0  = evoleq.phi0

    @fastmath @inbounds @threads for j in eachindex(yy)
        @inbounds for i in eachindex(xx)
            id  = Threads.threadid()
            aux = aux_acc[id]

            xi    = gauge.xi[1,i,j]
            xi_x  = Dx(gauge.xi, 1,i,j)
            xi_y  = Dy(gauge.xi, 1,i,j)
            xi_xx = Dxx(gauge.xi, 1,i,j)
            xi_yy = Dyy(gauge.xi, 1,i,j)
            xi_xy = Dx(Dy, gauge.xi, 1,i,j)

            @inbounds @simd for a in eachindex(uu)
                u     = uu[a]
                u2    = u * u
                u3    = u * u2
                u4    = u2 * u2

                B1    = bulk.B1[a,i,j]
                B2    = bulk.B2[a,i,j]
                G     = bulk.G[a,i,j]
                phi   = bulk.phi[a,i,j]
                S     = bulk.S[a,i,j]
                Fx    = bulk.Fx[a,i,j]
                Fy    = bulk.Fy[a,i,j]
                Sd    = bulk.Sd[a,i,j]

                # r derivatives

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

                # x and y derivatives

                B1_x       = Dx(bulk.B1, a,i,j)
                B2_x       = Dx(bulk.B2, a,i,j)
                G_x        = Dx(bulk.G,  a,i,j)
                phi_x      = Dx(bulk.phi,a,i,j)
                S_x        = Dx(bulk.S,  a,i,j)
                Fx_x       = Dx(bulk.Fx, a,i,j)
                Fy_x       = Dx(bulk.Fy, a,i,j)

                B1_y       = Dy(bulk.B1, a,i,j)
                B2_y       = Dy(bulk.B2, a,i,j)
                G_y        = Dy(bulk.G,  a,i,j)
                phi_y      = Dy(bulk.phi,a,i,j)
                S_y        = Dy(bulk.S,  a,i,j)
                Fx_y       = Dy(bulk.Fx, a,i,j)
                Fy_y       = Dy(bulk.Fy, a,i,j)

                B1p_x      = -u2 * Dx(Du_B1, a,i,j)
                B2p_x      = -u2 * Dx(Du_B2, a,i,j)
                Gp_x       = -u2 * Dx(Du_G,  a,i,j)
                phip_x     = -u2 * Dx(Du_phi,a,i,j)
                Sp_x       = -u2 * Dx(Du_S,  a,i,j)
                Fxp_x      = -u2 * Dx(Du_Fx,  a,i,j)
                Fyp_x      = -u2 * Dx(Du_Fy,  a,i,j)

                B1p_y      = -u2 * Dy(Du_B1, a,i,j)
                B2p_y      = -u2 * Dy(Du_B2, a,i,j)
                Gp_y       = -u2 * Dy(Du_G,  a,i,j)
                phip_y     = -u2 * Dy(Du_phi,a,i,j)
                Sp_y       = -u2 * Dy(Du_S,  a,i,j)
                Fxp_y      = -u2 * Dy(Du_Fx, a,i,j)
                Fyp_y      = -u2 * Dy(Du_Fy, a,i,j)

                B1_xx      = Dxx(bulk.B1, a,i,j)
                B2_xx      = Dxx(bulk.B2, a,i,j)
                G_xx       = Dxx(bulk.G,  a,i,j)
                phi_xx     = Dxx(bulk.phi,a,i,j)
                S_xx       = Dxx(bulk.S,  a,i,j)

                B1_yy      = Dyy(bulk.B1, a,i,j)
                B2_yy      = Dyy(bulk.B2, a,i,j)
                G_yy       = Dyy(bulk.G,  a,i,j)
                phi_yy     = Dyy(bulk.phi,a,i,j)
                S_yy       = Dyy(bulk.S,  a,i,j)

                B2_xy      = Dx(Dy, bulk.B2, a,i,j)
                G_xy       = Dx(Dy, bulk.G,  a,i,j)
                S_xy       = Dx(Dy, bulk.S,  a,i,j)

                vars = (
                    phi0, u, xi, xi_x, xi_y, xi_xx, xi_yy, xi_xy,
                    B1     ,    B2     ,    G      ,    phi    ,    S      ,    Fx     ,    Fy     ,  Sd,
                    B1p    ,    B2p    ,    Gp     ,    phip   ,    Sp     ,    Fxp    ,    Fyp    ,
                    B1pp   ,    B2pp   ,    Gpp    ,    phipp  ,    Spp    ,    Fxpp   ,    Fypp   ,
                    B1_x   ,    B2_x   ,    G_x    ,    phi_x  ,    S_x    ,    Fx_x   ,    Fy_x   ,
                    B1_y   ,    B2_y   ,    G_y    ,    phi_y  ,    S_y    ,    Fx_y   ,    Fy_y   ,
                    B1p_x  ,    B2p_x  ,    Gp_x   ,    phip_x ,    Sp_x   ,    Fxp_x  ,    Fyp_x  ,
                    B1p_y  ,    B2p_y  ,    Gp_y   ,    phip_y ,    Sp_y   ,    Fxp_y  ,    Fyp_y  ,
                    B1_xx  ,    B2_xx  ,    G_xx   ,    phi_xx ,    S_xx   ,
                    B1_yy  ,    B2_yy  ,    G_yy   ,    phi_yy ,    S_yy   ,
                                B2_xy  ,    G_xy   ,                S_xy
                )

                B2d_eq_coeff!(aux.ABCS, vars, sys.gridtype)

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

            solve_lin_system!(aux.A_mat, aux.b_vec)

            @inbounds @simd for aa in eachindex(uu)
                bulk.B2d[aa,i,j] = aux.b_vec[aa]
            end

        end
    end

    nothing
end

function solve_B1dGd!(bulk::Bulk, BC::Bulk, gauge::Gauge,
                      nested::Nested, evoleq::EvolEq)
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

    phi0  = evoleq.phi0

    Nu = length(uu)

    @fastmath @inbounds @threads for j in eachindex(yy)
        @inbounds for i in eachindex(xx)
            id  = Threads.threadid()
            aux = aux_acc[id]

            xi    = gauge.xi[1,i,j]
            xi_x  = Dx(gauge.xi, 1,i,j)
            xi_y  = Dy(gauge.xi, 1,i,j)
            xi_xx = Dxx(gauge.xi, 1,i,j)
            xi_yy = Dyy(gauge.xi, 1,i,j)
            xi_xy = Dx(Dy, gauge.xi, 1,i,j)

            @inbounds @simd for a in eachindex(uu)
                u     = uu[a]
                u2    = u * u
                u3    = u * u2
                u4    = u2 * u2

                B1    = bulk.B1[a,i,j]
                B2    = bulk.B2[a,i,j]
                G     = bulk.G[a,i,j]
                phi   = bulk.phi[a,i,j]
                S     = bulk.S[a,i,j]
                Fx    = bulk.Fx[a,i,j]
                Fy    = bulk.Fy[a,i,j]
                Sd    = bulk.Sd[a,i,j]

                # r derivatives

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

                # x and y derivatives

                B1_x       = Dx(bulk.B1, a,i,j)
                B2_x       = Dx(bulk.B2, a,i,j)
                G_x        = Dx(bulk.G,  a,i,j)
                phi_x      = Dx(bulk.phi,a,i,j)
                S_x        = Dx(bulk.S,  a,i,j)
                Fx_x       = Dx(bulk.Fx, a,i,j)
                Fy_x       = Dx(bulk.Fy, a,i,j)

                B1_y       = Dy(bulk.B1, a,i,j)
                B2_y       = Dy(bulk.B2, a,i,j)
                G_y        = Dy(bulk.G,  a,i,j)
                phi_y      = Dy(bulk.phi,a,i,j)
                S_y        = Dy(bulk.S,  a,i,j)
                Fx_y       = Dy(bulk.Fx, a,i,j)
                Fy_y       = Dy(bulk.Fy, a,i,j)

                B1p_x      = -u2 * Dx(Du_B1, a,i,j)
                B2p_x      = -u2 * Dx(Du_B2, a,i,j)
                Gp_x       = -u2 * Dx(Du_G,  a,i,j)
                phip_x     = -u2 * Dx(Du_phi,a,i,j)
                Sp_x       = -u2 * Dx(Du_S,  a,i,j)
                Fxp_x      = -u2 * Dx(Du_Fx,  a,i,j)
                Fyp_x      = -u2 * Dx(Du_Fy,  a,i,j)

                B1p_y      = -u2 * Dy(Du_B1, a,i,j)
                B2p_y      = -u2 * Dy(Du_B2, a,i,j)
                Gp_y       = -u2 * Dy(Du_G,  a,i,j)
                phip_y     = -u2 * Dy(Du_phi,a,i,j)
                Sp_y       = -u2 * Dy(Du_S,  a,i,j)
                Fxp_y      = -u2 * Dy(Du_Fx, a,i,j)
                Fyp_y      = -u2 * Dy(Du_Fy, a,i,j)

                B1_xx      = Dxx(bulk.B1, a,i,j)
                B2_xx      = Dxx(bulk.B2, a,i,j)
                G_xx       = Dxx(bulk.G,  a,i,j)
                phi_xx     = Dxx(bulk.phi,a,i,j)
                S_xx       = Dxx(bulk.S,  a,i,j)

                B1_yy      = Dyy(bulk.B1, a,i,j)
                B2_yy      = Dyy(bulk.B2, a,i,j)
                G_yy       = Dyy(bulk.G,  a,i,j)
                phi_yy     = Dyy(bulk.phi,a,i,j)
                S_yy       = Dyy(bulk.S,  a,i,j)

                B2_xy      = Dx(Dy, bulk.B2, a,i,j)
                G_xy       = Dx(Dy, bulk.G,  a,i,j)
                S_xy       = Dx(Dy, bulk.S,  a,i,j)

                vars = (
                    phi0, u, xi, xi_x, xi_y, xi_xx, xi_yy, xi_xy,
                    B1     ,    B2     ,    G      ,    phi    ,    S      ,    Fx     ,    Fy     ,  Sd,
                    B1p    ,    B2p    ,    Gp     ,    phip   ,    Sp     ,    Fxp    ,    Fyp    ,
                    B1pp   ,    B2pp   ,    Gpp    ,    phipp  ,    Spp    ,    Fxpp   ,    Fypp   ,
                    B1_x   ,    B2_x   ,    G_x    ,    phi_x  ,    S_x    ,    Fx_x   ,    Fy_x   ,
                    B1_y   ,    B2_y   ,    G_y    ,    phi_y  ,    S_y    ,    Fx_y   ,    Fy_y   ,
                    B1p_x  ,    B2p_x  ,    Gp_x   ,    phip_x ,    Sp_x   ,    Fxp_x  ,    Fyp_x  ,
                    B1p_y  ,    B2p_y  ,    Gp_y   ,    phip_y ,    Sp_y   ,    Fxp_y  ,    Fyp_y  ,
                    B1_xx  ,    B2_xx  ,    G_xx   ,    phi_xx ,    S_xx   ,
                    B1_yy  ,    B2_yy  ,    G_yy   ,    phi_yy ,    S_yy   ,
                                B2_xy  ,    G_xy   ,                S_xy
                )

                B1dGd_eq_coeff!(aux.AA, aux.BB, aux.CC, aux.SS, vars, sys.gridtype)

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

            solve_lin_system!(aux.A_mat2, aux.b_vec2)

            @inbounds @simd for aa in eachindex(uu)
                bulk.B1d[aa,i,j] = aux.b_vec2[aa]
                bulk.Gd[aa,i,j]  = aux.b_vec2[aa+Nu]
            end

        end
    end

    nothing
end

function solve_phid!(bulk::Bulk, BC::Bulk, gauge::Gauge,
                     nested::Nested, evoleq::EvolEq)
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

    potential = evoleq.potential
    phi0      = evoleq.phi0

    # if phi0 = 0 set phid to zero and return
    if abs(phi0) < 1e-9
        fill!(bulk.phid, 0)
        return
    end

    @fastmath @inbounds @threads for j in eachindex(yy)
        @inbounds for i in eachindex(xx)
            id  = Threads.threadid()
            aux = aux_acc[id]

            xi    = gauge.xi[1,i,j]
            xi_x  = Dx(gauge.xi, 1,i,j)
            xi_y  = Dy(gauge.xi, 1,i,j)
            xi_xx = Dxx(gauge.xi, 1,i,j)
            xi_yy = Dyy(gauge.xi, 1,i,j)
            xi_xy = Dx(Dy, gauge.xi, 1,i,j)

            @inbounds @simd for a in eachindex(uu)
                u     = uu[a]
                u2    = u * u
                u3    = u * u2
                u4    = u2 * u2

                B1    = bulk.B1[a,i,j]
                B2    = bulk.B2[a,i,j]
                G     = bulk.G[a,i,j]
                phi   = bulk.phi[a,i,j]
                S     = bulk.S[a,i,j]
                Fx    = bulk.Fx[a,i,j]
                Fy    = bulk.Fy[a,i,j]
                Sd    = bulk.Sd[a,i,j]

                # r derivatives

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

                # x and y derivatives

                B1_x       = Dx(bulk.B1, a,i,j)
                B2_x       = Dx(bulk.B2, a,i,j)
                G_x        = Dx(bulk.G,  a,i,j)
                phi_x      = Dx(bulk.phi,a,i,j)
                S_x        = Dx(bulk.S,  a,i,j)
                Fx_x       = Dx(bulk.Fx, a,i,j)
                Fy_x       = Dx(bulk.Fy, a,i,j)

                B1_y       = Dy(bulk.B1, a,i,j)
                B2_y       = Dy(bulk.B2, a,i,j)
                G_y        = Dy(bulk.G,  a,i,j)
                phi_y      = Dy(bulk.phi,a,i,j)
                S_y        = Dy(bulk.S,  a,i,j)
                Fx_y       = Dy(bulk.Fx, a,i,j)
                Fy_y       = Dy(bulk.Fy, a,i,j)

                B1p_x      = -u2 * Dx(Du_B1, a,i,j)
                B2p_x      = -u2 * Dx(Du_B2, a,i,j)
                Gp_x       = -u2 * Dx(Du_G,  a,i,j)
                phip_x     = -u2 * Dx(Du_phi,a,i,j)
                Sp_x       = -u2 * Dx(Du_S,  a,i,j)
                Fxp_x      = -u2 * Dx(Du_Fx,  a,i,j)
                Fyp_x      = -u2 * Dx(Du_Fy,  a,i,j)

                B1p_y      = -u2 * Dy(Du_B1, a,i,j)
                B2p_y      = -u2 * Dy(Du_B2, a,i,j)
                Gp_y       = -u2 * Dy(Du_G,  a,i,j)
                phip_y     = -u2 * Dy(Du_phi,a,i,j)
                Sp_y       = -u2 * Dy(Du_S,  a,i,j)
                Fxp_y      = -u2 * Dy(Du_Fx, a,i,j)
                Fyp_y      = -u2 * Dy(Du_Fy, a,i,j)

                B1_xx      = Dxx(bulk.B1, a,i,j)
                B2_xx      = Dxx(bulk.B2, a,i,j)
                G_xx       = Dxx(bulk.G,  a,i,j)
                phi_xx     = Dxx(bulk.phi,a,i,j)
                S_xx       = Dxx(bulk.S,  a,i,j)

                B1_yy      = Dyy(bulk.B1, a,i,j)
                B2_yy      = Dyy(bulk.B2, a,i,j)
                G_yy       = Dyy(bulk.G,  a,i,j)
                phi_yy     = Dyy(bulk.phi,a,i,j)
                S_yy       = Dyy(bulk.S,  a,i,j)

                B2_xy      = Dx(Dy, bulk.B2, a,i,j)
                G_xy       = Dx(Dy, bulk.G,  a,i,j)
                phi_xy     = Dx(Dy, bulk.phi,a,i,j)
                S_xy       = Dx(Dy, bulk.S,  a,i,j)

                vars = (
                    potential, phi0, u, xi, xi_x, xi_y, xi_xx, xi_yy, xi_xy,
                    B1     ,    B2     ,    G      ,    phi    ,    S      ,    Fx     ,    Fy     ,  Sd,
                    B1p    ,    B2p    ,    Gp     ,    phip   ,    Sp     ,    Fxp    ,    Fyp    ,
                    B1pp   ,    B2pp   ,    Gpp    ,    phipp  ,    Spp    ,    Fxpp   ,    Fypp   ,
                    B1_x   ,    B2_x   ,    G_x    ,    phi_x  ,    S_x    ,    Fx_x   ,    Fy_x   ,
                    B1_y   ,    B2_y   ,    G_y    ,    phi_y  ,    S_y    ,    Fx_y   ,    Fy_y   ,
                    B1p_x  ,    B2p_x  ,    Gp_x   ,    phip_x ,    Sp_x   ,    Fxp_x  ,    Fyp_x  ,
                    B1p_y  ,    B2p_y  ,    Gp_y   ,    phip_y ,    Sp_y   ,    Fxp_y  ,    Fyp_y  ,
                    B1_xx  ,    B2_xx  ,    G_xx   ,    phi_xx ,    S_xx   ,
                    B1_yy  ,    B2_yy  ,    G_yy   ,    phi_yy ,    S_yy   ,
                                B2_xy  ,    G_xy   ,    phi_xy,     S_xy
                )

                phid_eq_coeff!(aux.ABCS, vars, sys.gridtype)

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

            solve_lin_system!(aux.A_mat, aux.b_vec)

            @inbounds @simd for aa in eachindex(uu)
                bulk.phid[aa,i,j] = aux.b_vec[aa]
            end

        end
    end

    nothing
end

function solve_A!(bulk::Bulk, BC::Bulk, dBC::Bulk, gauge::Gauge,
                  nested::Nested, evoleq::EvolEq)
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

    potential = evoleq.potential
    phi0      = evoleq.phi0

    @fastmath @inbounds @threads for j in eachindex(yy)
        @inbounds for i in eachindex(xx)
            id  = Threads.threadid()
            aux = aux_acc[id]

            xi    = gauge.xi[1,i,j]
            xi_x  = Dx(gauge.xi, 1,i,j)
            xi_y  = Dy(gauge.xi, 1,i,j)
            xi_xx = Dxx(gauge.xi, 1,i,j)
            xi_yy = Dyy(gauge.xi, 1,i,j)
            xi_xy = Dx(Dy, gauge.xi, 1,i,j)

            @inbounds @simd for a in eachindex(uu)
                u     = uu[a]
                u2    = u * u
                u3    = u * u2
                u4    = u2 * u2

                B1    = bulk.B1[a,i,j]
                B2    = bulk.B2[a,i,j]
                G     = bulk.G[a,i,j]
                phi   = bulk.phi[a,i,j]
                S     = bulk.S[a,i,j]
                Fx    = bulk.Fx[a,i,j]
                Fy    = bulk.Fy[a,i,j]
                Sd    = bulk.Sd[a,i,j]
                B1d   = bulk.B1d[a,i,j]
                B2d   = bulk.B2d[a,i,j]
                Gd    = bulk.Gd[a,i,j]
                phid  = bulk.phid[a,i,j]

                # r derivatives

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

                # x and y derivatives

                B1_x       = Dx(bulk.B1, a,i,j)
                B2_x       = Dx(bulk.B2, a,i,j)
                G_x        = Dx(bulk.G,  a,i,j)
                phi_x      = Dx(bulk.phi,a,i,j)
                S_x        = Dx(bulk.S,  a,i,j)
                Fx_x       = Dx(bulk.Fx, a,i,j)
                Fy_x       = Dx(bulk.Fy, a,i,j)

                B1_y       = Dy(bulk.B1, a,i,j)
                B2_y       = Dy(bulk.B2, a,i,j)
                G_y        = Dy(bulk.G,  a,i,j)
                phi_y      = Dy(bulk.phi,a,i,j)
                S_y        = Dy(bulk.S,  a,i,j)
                Fx_y       = Dy(bulk.Fx, a,i,j)
                Fy_y       = Dy(bulk.Fy, a,i,j)

                B1p_x      = -u2 * Dx(Du_B1, a,i,j)
                B2p_x      = -u2 * Dx(Du_B2, a,i,j)
                Gp_x       = -u2 * Dx(Du_G,  a,i,j)
                phip_x     = -u2 * Dx(Du_phi,a,i,j)
                Sp_x       = -u2 * Dx(Du_S,  a,i,j)
                Fxp_x      = -u2 * Dx(Du_Fx,  a,i,j)
                Fyp_x      = -u2 * Dx(Du_Fy,  a,i,j)

                B1p_y      = -u2 * Dy(Du_B1, a,i,j)
                B2p_y      = -u2 * Dy(Du_B2, a,i,j)
                Gp_y       = -u2 * Dy(Du_G,  a,i,j)
                phip_y     = -u2 * Dy(Du_phi,a,i,j)
                Sp_y       = -u2 * Dy(Du_S,  a,i,j)
                Fxp_y      = -u2 * Dy(Du_Fx, a,i,j)
                Fyp_y      = -u2 * Dy(Du_Fy, a,i,j)

                B1_xx      = Dxx(bulk.B1, a,i,j)
                B2_xx      = Dxx(bulk.B2, a,i,j)
                G_xx       = Dxx(bulk.G,  a,i,j)
                phi_xx     = Dxx(bulk.phi,a,i,j)
                S_xx       = Dxx(bulk.S,  a,i,j)

                B1_yy      = Dyy(bulk.B1, a,i,j)
                B2_yy      = Dyy(bulk.B2, a,i,j)
                G_yy       = Dyy(bulk.G,  a,i,j)
                phi_yy     = Dyy(bulk.phi,a,i,j)
                S_yy       = Dyy(bulk.S,  a,i,j)

                B2_xy      = Dx(Dy, bulk.B2, a,i,j)
                G_xy       = Dx(Dy, bulk.G,  a,i,j)
                phi_xy     = Dx(Dy, bulk.phi,a,i,j)
                S_xy       = Dx(Dy, bulk.S,  a,i,j)

                vars = (
                    potential, phi0, u, xi, xi_x, xi_y, xi_xx, xi_yy, xi_xy,
                    B1   , B2   , G   , phi   , S    , Fx    , Fy    , Sd, B1d, B2d, Gd, phid,
                    B1p  , B2p  , Gp  , phip  , Sp   , Fxp   , Fyp   ,
                    B1pp , B2pp , Gpp , phipp , Spp  , Fxpp  , Fypp  ,
                    B1_x , B2_x , G_x , phi_x , S_x  , Fx_x  , Fy_x  ,
                    B1_y , B2_y , G_y , phi_y , S_y  , Fx_y  , Fy_y  ,
                    B1p_x, B2p_x, Gp_x, phip_x, Sp_x , Fxp_x , Fyp_x ,
                    B1p_y, B2p_y, Gp_y, phip_y, Sp_y , Fxp_y , Fyp_y ,
                    B1_xx, B2_xx, G_xx, phi_xx, S_xx ,
                    B1_yy, B2_yy, G_yy, phi_yy, S_yy ,
                           B2_xy, G_xy, phi_xy, S_xy
                )

                A_eq_coeff!(aux.ABCS, vars, sys.gridtype)

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

            solve_lin_system!(aux.A_mat, aux.b_vec)

            @inbounds @simd for aa in eachindex(uu)
                bulk.A[aa,i,j] = aux.b_vec[aa]
            end

        end
    end

    nothing
end

function solve_nested!(bulk::Bulk, BC::Bulk, dBC::Bulk, gauge::Gauge,
                       nested::Nested, evoleq::EvolEq)
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
    solve_S!(bulk, BC, dBC, gauge, nested, evoleq)

    # take u-derivatives of S
    @sync begin
        @spawn mul!(Du_S,   Du,  bulk.S)
        @spawn mul!(Duu_S,  Duu, bulk.S)
    end

    # solve for Fx and Fy
    solve_Fxy!(bulk, BC, dBC, gauge, nested, evoleq)

    # take u-derivatives of Fx and Fy
    @sync begin
        @spawn mul!(Du_Fx,   Du,  bulk.Fx)
        @spawn mul!(Du_Fy,   Du,  bulk.Fy)
        @spawn mul!(Duu_Fx,  Duu, bulk.Fx)
        @spawn mul!(Duu_Fy,  Duu, bulk.Fy)
    end

    # solve for Sd
    solve_Sd!(bulk, BC, gauge, nested, evoleq)

    # solving for B2d, (B1d,Gd) and phid are independent processes. we can
    # therefore @spawn, here
    @sync begin
        @spawn solve_B2d!(bulk, BC, gauge, nested, evoleq)
        @spawn solve_B1dGd!(bulk, BC, gauge, nested, evoleq)
        @spawn solve_phid!(bulk, BC, gauge, nested, evoleq)
    end

    # solve for A
    solve_A!(bulk, BC, dBC, gauge, nested, evoleq)

    nothing
end


function syncBCs!(BC::Bulk, dBC::Bulk, bulk::Bulk, nested::Nested)
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

function set_innerBCs!(BC::Bulk, dBC::Bulk, bulk::Bulk,
                       boundary::Boundary, gauge::Gauge,
                       nested::Nested, evoleq::EvolEq)
    _, Nx, Ny = size(nested.sys)

    Dx  = nested.sys.Dx
    Dy  = nested.sys.Dy

    phi0  = evoleq.phi0
    phi02 = phi0 * phi0
    phi03 = phi0 * phi02
    phi04 = phi02 * phi02

    @fastmath @inbounds @threads for j in 1:Ny
        @inbounds @simd for i in 1:Nx
            xi      = gauge.xi[1,i,j]
            xi3     = xi*xi*xi

            phi     = bulk.phi[1,i,j]
            phi_u   = nested.Du_phi[1,i,j]

            phi_x   = Dx(bulk.phi, 1,i,j)
            phi_y   = Dy(bulk.phi, 1,i,j)

            xi_x    = Dx(gauge.xi, 1,i,j)
            xi_y    = Dy(gauge.xi, 1,i,j)

            b14     = bulk.B1[1,i,j]
            b24     = bulk.B2[1,i,j]
            g4      = bulk.G[1,i,j]

            a4      = boundary.a4[1,i,j]

            fx2     = boundary.fx2[1,i,j]
            fy2     = boundary.fy2[1,i,j]

            fx2_x   = Dx(boundary.fx2, 1,i,j)
            fy2_y   = Dy(boundary.fy2, 1,i,j)

            b14_x   = Dx(bulk.B1,1,i,j)
            b24_x   = Dx(bulk.B2,1,i,j)
            phi_x   = Dx(bulk.phi,1,i,j)
            g4_x    = Dx(bulk.G,1,i,j)

            b14_y   = Dy(bulk.B1,1,i,j)
            b24_y   = Dy(bulk.B2,1,i,j)
            phi_y   = Dy(bulk.phi,1,i,j)
            g4_y    = Dy(bulk.G,1,i,j)

            phi2    = phi03 * phi - phi0 * xi * xi
            phi2_x  = phi03 * phi_x - 2 * phi0 * xi * xi_x
            phi2_y  = phi03 * phi_y - 2 * phi0 * xi * xi_y

            BC.S[i,j]  = phi04 * (1 - 18 * phi) / 54
            dBC.S[i,j] = phi02 * (-12 * xi3 +
                                  phi02 * xi * (18 * phi - 5) -
                                  24 * phi02 * phi_u) / 90

            BC.Fx[i,j]  = fx2
            dBC.Fx[i,j] = -2 * fx2 * xi - 12 / 15 * (b14_x + b24_x - g4_y) +
                4/15 * phi0 * phi2_x

            BC.Fy[i,j]  = fy2
            dBC.Fy[i,j] = -2 * fy2 * xi - 12 / 15 * (-b14_y + b24_y - g4_x) +
                4/15 * phi0 * phi2_y

            BC.Sd[i,j] = a4 / 2

            BC.B2d[i,j] = -2 * b24
            BC.B1d[i,j] = -2 * b14
            BC.Gd[i,j]  = -2 * g4

            BC.A[i,j]  = a4
            dBC.A[i,j] = -2 * xi * a4 - 2 * phi0 * xi * phi2 -
                2/3 * (phi02 * xi3 + fx2_x + fy2_y + phi04 * phi_u)

        end
    end

    # separate the phid case, since the if statement may prevent loop
    # vectorization

    if abs(phi0) < 1e-9
        fill!(BC.phid, 0)
        return
    end

    @fastmath @inbounds for j in 1:Ny
        @inbounds @simd for i in 1:Nx
            xi      = gauge.xi[1,i,j]
            phi     = bulk.phi[1,i,j]
            phi2    = phi03 * phi - phi0 * xi * xi

            BC.phid[i,j] = 1/3 - 3/2 * phi2 / phi03
        end
    end

    nothing
end

function set_outerBCs!(BC::Bulk, dBC::Bulk, bulk::Bulk,
                       gauge::Gauge, nested::Nested, evoleq::EvolEq)
    Nu, Nx, Ny = size(nested.sys)
    Du = nested.sys.Du

    phi0 = evoleq.phi0

    # we are here assuming that the inner and outer grids merely touch at the
    # interface, so we pass the values at this point without any interpolation
    u0 = nested.sys.ucoord[end]

    @fastmath @inbounds @threads for j in 1:Ny
        @inbounds @simd for i in 1:Nx
            xi     = gauge.xi[1,i,j]

            S      = bulk.S[end,i,j]
            Fx     = bulk.Fx[end,i,j]
            Fy     = bulk.Fy[end,i,j]
            Sd     = bulk.Sd[end,i,j]
            B1d    = bulk.B1d[end,i,j]
            B2d    = bulk.B2d[end,i,j]
            Gd     = bulk.Gd[end,i,j]
            phid   = bulk.phid[end,i,j]
            A      = bulk.A[end,i,j]

            S_u    = nested.Du_S[end,i,j]
            Fx_u   = nested.Du_Fx[end,i,j]
            Fy_u   = nested.Du_Fy[end,i,j]
            A_u    = Du(bulk.A, Nu,i,j)

            BC.S[i,j]  = S_inner_to_outer(S, u0, xi, phi0)
            dBC.S[i,j] = S_u_inner_to_outer(S_u, S, u0, xi, phi0)

            BC.Fx[i,j]  = F_inner_to_outer(Fx, u0)
            BC.Fy[i,j]  = F_inner_to_outer(Fy, u0)
            dBC.Fx[i,j] = F_u_inner_to_outer(Fx_u, Fx, u0)
            dBC.Fy[i,j] = F_u_inner_to_outer(Fy_u, Fy, u0)

            BC.Sd[i,j]  = Sd_inner_to_outer(Sd, u0, xi, phi0)

            # B1d, B2d, and Gd transform in the same way
            BC.B2d[i,j] = Bd_inner_to_outer(B2d, u0)
            BC.B1d[i,j] = Bd_inner_to_outer(B1d, u0)
            BC.Gd[i,j]  = Bd_inner_to_outer(Gd, u0)

            BC.phid[i,j] = phid_inner_to_outer(phid, u0, phi0)

            BC.A[i,j]  = A_inner_to_outer(A, u0, xi, phi0)
            dBC.A[i,j] = A_u_inner_to_outer(A_u, A, u0, xi, phi0)
        end
    end

    nothing
end


# We assume that the first entry on these arrays is the inner grid, and that
# there is only one domain spanning this grid. If we ever change this
# construction we must remember to make the appropriate changes here.
function solve_nested!(bulks, BCs, dBCs, boundary::Boundary,
                       gauge::Gauge, nesteds, evoleq::EvolEq)
    Nsys = length(nesteds)

    set_innerBCs!(BCs[1], dBCs[1], bulks[1], boundary, gauge, nesteds[1], evoleq)

    solve_nested!(bulks[1], BCs[1], dBCs[1], gauge, nesteds[1], evoleq)

    set_outerBCs!(BCs[2], dBCs[2], bulks[1], gauge, nesteds[1], evoleq)

    for i in 2:Nsys-1
        solve_nested!(bulks[i], BCs[i], dBCs[i], gauge, nesteds[i], evoleq)
        syncBCs!(BCs[i+1], dBCs[i+1], bulks[i], nesteds[i])
    end
    solve_nested!(bulks[Nsys], BCs[Nsys], dBCs[Nsys], gauge, nesteds[Nsys], evoleq)

    nothing
end

function nested_solver(systems::SystemPartition, evoleq::EvolEq)
    sys1  = systems[1]
    T     = Jecco.coord_eltype(sys1.ucoord)
    _, Nx, Ny = size(sys1)

    nesteds = Tuple([Nested(sys) for sys in systems])
    BCs     = Tuple([Bulk{T}(undef, Nx, Ny) for sys in systems])
    dBCs    = Tuple([Bulk{T}(undef, Nx, Ny) for sys in systems])

    function (bulks, boundary::Boundary, gauge::Gauge)
        solve_nested!(bulks, BCs, dBCs, boundary, gauge, nesteds, evoleq)
        nothing
    end
end

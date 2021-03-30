
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
function Aux{T}(sys::System) where {T}
    Nu, Nx, Ny = size(sys)
    nt = Threads.nthreads()
    # pre-allocate thread-local aux quantities
    [Aux{T}(Nu) for _ in 1:nt]
end

struct BC{T}
    S    :: Array{T,2}
    Fx   :: Array{T,2}
    Fy   :: Array{T,2}
    B1d  :: Array{T,2}
    B2d  :: Array{T,2}
    Gd   :: Array{T,2}
    phid :: Array{T,2}
    Sd   :: Array{T,2}
    A    :: Array{T,2}
    S_u  :: Array{T,2}
    Fx_u :: Array{T,2}
    Fy_u :: Array{T,2}
    A_u  :: Array{T,2}
end
function BC{T}(Nx::Int, Ny::Int) where {T<:Real}
    S    = Array{T}(undef, Nx, Ny)
    Fx   = Array{T}(undef, Nx, Ny)
    Fy   = Array{T}(undef, Nx, Ny)
    B1d  = Array{T}(undef, Nx, Ny)
    B2d  = Array{T}(undef, Nx, Ny)
    Gd   = Array{T}(undef, Nx, Ny)
    phid = Array{T}(undef, Nx, Ny)
    Sd   = Array{T}(undef, Nx, Ny)
    A    = Array{T}(undef, Nx, Ny)
    S_u  = Array{T}(undef, Nx, Ny)
    Fx_u = Array{T}(undef, Nx, Ny)
    Fy_u = Array{T}(undef, Nx, Ny)
    A_u  = Array{T}(undef, Nx, Ny)
    BC{T}(S, Fx, Fy, B1d, B2d, Gd, phid, Sd, A, S_u, Fx_u, Fy_u, A_u)
end

function BCs(systems::SystemPartition)
    sys1  = systems[1]
    T     = Jecco.coord_eltype(sys1.ucoord)
    _, Nx, Ny = size(sys1)

    [BC{T}(Nx, Ny) for sys in systems]
end


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

function solve_S!(bulk::Bulk, bc::BC, gauge::Gauge, deriv::BulkDeriv, aux_acc,
                  sys::System, evoleq::AffineNull)
    Nu, Nx, Ny = size(sys)

    B1GF    = getB1(bulk)
    B2GF    = getB2(bulk)
    GGF     = getG(bulk)
    phiGF   = getphi(bulk)
    SGF     = getS(bulk)
    xiGF    = getxi(gauge)

    Du_B1   = deriv.Du_B1
    Du_B2   = deriv.Du_B2
    Du_G    = deriv.Du_G
    Du_phi  = deriv.Du_phi
    # Duu_B1  = deriv.Duu_B1
    # Duu_B2  = deriv.Duu_B2
    # Duu_G   = deriv.Duu_G

    Du  = sys.Du
    Duu = sys.Duu
    # Dx  = sys.Dx
    # Dxx = sys.Dxx
    # Dy  = sys.Dy
    # Dyy = sys.Dyy

    phi0  = evoleq.phi0

    @fastmath @inbounds @threads for j in 1:Ny
        @inbounds for i in 1:Nx
            id  = Threads.threadid()
            aux = aux_acc[id]

            xi  = xiGF[1,i,j]

            @inbounds @simd for a in 1:Nu
                u     = sys.ucoord[a]

                B1    = B1GF[a,i,j]
                B2    = B2GF[a,i,j]
                G     = GGF[a,i,j]
                phi   = phiGF[a,i,j]

                B1p   = -u*u * Du_B1[a,i,j]
                B2p   = -u*u * Du_B2[a,i,j]
                Gp    = -u*u * Du_G[a,i,j]
                phip  = -u*u * Du_phi[a,i,j]

                vars = (phi0, u, xi, B1, B1p, B2, B2p, G, Gp, phi, phip)

                S_eq_coeff!(aux.ABCS, vars, sys.gridtype)

                aux.b_vec[a]   = -aux.ABCS[4]
                @inbounds @simd for aa in 1:Nu
                    aux.A_mat[a,aa] = aux.ABCS[1] * Duu[a,aa] + aux.ABCS[2] * Du[a,aa]
                end
                aux.A_mat[a,a] += aux.ABCS[3]
            end

            # BC

            aux.b_vec[1]    = bc.S[i,j]
            aux.A_mat[1,:] .= 0.0
            aux.A_mat[1,1]  = 1.0

            aux.b_vec[end]    = bc.S_u[i,j]
            @inbounds @simd for aa in 1:Nu
                aux.A_mat[end,aa]  = Du[1,aa]
            end

            solve_lin_system!(aux.A_mat, aux.b_vec)

            @inbounds @simd for aa in 1:Nu
                SGF[aa,i,j] = aux.b_vec[aa]
            end

        end
    end

    nothing
end

function solve_Fxy!(bulk::Bulk, bc::BC, gauge::Gauge, deriv::BulkDeriv, aux_acc,
                    sys::System, evoleq::AffineNull)
    Nu, Nx, Ny = size(sys)

    B1GF    = getB1(bulk)
    B2GF    = getB2(bulk)
    GGF     = getG(bulk)
    phiGF   = getphi(bulk)
    SGF     = getS(bulk)
    FxGF    = getFx(bulk)
    FyGF    = getFy(bulk)
    xiGF    = getxi(gauge)

    Du_B1   = deriv.Du_B1
    Du_B2   = deriv.Du_B2
    Du_G    = deriv.Du_G
    Du_phi  = deriv.Du_phi
    Du_S    = deriv.Du_S
    Duu_B1  = deriv.Duu_B1
    Duu_B2  = deriv.Duu_B2
    Duu_G   = deriv.Duu_G
    Duu_phi = deriv.Duu_phi
    Duu_S   = deriv.Duu_S

    Du  = sys.Du
    Duu = sys.Duu
    Dx  = sys.Dx
    Dxx = sys.Dxx
    Dy  = sys.Dy
    Dyy = sys.Dyy

    phi0  = evoleq.phi0

    @fastmath @inbounds @threads for j in 1:Ny
        @inbounds for i in 1:Nx
            id  = Threads.threadid()
            aux = aux_acc[id]

            xi    = xiGF[1,i,j]
            xi_x  = Dx(xiGF, 1,i,j)
            xi_y  = Dy(xiGF, 1,i,j)

            @inbounds @simd for a in 1:Nu
                u     = sys.ucoord[a]
                u2    = u * u
                u3    = u * u2
                u4    = u2 * u2

                B1    = B1GF[a,i,j]
                B1p   = -u2 * Du_B1[a,i,j]
                B1_x  = Dx(B1GF, a,i,j)
                B1_y  = Dy(B1GF, a,i,j)
                B1pp  = 2*u3 * Du_B1[a,i,j] + u4 * Duu_B1[a,i,j]
                B1p_x = -u2 * Dx(Du_B1, a,i,j)
                B1p_y = -u2 * Dy(Du_B1, a,i,j)

                B2    = B2GF[a,i,j]
                B2p   = -u2 * Du_B2[a,i,j]
                B2_x  = Dx(B2GF, a,i,j)
                B2_y  = Dy(B2GF, a,i,j)
                B2pp  = 2*u3 * Du_B2[a,i,j] + u4 * Duu_B2[a,i,j]
                B2p_x = -u2 * Dx(Du_B2, a,i,j)
                B2p_y = -u2 * Dy(Du_B2, a,i,j)

                G     = GGF[a,i,j]
                Gp    = -u2 * Du_G[a,i,j]
                G_x   = Dx(GGF, a,i,j)
                G_y   = Dy(GGF, a,i,j)
                Gpp   = 2*u3 * Du_G[a,i,j] + u4 * Duu_G[a,i,j]
                Gp_x  = -u2 * Dx(Du_G, a,i,j)
                Gp_y  = -u2 * Dy(Du_G, a,i,j)

                phi   = phiGF[a,i,j]
                phip  = -u2 * Du_phi[a,i,j]
                phi_x = Dx(phiGF, a,i,j)
                phi_y = Dy(phiGF, a,i,j)
                # phipp   = 2*u3 * Du_phi[a,i,j] + u4 * Duu_phi[a,i,j]
                # phip_x  = -u2 * Dx(Du_phi, a,i,j)
                # phip_y  = -u2 * Dy(Du_phi, a,i,j)

                S     = SGF[a,i,j]
                Sp    = -u2 * Du_S[a,i,j]
                S_x   = Dx(SGF, a,i,j)
                S_y   = Dy(SGF, a,i,j)
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
                @inbounds @simd for aa in 1:Nu
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

            aux.b_vec2[1]     = bc.Fx[i,j]
            aux.A_mat2[1,:]  .= 0.0
            aux.A_mat2[1,1]   = 1.0

            aux.b_vec2[1+Nu]      = bc.Fy[i,j]
            aux.A_mat2[1+Nu,:]   .= 0.0
            aux.A_mat2[1+Nu,1+Nu] = 1.0

            aux.b_vec2[Nu]   = bc.Fx_u[i,j]
            aux.b_vec2[2*Nu] = bc.Fy_u[i,j]

            aux.A_mat2[Nu,:]   .= 0.0
            aux.A_mat2[2*Nu,:] .= 0.0
            @inbounds @simd for aa in 1:Nu
                aux.A_mat2[Nu,aa]      = Du[1,aa]
                aux.A_mat2[2*Nu,aa+Nu] = Du[1,aa]
            end

            solve_lin_system!(aux.A_mat2, aux.b_vec2)

            @inbounds @simd for aa in 1:Nu
                FxGF[aa,i,j] = aux.b_vec2[aa]
                FyGF[aa,i,j] = aux.b_vec2[aa+Nu]
            end

        end
    end

    nothing
end

function solve_Sd!(bulk::Bulk, bc::BC, gauge::Gauge, deriv::BulkDeriv, aux_acc,
                   sys::System, evoleq::AffineNull)
    Nu, Nx, Ny = size(sys)

    B1GF    = getB1(bulk)
    B2GF    = getB2(bulk)
    GGF     = getG(bulk)
    phiGF   = getphi(bulk)
    SGF     = getS(bulk)
    FxGF    = getFx(bulk)
    FyGF    = getFy(bulk)
    SdGF    = getSd(bulk)
    xiGF    = getxi(gauge)

    Du_B1   = deriv.Du_B1
    Du_B2   = deriv.Du_B2
    Du_G    = deriv.Du_G
    Du_phi  = deriv.Du_phi
    Du_S    = deriv.Du_S
    Du_Fx   = deriv.Du_Fx
    Du_Fy   = deriv.Du_Fy
    Duu_B1  = deriv.Duu_B1
    Duu_B2  = deriv.Duu_B2
    Duu_G   = deriv.Duu_G
    Duu_phi = deriv.Duu_phi
    Duu_S   = deriv.Duu_S
    Duu_Fx  = deriv.Duu_Fx
    Duu_Fy  = deriv.Duu_Fy

    Du  = sys.Du
    Duu = sys.Duu
    Dx  = sys.Dx
    Dxx = sys.Dxx
    Dy  = sys.Dy
    Dyy = sys.Dyy

    potential = evoleq.potential
    phi0      = evoleq.phi0

    @fastmath @inbounds @threads for j in 1:Ny
        @inbounds for i in 1:Nx
            id   = Threads.threadid()
            aux  = aux_acc[id]

            xi    = xiGF[1,i,j]
            xi_x  = Dx(xiGF, 1,i,j)
            xi_y  = Dy(xiGF, 1,i,j)
            xi_xx = Dxx(xiGF, 1,i,j)
            xi_yy = Dyy(xiGF, 1,i,j)
            xi_xy = Dx(Dy, xiGF, 1,i,j)

            @inbounds @simd for a in 1:Nu
                u     = sys.ucoord[a]
                u2    = u * u
                u3    = u * u2
                u4    = u2 * u2

                B1    = B1GF[a,i,j]
                B2    = B2GF[a,i,j]
                G     = GGF[a,i,j]
                phi   = phiGF[a,i,j]
                S     = SGF[a,i,j]
                Fx    = FxGF[a,i,j]
                Fy    = FyGF[a,i,j]

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

                B1_x       = Dx(B1GF, a,i,j)
                B2_x       = Dx(B2GF, a,i,j)
                G_x        = Dx(GGF,  a,i,j)
                phi_x      = Dx(phiGF,a,i,j)
                S_x        = Dx(SGF,  a,i,j)
                Fx_x       = Dx(FxGF, a,i,j)
                Fy_x       = Dx(FyGF, a,i,j)

                B1_y       = Dy(B1GF, a,i,j)
                B2_y       = Dy(B2GF, a,i,j)
                G_y        = Dy(GGF,  a,i,j)
                phi_y      = Dy(phiGF,a,i,j)
                S_y        = Dy(SGF,  a,i,j)
                Fx_y       = Dy(FxGF, a,i,j)
                Fy_y       = Dy(FyGF, a,i,j)

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

                B1_xx      = Dxx(B1GF, a,i,j)
                B2_xx      = Dxx(B2GF, a,i,j)
                G_xx       = Dxx(GGF,  a,i,j)
                phi_xx     = Dxx(phiGF,a,i,j)
                S_xx       = Dxx(SGF,  a,i,j)

                B1_yy      = Dyy(B1GF, a,i,j)
                B2_yy      = Dyy(B2GF, a,i,j)
                G_yy       = Dyy(GGF,  a,i,j)
                phi_yy     = Dyy(phiGF,a,i,j)
                S_yy       = Dyy(SGF,  a,i,j)

                B2_xy      = Dx(Dy, B2GF, a,i,j)
                G_xy       = Dx(Dy, GGF,  a,i,j)
                S_xy       = Dx(Dy, SGF,  a,i,j)

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
                @inbounds @simd for aa in 1:Nu
                    aux.A_mat[a,aa] = aux.ABCS[1] * Duu[a,aa] + aux.ABCS[2] * Du[a,aa]
                end
                aux.A_mat[a,a] += aux.ABCS[3]
            end

            # BC (first order equation)

            aux.b_vec[1]    = bc.Sd[i,j]
            aux.A_mat[1,:] .= 0.0
            aux.A_mat[1,1]  = 1.0

            solve_lin_system!(aux.A_mat, aux.b_vec)

            @inbounds @simd for aa in 1:Nu
                SdGF[aa,i,j] = aux.b_vec[aa]
            end

        end
    end

    nothing
end

function solve_B2d!(bulk::Bulk, bc::BC, gauge::Gauge, deriv::BulkDeriv, aux_acc,
                    sys::System, evoleq::AffineNull)
    Nu, Nx, Ny = size(sys)

    B1GF    = getB1(bulk)
    B2GF    = getB2(bulk)
    GGF     = getG(bulk)
    phiGF   = getphi(bulk)
    SGF     = getS(bulk)
    FxGF    = getFx(bulk)
    FyGF    = getFy(bulk)
    SdGF    = getSd(bulk)
    B2dGF   = getB2d(bulk)
    phidGF  = getphid(bulk)
    xiGF    = getxi(gauge)

    Du_B1   = deriv.Du_B1
    Du_B2   = deriv.Du_B2
    Du_G    = deriv.Du_G
    Du_phi  = deriv.Du_phi
    Du_S    = deriv.Du_S
    Du_Fx   = deriv.Du_Fx
    Du_Fy   = deriv.Du_Fy
    Duu_B1  = deriv.Duu_B1
    Duu_B2  = deriv.Duu_B2
    Duu_G   = deriv.Duu_G
    Duu_phi = deriv.Duu_phi
    Duu_S   = deriv.Duu_S
    Duu_Fx  = deriv.Duu_Fx
    Duu_Fy  = deriv.Duu_Fy

    Du  = sys.Du
    Duu = sys.Duu
    Dx  = sys.Dx
    Dxx = sys.Dxx
    Dy  = sys.Dy
    Dyy = sys.Dyy

    phi0  = evoleq.phi0

    @fastmath @inbounds @threads for j in 1:Ny
        @inbounds for i in 1:Nx
            id  = Threads.threadid()
            aux = aux_acc[id]

            xi    = xiGF[1,i,j]
            xi_x  = Dx(xiGF, 1,i,j)
            xi_y  = Dy(xiGF, 1,i,j)
            xi_xx = Dxx(xiGF, 1,i,j)
            xi_yy = Dyy(xiGF, 1,i,j)
            xi_xy = Dx(Dy, xiGF, 1,i,j)

            @inbounds @simd for a in 1:Nu
                u     = sys.ucoord[a]
                u2    = u * u
                u3    = u * u2
                u4    = u2 * u2

                B1    = B1GF[a,i,j]
                B2    = B2GF[a,i,j]
                G     = GGF[a,i,j]
                phi   = phiGF[a,i,j]
                S     = SGF[a,i,j]
                Fx    = FxGF[a,i,j]
                Fy    = FyGF[a,i,j]
                Sd    = SdGF[a,i,j]

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

                B1_x       = Dx(B1GF, a,i,j)
                B2_x       = Dx(B2GF, a,i,j)
                G_x        = Dx(GGF,  a,i,j)
                phi_x      = Dx(phiGF,a,i,j)
                S_x        = Dx(SGF,  a,i,j)
                Fx_x       = Dx(FxGF, a,i,j)
                Fy_x       = Dx(FyGF, a,i,j)

                B1_y       = Dy(B1GF, a,i,j)
                B2_y       = Dy(B2GF, a,i,j)
                G_y        = Dy(GGF,  a,i,j)
                phi_y      = Dy(phiGF,a,i,j)
                S_y        = Dy(SGF,  a,i,j)
                Fx_y       = Dy(FxGF, a,i,j)
                Fy_y       = Dy(FyGF, a,i,j)

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

                B1_xx      = Dxx(B1GF, a,i,j)
                B2_xx      = Dxx(B2GF, a,i,j)
                G_xx       = Dxx(GGF,  a,i,j)
                phi_xx     = Dxx(phiGF,a,i,j)
                S_xx       = Dxx(SGF,  a,i,j)

                B1_yy      = Dyy(B1GF, a,i,j)
                B2_yy      = Dyy(B2GF, a,i,j)
                G_yy       = Dyy(GGF,  a,i,j)
                phi_yy     = Dyy(phiGF,a,i,j)
                S_yy       = Dyy(SGF,  a,i,j)

                B2_xy      = Dx(Dy, B2GF, a,i,j)
                G_xy       = Dx(Dy, GGF,  a,i,j)
                S_xy       = Dx(Dy, SGF,  a,i,j)

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
                @inbounds @simd for aa in 1:Nu
                    aux.A_mat[a,aa] = aux.ABCS[1] * Duu[a,aa] + aux.ABCS[2] * Du[a,aa]
                end
                aux.A_mat[a,a] += aux.ABCS[3]
            end

            # BC (first order equation)

            aux.b_vec[1]    = bc.B2d[i,j]
            aux.A_mat[1,:] .= 0.0
            aux.A_mat[1,1]  = 1.0

            solve_lin_system!(aux.A_mat, aux.b_vec)

            @inbounds @simd for aa in 1:Nu
                B2dGF[aa,i,j] = aux.b_vec[aa]
            end

        end
    end

    nothing
end

function solve_B1dGd!(bulk::Bulk, bc::BC, gauge::Gauge, deriv::BulkDeriv, aux_acc,
                      sys::System, evoleq::AffineNull)
    Nu, Nx, Ny = size(sys)

    B1GF    = getB1(bulk)
    B2GF    = getB2(bulk)
    GGF     = getG(bulk)
    phiGF   = getphi(bulk)
    SGF     = getS(bulk)
    FxGF    = getFx(bulk)
    FyGF    = getFy(bulk)
    SdGF    = getSd(bulk)
    B1dGF   = getB1d(bulk)
    GdGF    = getGd(bulk)
    xiGF    = getxi(gauge)

    Du_B1   = deriv.Du_B1
    Du_B2   = deriv.Du_B2
    Du_G    = deriv.Du_G
    Du_phi  = deriv.Du_phi
    Du_S    = deriv.Du_S
    Du_Fx   = deriv.Du_Fx
    Du_Fy   = deriv.Du_Fy
    Duu_B1  = deriv.Duu_B1
    Duu_B2  = deriv.Duu_B2
    Duu_G   = deriv.Duu_G
    Duu_phi = deriv.Duu_phi
    Duu_S   = deriv.Duu_S
    Duu_Fx  = deriv.Duu_Fx
    Duu_Fy  = deriv.Duu_Fy

    Du  = sys.Du
    Duu = sys.Duu
    Dx  = sys.Dx
    Dxx = sys.Dxx
    Dy  = sys.Dy
    Dyy = sys.Dyy

    phi0  = evoleq.phi0

    @fastmath @inbounds @threads for j in 1:Ny
        @inbounds for i in 1:Nx
            id  = Threads.threadid()
            aux = aux_acc[id]

            xi    = xiGF[1,i,j]
            xi_x  = Dx(xiGF, 1,i,j)
            xi_y  = Dy(xiGF, 1,i,j)
            xi_xx = Dxx(xiGF, 1,i,j)
            xi_yy = Dyy(xiGF, 1,i,j)
            xi_xy = Dx(Dy, xiGF, 1,i,j)

            @inbounds @simd for a in 1:Nu
                u     = sys.ucoord[a]
                u2    = u * u
                u3    = u * u2
                u4    = u2 * u2

                B1    = B1GF[a,i,j]
                B2    = B2GF[a,i,j]
                G     = GGF[a,i,j]
                phi   = phiGF[a,i,j]
                S     = SGF[a,i,j]
                Fx    = FxGF[a,i,j]
                Fy    = FyGF[a,i,j]
                Sd    = SdGF[a,i,j]

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

                B1_x       = Dx(B1GF, a,i,j)
                B2_x       = Dx(B2GF, a,i,j)
                G_x        = Dx(GGF,  a,i,j)
                phi_x      = Dx(phiGF,a,i,j)
                S_x        = Dx(SGF,  a,i,j)
                Fx_x       = Dx(FxGF, a,i,j)
                Fy_x       = Dx(FyGF, a,i,j)

                B1_y       = Dy(B1GF, a,i,j)
                B2_y       = Dy(B2GF, a,i,j)
                G_y        = Dy(GGF,  a,i,j)
                phi_y      = Dy(phiGF,a,i,j)
                S_y        = Dy(SGF,  a,i,j)
                Fx_y       = Dy(FxGF, a,i,j)
                Fy_y       = Dy(FyGF, a,i,j)

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

                B1_xx      = Dxx(B1GF, a,i,j)
                B2_xx      = Dxx(B2GF, a,i,j)
                G_xx       = Dxx(GGF,  a,i,j)
                phi_xx     = Dxx(phiGF,a,i,j)
                S_xx       = Dxx(SGF,  a,i,j)

                B1_yy      = Dyy(B1GF, a,i,j)
                B2_yy      = Dyy(B2GF, a,i,j)
                G_yy       = Dyy(GGF,  a,i,j)
                phi_yy     = Dyy(phiGF,a,i,j)
                S_yy       = Dyy(SGF,  a,i,j)

                B2_xy      = Dx(Dy, B2GF, a,i,j)
                G_xy       = Dx(Dy, GGF,  a,i,j)
                S_xy       = Dx(Dy, SGF,  a,i,j)

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
                @inbounds @simd for aa in 1:Nu
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

            aux.b_vec2[1]     = bc.B1d[i,j]
            aux.A_mat2[1,:]  .= 0.0
            aux.A_mat2[1,1]   = 1.0

            aux.b_vec2[1+Nu]      = bc.Gd[i,j]
            aux.A_mat2[1+Nu,:]   .= 0.0
            aux.A_mat2[1+Nu,1+Nu] = 1.0

            solve_lin_system!(aux.A_mat2, aux.b_vec2)

            @inbounds @simd for aa in 1:Nu
                B1dGF[aa,i,j] = aux.b_vec2[aa]
                GdGF[aa,i,j]  = aux.b_vec2[aa+Nu]
            end

        end
    end

    nothing
end

function solve_phid!(bulk::Bulk, bc::BC, gauge::Gauge, deriv::BulkDeriv, aux_acc,
                     sys::System, evoleq::AffineNull)
    Nu, Nx, Ny = size(sys)

    B1GF    = getB1(bulk)
    B2GF    = getB2(bulk)
    GGF     = getG(bulk)
    phiGF   = getphi(bulk)
    SGF     = getS(bulk)
    FxGF    = getFx(bulk)
    FyGF    = getFy(bulk)
    SdGF    = getSd(bulk)
    phidGF  = getphid(bulk)
    xiGF    = getxi(gauge)

    Du_B1   = deriv.Du_B1
    Du_B2   = deriv.Du_B2
    Du_G    = deriv.Du_G
    Du_phi  = deriv.Du_phi
    Du_S    = deriv.Du_S
    Du_Fx   = deriv.Du_Fx
    Du_Fy   = deriv.Du_Fy
    Duu_B1  = deriv.Duu_B1
    Duu_B2  = deriv.Duu_B2
    Duu_G   = deriv.Duu_G
    Duu_phi = deriv.Duu_phi
    Duu_S   = deriv.Duu_S
    Duu_Fx  = deriv.Duu_Fx
    Duu_Fy  = deriv.Duu_Fy

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
        fill!(phidGF, 0)
        return
    end

    @fastmath @inbounds @threads for j in 1:Ny
        @inbounds for i in 1:Nx
            id  = Threads.threadid()
            aux = aux_acc[id]

            xi    = xiGF[1,i,j]
            xi_x  = Dx(xiGF, 1,i,j)
            xi_y  = Dy(xiGF, 1,i,j)
            xi_xx = Dxx(xiGF, 1,i,j)
            xi_yy = Dyy(xiGF, 1,i,j)
            xi_xy = Dx(Dy, xiGF, 1,i,j)

            @inbounds @simd for a in 1:Nu
                u     = sys.ucoord[a]
                u2    = u * u
                u3    = u * u2
                u4    = u2 * u2

                B1    = B1GF[a,i,j]
                B2    = B2GF[a,i,j]
                G     = GGF[a,i,j]
                phi   = phiGF[a,i,j]
                S     = SGF[a,i,j]
                Fx    = FxGF[a,i,j]
                Fy    = FyGF[a,i,j]
                Sd    = SdGF[a,i,j]

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

                B1_x       = Dx(B1GF, a,i,j)
                B2_x       = Dx(B2GF, a,i,j)
                G_x        = Dx(GGF,  a,i,j)
                phi_x      = Dx(phiGF,a,i,j)
                S_x        = Dx(SGF,  a,i,j)
                Fx_x       = Dx(FxGF, a,i,j)
                Fy_x       = Dx(FyGF, a,i,j)

                B1_y       = Dy(B1GF, a,i,j)
                B2_y       = Dy(B2GF, a,i,j)
                G_y        = Dy(GGF,  a,i,j)
                phi_y      = Dy(phiGF,a,i,j)
                S_y        = Dy(SGF,  a,i,j)
                Fx_y       = Dy(FxGF, a,i,j)
                Fy_y       = Dy(FyGF, a,i,j)

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

                B1_xx      = Dxx(B1GF, a,i,j)
                B2_xx      = Dxx(B2GF, a,i,j)
                G_xx       = Dxx(GGF,  a,i,j)
                phi_xx     = Dxx(phiGF,a,i,j)
                S_xx       = Dxx(SGF,  a,i,j)

                B1_yy      = Dyy(B1GF, a,i,j)
                B2_yy      = Dyy(B2GF, a,i,j)
                G_yy       = Dyy(GGF,  a,i,j)
                phi_yy     = Dyy(phiGF,a,i,j)
                S_yy       = Dyy(SGF,  a,i,j)

                B2_xy      = Dx(Dy, B2GF, a,i,j)
                G_xy       = Dx(Dy, GGF,  a,i,j)
                phi_xy     = Dx(Dy, phiGF,a,i,j)
                S_xy       = Dx(Dy, SGF,  a,i,j)

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
                @inbounds @simd for aa in 1:Nu
                    aux.A_mat[a,aa] = aux.ABCS[1] * Duu[a,aa] + aux.ABCS[2] * Du[a,aa]
                end
                aux.A_mat[a,a] += aux.ABCS[3]
            end

            # BC (first order equation)

            aux.b_vec[1]    = bc.phid[i,j]
            aux.A_mat[1,:] .= 0.0
            aux.A_mat[1,1]  = 1.0

            solve_lin_system!(aux.A_mat, aux.b_vec)

            @inbounds @simd for aa in 1:Nu
                phidGF[aa,i,j] = aux.b_vec[aa]
            end

        end
    end

    nothing
end

function solve_A!(bulk::Bulk, bc::BC, gauge::Gauge, deriv::BulkDeriv, aux_acc,
                  sys::System, evoleq::AffineNull)
    Nu, Nx, Ny = size(sys)

    B1GF    = getB1(bulk)
    B2GF    = getB2(bulk)
    GGF     = getG(bulk)
    phiGF   = getphi(bulk)
    SGF     = getS(bulk)
    FxGF    = getFx(bulk)
    FyGF    = getFy(bulk)
    SdGF    = getSd(bulk)
    B1dGF   = getB1d(bulk)
    B2dGF   = getB2d(bulk)
    GdGF    = getGd(bulk)
    phidGF  = getphid(bulk)
    AGF     = getA(bulk)
    xiGF    = getxi(gauge)

    Du_B1   = deriv.Du_B1
    Du_B2   = deriv.Du_B2
    Du_G    = deriv.Du_G
    Du_phi  = deriv.Du_phi
    Du_S    = deriv.Du_S
    Du_Fx   = deriv.Du_Fx
    Du_Fy   = deriv.Du_Fy
    Duu_B1  = deriv.Duu_B1
    Duu_B2  = deriv.Duu_B2
    Duu_G   = deriv.Duu_G
    Duu_phi = deriv.Duu_phi
    Duu_S   = deriv.Duu_S
    Duu_Fx  = deriv.Duu_Fx
    Duu_Fy  = deriv.Duu_Fy

    Du  = sys.Du
    Duu = sys.Duu
    Dx  = sys.Dx
    Dxx = sys.Dxx
    Dy  = sys.Dy
    Dyy = sys.Dyy

    potential = evoleq.potential
    phi0      = evoleq.phi0

    @fastmath @inbounds @threads for j in 1:Ny
        @inbounds for i in 1:Nx
            id  = Threads.threadid()
            aux = aux_acc[id]

            xi    = xiGF[1,i,j]
            xi_x  = Dx(xiGF, 1,i,j)
            xi_y  = Dy(xiGF, 1,i,j)
            xi_xx = Dxx(xiGF, 1,i,j)
            xi_yy = Dyy(xiGF, 1,i,j)
            xi_xy = Dx(Dy, xiGF, 1,i,j)

            @inbounds @simd for a in 1:Nu
                u     = sys.ucoord[a]
                u2    = u * u
                u3    = u * u2
                u4    = u2 * u2

                B1    = B1GF[a,i,j]
                B2    = B2GF[a,i,j]
                G     = GGF[a,i,j]
                phi   = phiGF[a,i,j]
                S     = SGF[a,i,j]
                Fx    = FxGF[a,i,j]
                Fy    = FyGF[a,i,j]
                Sd    = SdGF[a,i,j]
                B1d   = B1dGF[a,i,j]
                B2d   = B2dGF[a,i,j]
                Gd    = GdGF[a,i,j]
                phid  = phidGF[a,i,j]

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

                B1_x       = Dx(B1GF, a,i,j)
                B2_x       = Dx(B2GF, a,i,j)
                G_x        = Dx(GGF,  a,i,j)
                phi_x      = Dx(phiGF,a,i,j)
                S_x        = Dx(SGF,  a,i,j)
                Fx_x       = Dx(FxGF, a,i,j)
                Fy_x       = Dx(FyGF, a,i,j)

                B1_y       = Dy(B1GF, a,i,j)
                B2_y       = Dy(B2GF, a,i,j)
                G_y        = Dy(GGF,  a,i,j)
                phi_y      = Dy(phiGF,a,i,j)
                S_y        = Dy(SGF,  a,i,j)
                Fx_y       = Dy(FxGF, a,i,j)
                Fy_y       = Dy(FyGF, a,i,j)

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

                B1_xx      = Dxx(B1GF, a,i,j)
                B2_xx      = Dxx(B2GF, a,i,j)
                G_xx       = Dxx(GGF,  a,i,j)
                phi_xx     = Dxx(phiGF,a,i,j)
                S_xx       = Dxx(SGF,  a,i,j)

                B1_yy      = Dyy(B1GF, a,i,j)
                B2_yy      = Dyy(B2GF, a,i,j)
                G_yy       = Dyy(GGF,  a,i,j)
                phi_yy     = Dyy(phiGF,a,i,j)
                S_yy       = Dyy(SGF,  a,i,j)

                B2_xy      = Dx(Dy, B2GF, a,i,j)
                G_xy       = Dx(Dy, GGF,  a,i,j)
                phi_xy     = Dx(Dy, phiGF,a,i,j)
                S_xy       = Dx(Dy, SGF,  a,i,j)

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
                @inbounds @simd for aa in 1:Nu
                    aux.A_mat[a,aa] = aux.ABCS[1] * Duu[a,aa] + aux.ABCS[2] * Du[a,aa]
                end
                aux.A_mat[a,a] += aux.ABCS[3]
            end

            # BC

            aux.b_vec[1]    = bc.A[i,j]
            aux.A_mat[1,:] .= 0.0
            aux.A_mat[1,1]  = 1.0

            aux.b_vec[end]    = bc.A_u[i,j]
            @inbounds @simd for aa in 1:Nu
                aux.A_mat[end,aa]  = Du[1,aa]
            end

            solve_lin_system!(aux.A_mat, aux.b_vec)

            @inbounds @simd for aa in 1:Nu
                AGF[aa,i,j] = aux.b_vec[aa]
            end

        end
    end

    nothing
end

function solve_nested!(bulkconstrain::BulkConstrained, bulkevol::BulkEvolved, bc::BC,
                       gauge::Gauge, deriv::BulkDeriv, aux_acc,
                       sys::System, evoleq::AffineNull)
    Du_B1   = deriv.Du_B1
    Du_B2   = deriv.Du_B2
    Du_G    = deriv.Du_G
    Du_phi  = deriv.Du_phi
    Du_S    = deriv.Du_S
    Du_Fx   = deriv.Du_Fx
    Du_Fy   = deriv.Du_Fy
    Du_A    = deriv.Du_A
    Duu_B1  = deriv.Duu_B1
    Duu_B2  = deriv.Duu_B2
    Duu_G   = deriv.Duu_G
    Duu_phi = deriv.Duu_phi
    Duu_S   = deriv.Duu_S
    Duu_Fx  = deriv.Duu_Fx
    Duu_Fy  = deriv.Duu_Fy
    Duu_A   = deriv.Duu_A

    Du  = sys.Du
    Duu = sys.Duu
    Dx  = sys.Dx
    Dxx = sys.Dxx
    Dy  = sys.Dy
    Dyy = sys.Dyy

    bulk = Bulk(bulkevol, bulkconstrain)

    SGF     = getS(bulk)
    FxGF    = getFx(bulk)
    FyGF    = getFy(bulk)
    AGF     = getA(bulk)

    # solve for S
    solve_S!(bulk, bc, gauge, deriv, aux_acc, sys, evoleq)

    # take u-derivatives of S
    @sync begin
        @spawn mul!(Du_S,   Du,  SGF)
        @spawn mul!(Duu_S,  Duu, SGF)
    end

    # solve for Fx and Fy
    solve_Fxy!(bulk, bc, gauge, deriv, aux_acc, sys, evoleq)

    # take u-derivatives of Fx and Fy
    @sync begin
        @spawn mul!(Du_Fx,   Du,  FxGF)
        @spawn mul!(Du_Fy,   Du,  FyGF)
        @spawn mul!(Duu_Fx,  Duu, FxGF)
        @spawn mul!(Duu_Fy,  Duu, FyGF)
    end

    # solve for Sd
    solve_Sd!(bulk, bc, gauge, deriv, aux_acc, sys, evoleq)

    # solving for B2d, (B1d,Gd) and phid are independent processes. we can
    # therefore @spawn, here
    @sync begin
        @spawn solve_B2d!(bulk, bc, gauge, deriv, aux_acc, sys, evoleq)
        @spawn solve_B1dGd!(bulk, bc, gauge, deriv, aux_acc, sys, evoleq)
        @spawn solve_phid!(bulk, bc, gauge, deriv, aux_acc, sys, evoleq)
    end

    # solve for A
    solve_A!(bulk, bc, gauge, deriv, aux_acc, sys, evoleq)

    # take u-derivatives of A. they will be needed for syncing the domains and
    # also for the xi_t function
    mul!(Du_A,   Du,  AGF)
    mul!(Duu_A,  Duu, AGF)

    nothing
end


function syncBCs!(bc::BC, bulk::BulkConstrained, deriv::BulkDeriv)
    SGF     = getS(bulk)
    FxGF    = getFx(bulk)
    FyGF    = getFy(bulk)
    SdGF    = getSd(bulk)
    B1dGF   = getB1d(bulk)
    B2dGF   = getB2d(bulk)
    GdGF    = getGd(bulk)
    phidGF  = getphid(bulk)
    AGF     = getA(bulk)

    Nu, Nx, Ny = size(SGF)

    @fastmath @inbounds @threads for j in 1:Ny
        @inbounds @simd for i in 1:Nx
            bc.S[i,j]    = SGF[end,i,j]
            bc.Fx[i,j]   = FxGF[end,i,j]
            bc.Fy[i,j]   = FyGF[end,i,j]
            bc.Sd[i,j]   = SdGF[end,i,j]
            bc.B1d[i,j]  = B1dGF[end,i,j]
            bc.B2d[i,j]  = B2dGF[end,i,j]
            bc.Gd[i,j]   = GdGF[end,i,j]
            bc.phid[i,j] = phidGF[end,i,j]
            bc.A[i,j]    = AGF[end,i,j]

            bc.S_u[i,j]   = deriv.Du_S[end,i,j]
            bc.Fx_u[i,j]  = deriv.Du_Fx[end,i,j]
            bc.Fy_u[i,j]  = deriv.Du_Fy[end,i,j]
            bc.A_u[i,j]   = deriv.Du_A[end,i,j]
        end
    end

    nothing
end

function set_innerBCs!(bc::BC, bulk::BulkEvolved, boundary::Boundary,
                       gauge::Gauge, deriv::BulkDeriv, sys::System{Inner},
                       evoleq::AffineNull)
    _, Nx, Ny = size(sys)

    B1GF    = getB1(bulk)
    B2GF    = getB2(bulk)
    GGF     = getG(bulk)
    phiGF   = getphi(bulk)
    a4GF    = geta4(boundary)
    fx2GF   = getfx2(boundary)
    fy2GF   = getfy2(boundary)
    xiGF    = getxi(gauge)

    Dx  = sys.Dx
    Dy  = sys.Dy

    phi0  = evoleq.phi0
    phi02 = phi0 * phi0
    phi03 = phi0 * phi02
    phi04 = phi02 * phi02

    @fastmath @inbounds @threads for j in 1:Ny
        @inbounds @simd for i in 1:Nx
            xi      = xiGF[1,i,j]
            xi3     = xi*xi*xi

            phi     = phiGF[1,i,j]
            phi_u   = deriv.Du_phi[1,i,j]

            xi_x    = Dx(xiGF, 1,i,j)
            xi_y    = Dy(xiGF, 1,i,j)

            b14     = B1GF[1,i,j]
            b24     = B2GF[1,i,j]
            g4      = GGF[1,i,j]

            a4      = a4GF[1,i,j]

            fx2     = fx2GF[1,i,j]
            fy2     = fy2GF[1,i,j]

            fx2_x   = Dx(fx2GF, 1,i,j)
            fy2_y   = Dy(fy2GF, 1,i,j)

            b14_x   = Dx(B1GF,1,i,j)
            b24_x   = Dx(B2GF,1,i,j)
            phi_x   = Dx(phiGF,1,i,j)
            g4_x    = Dx(GGF,1,i,j)

            b14_y   = Dy(B1GF,1,i,j)
            b24_y   = Dy(B2GF,1,i,j)
            phi_y   = Dy(phiGF,1,i,j)
            g4_y    = Dy(GGF,1,i,j)

            phi2    = phi03 * phi - phi0 * xi * xi
            phi2_x  = phi03 * phi_x - 2 * phi0 * xi * xi_x
            phi2_y  = phi03 * phi_y - 2 * phi0 * xi * xi_y

            bc.S[i,j]   = phi04 * (1 - 18 * phi) / 54
            bc.S_u[i,j] = phi02 * (-12 * xi3 +
                                   phi02 * xi * (18 * phi - 5) -
                                   24 * phi02 * phi_u) / 90

            bc.Fx[i,j]   = fx2
            bc.Fx_u[i,j] = -2 * fx2 * xi - 12 / 15 * (b14_x + b24_x - g4_y) +
                4/15 * phi0 * phi2_x

            bc.Fy[i,j]   = fy2
            bc.Fy_u[i,j] = -2 * fy2 * xi - 12 / 15 * (-b14_y + b24_y - g4_x) +
                4/15 * phi0 * phi2_y

            bc.Sd[i,j] = a4 / 2 - 5 * phi04/36 + phi0 * phi2/2

            bc.B2d[i,j] = -2 * b24
            bc.B1d[i,j] = -2 * b14
            bc.Gd[i,j]  = -2 * g4

            bc.A[i,j]   = a4
            bc.A_u[i,j] = -2 * xi * a4 - 2 * phi0 * xi * phi2 -
                2/3 * (phi02 * xi3 + fx2_x + fy2_y + phi04 * phi_u)
        end
    end

    # separate the phid case, since the if statement may prevent loop
    # vectorization

    if abs(phi0) < 1e-9
        fill!(bc.phid, 0)
        return
    end

    @fastmath @inbounds for j in 1:Ny
        @inbounds @simd for i in 1:Nx
            xi      = xiGF[1,i,j]
            phi     = phiGF[1,i,j]
            phi2    = phi03 * phi - phi0 * xi * xi

            bc.phid[i,j] = 1/3 - 3/2 * phi2 / phi03
        end
    end

    nothing
end

function set_outerBCs!(bc::BC, bulk::BulkConstrained, gauge::Gauge,
                       deriv::BulkDeriv, sys::System, evoleq::AffineNull)
    _, Nx, Ny = size(sys)

    SGF     = getS(bulk)
    FxGF    = getFx(bulk)
    FyGF    = getFy(bulk)
    SdGF    = getSd(bulk)
    B1dGF   = getB1d(bulk)
    B2dGF   = getB2d(bulk)
    GdGF    = getGd(bulk)
    phidGF  = getphid(bulk)
    AGF     = getA(bulk)
    xiGF    = getxi(gauge)

    phi0 = evoleq.phi0

    # we are here assuming that the inner and outer grids merely touch at the
    # interface, so we pass the values at this point without any interpolation
    u0 = sys.ucoord[end]

    @fastmath @inbounds @threads for j in 1:Ny
        @inbounds @simd for i in 1:Nx
            xi     = xiGF[1,i,j]

            S      = SGF[end,i,j]
            Fx     = FxGF[end,i,j]
            Fy     = FyGF[end,i,j]
            Sd     = SdGF[end,i,j]
            B1d    = B1dGF[end,i,j]
            B2d    = B2dGF[end,i,j]
            Gd     = GdGF[end,i,j]
            phid   = phidGF[end,i,j]
            A      = AGF[end,i,j]

            S_u    = deriv.Du_S[end,i,j]
            Fx_u   = deriv.Du_Fx[end,i,j]
            Fy_u   = deriv.Du_Fy[end,i,j]
            A_u    = deriv.Du_A[end,i,j]

            bc.S[i,j]   = S_inner_to_outer(S, u0, xi, phi0)
            bc.S_u[i,j] = S_u_inner_to_outer(S_u, S, u0, xi, phi0)

            bc.Fx[i,j]   = F_inner_to_outer(Fx, u0)
            bc.Fy[i,j]   = F_inner_to_outer(Fy, u0)
            bc.Fx_u[i,j] = F_u_inner_to_outer(Fx_u, Fx, u0)
            bc.Fy_u[i,j] = F_u_inner_to_outer(Fy_u, Fy, u0)

            bc.Sd[i,j]  = Sd_inner_to_outer(Sd, u0, xi, phi0)

            # B1d, B2d, and Gd transform in the same way
            bc.B2d[i,j] = Bd_inner_to_outer(B2d, u0)
            bc.B1d[i,j] = Bd_inner_to_outer(B1d, u0)
            bc.Gd[i,j]  = Bd_inner_to_outer(Gd, u0)

            bc.phid[i,j] = phid_inner_to_outer(phid, u0, phi0)

            bc.A[i,j]   = A_inner_to_outer(A, u0, xi, phi0)
            bc.A_u[i,j] = A_u_inner_to_outer(A_u, A, u0, xi, phi0)
        end
    end

    nothing
end


# We assume that the first entry on these arrays is the inner grid, and that
# there is only one domain spanning this grid. If we ever change this
# construction we must remember to make the appropriate changes here.
function solve_nesteds!(bulkconstrains, bulkevols, boundary::Boundary, gauge::Gauge,
                        bcs, derivs, aux_accs,
                        systems::SystemPartition, evoleq::AffineNull)
    Nsys = length(systems)

    # take all u-derivatives of the bulkevols functions
    @sync begin
        @inbounds for i in 1:Nsys
            deriv  = derivs[i]
            sys    = systems[i]
            bulk   = bulkevols[i]
            B1GF   = getB1(bulk)
            B2GF   = getB2(bulk)
            GGF    = getG(bulk)
            phiGF  = getphi(bulk)

            @spawn mul!(deriv.Du_B1,  sys.Du,  B1GF)
            @spawn mul!(deriv.Du_B2,  sys.Du,  B2GF)
            @spawn mul!(deriv.Du_G,   sys.Du,  GGF)
            @spawn mul!(deriv.Du_phi, sys.Du,  phiGF)
            @spawn mul!(deriv.Duu_B1, sys.Duu, B1GF)
            @spawn mul!(deriv.Duu_B2, sys.Duu, B2GF)
            @spawn mul!(deriv.Duu_G,  sys.Duu, GGF)
            @spawn mul!(deriv.Duu_phi,sys.Duu, phiGF)
        end
    end

    set_innerBCs!(bcs[1], bulkevols[1], boundary, gauge, derivs[1], systems[1], evoleq)

    solve_nested!(bulkconstrains[1], bulkevols[1], bcs[1], gauge,
                  derivs[1], aux_accs[1], systems[1], evoleq)

    set_outerBCs!(bcs[2], bulkconstrains[1], gauge, derivs[1], systems[1], evoleq)

    @inbounds for i in 2:Nsys-1
        solve_nested!(bulkconstrains[i], bulkevols[i], bcs[i], gauge,
                      derivs[i], aux_accs[i], systems[i], evoleq)
        syncBCs!(bcs[i+1], bulkconstrains[i], derivs[i])
    end
    solve_nested!(bulkconstrains[Nsys], bulkevols[Nsys], bcs[Nsys], gauge,
                  derivs[Nsys], aux_accs[Nsys], systems[Nsys], evoleq)

    nothing
end

# for testing only: set all constrained variables to zero
function solve_nesteds!(bulkconstrains, bulkevols, boundary::Boundary, gauge::Gauge,
                        bcs, derivs, aux_accs,
                        systems, evoleq::EvolTest0)
    Nsys = length(systems)

    for i in 1:Nsys
        fill!(bulkconstrains[i], 0)
    end
    nothing
end


struct Nested{S,C,D,B,A}
    systems        :: S
    bulkconstrains :: C
    derivs         :: D
    bcs            :: B
    aux_accs       :: A
end

function Nested(systems::SystemPartition, bulkconstrains::BulkPartition,
                derivs::BulkPartition)
    T        = Jecco.coord_eltype(systems[1].ucoord)
    bcs      = BCs(systems)
    aux_accs = [Aux{T}(sys) for sys in systems]

    Nested{typeof(systems),typeof(bulkconstrains),typeof(derivs),
           typeof(bcs),typeof(aux_accs)}(systems, bulkconstrains, derivs, bcs,
                                         aux_accs)
end

function Nested(systems::SystemPartition, bulkconstrains::BulkPartition)
    T        = Jecco.coord_eltype(systems[1].ucoord)
    bcs      = BCs(systems)
    aux_accs = [Aux{T}(sys) for sys in systems]
    derivs   = BulkDerivPartition(systems)

    Nested(systems, bulkconstrains, derivs)
end

function Nested(systems::SystemPartition)
    T        = Jecco.coord_eltype(systems[1].ucoord)
    bcs      = BCs(systems)
    aux_accs = [Aux{T}(sys) for sys in systems]

    bulkconstrains = BulkConstrainedPartition(systems)
    derivs         = BulkDerivPartition(systems)

    Nested(systems, bulkconstrains, derivs)
end

function (nested::Nested)(bulkevols::BulkPartition, boundary::Boundary,
                          gauge::Gauge, evoleq::EvolutionEquations)
    solve_nesteds!(nested.bulkconstrains, bulkevols, boundary, gauge, nested.bcs,
                   nested.derivs, nested.aux_accs, nested.systems, evoleq)
end

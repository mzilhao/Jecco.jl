
using Jecco
using Jecco.AdS5_3_1

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
)


systems = Jecco.AdS5_3_1.create_systems(par_grid)

sys = systems[1]

Nu, Nx, Ny = size(sys.grid)
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


uu   = nested.uu
xx   = nested.xx
yy   = nested.yy

aux_acc = nested.aux_acc

Du  = sys.Du
Duu = sys.Duu
Dx  = sys.Dx
Dxx = sys.Dxx
Dy  = sys.Dy
Dyy = sys.Dyy

# mul!(Du_B1,  Du, bulk.B1)
# mul!(Du_B2,  Du, bulk.B2)
# mul!(Du_G,   Du, bulk.G)
# mul!(Du_phi, Du, bulk.phi)


AA = zeros(2,2)
BB = zeros(2,2)
CC = zeros(2,2)
SS = zeros(2)


j = 4
i = 8

a = 2

# id  = Threads.threadid()
id  = 1
# aux = aux_acc[id]

vars = Jecco.AdS5_3_1.FxyVars{Float64}()


u          = uu[a]
u2         = u * u
u3         = u * u2
u4         = u2 * u2

vars.u     = u

# maybe it's worth to make some structs (or macros), here...

vars.B1    = bulk.B1[a,i,j]
vars.B1p   = -u2 * Du(bulk.B1, a,i,j)
vars.B1_x  = Dx(bulk.B1, a,i,j)
vars.B1_y  = Dy(bulk.B1, a,i,j)
vars.B1pp  = 2*u3 * Du(bulk.B1, a,i,j) + u4 * Duu(bulk.B1, a,i,j)
vars.B1p_x = -u2 * Du(Dx, bulk.B1, a,i,j)
vars.B1p_y = -u2 * Du(Dy, bulk.B1, a,i,j)

vars.B2    = bulk.B2[a,i,j]
vars.B2p   = -u2 * Du(bulk.B2, a,i,j)
vars.B2_x  = Dx(bulk.B2, a,i,j)
vars.B2_y  = Dy(bulk.B2, a,i,j)
vars.B2pp  = 2*u3 * Du(bulk.B2, a,i,j) + u4 * Duu(bulk.B2, a,i,j)
vars.B2p_x = -u2 * Du(Dx, bulk.B2, a,i,j)
vars.B2p_y = -u2 * Du(Dy, bulk.B2, a,i,j)

vars.G    = bulk.G[a,i,j]
vars.Gp   = -u2 * Du(bulk.G, a,i,j)
vars.G_x  = Dx(bulk.G, a,i,j)
vars.G_y  = Dy(bulk.G, a,i,j)
vars.Gpp  = 2*u3 * Du(bulk.G, a,i,j) + u4 * Duu(bulk.G, a,i,j)
vars.Gp_x = -u2 * Du(Dx, bulk.G, a,i,j)
vars.Gp_y = -u2 * Du(Dy, bulk.G, a,i,j)

vars.phi    = bulk.phi[a,i,j]
vars.phip   = -u2 * Du(bulk.phi, a,i,j)
vars.phi_x  = Dx(bulk.phi, a,i,j)
vars.phi_y  = Dy(bulk.phi, a,i,j)
# vars.phipp  = 2*u3 * Du(bulk.phi, a,i,j) + u4 * Duu(bulk.phi, a,i,j)
# vars.phip_x = -u2 * Du(Dx, bulk.phi, a,i,j)
# vars.phip_y = -u2 * Du(Dy, bulk.phi, a,i,j)

vars.S    = bulk.S[a,i,j]
vars.Sp   = -u2 * Du(bulk.S, a,i,j)
vars.S_x  = Dx(bulk.S, a,i,j)
vars.S_y  = Dy(bulk.S, a,i,j)
vars.Spp  = 2*u3 * Du(bulk.S, a,i,j) + u4 * Duu(bulk.S, a,i,j)
vars.Sp_x = -u2 * Du(Dx, bulk.S, a,i,j)
vars.Sp_y = -u2 * Du(Dy, bulk.S, a,i,j)

Jecco.AdS5_3_1.Fxy_outer_eq_coeff!(AA, BB, CC, SS, vars)


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

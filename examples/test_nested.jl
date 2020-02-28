
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

mul!(Du_B1,  Du, bulk.B1)
mul!(Du_B2,  Du, bulk.B2)
mul!(Du_G,   Du, bulk.G)
mul!(Du_phi, Du, bulk.phi)


AA = zeros(4,4)
BB = zeros(4,4)
CC = zeros(4,4)
SS = zeros(2)


j = 4
i = 8

a = 2

# id  = Threads.threadid()
id  = 1
# aux = aux_acc[id]

vars = Jecco.AdS5_3_1.FxyVars{Float64}()


u              = uu[a]
vars.u     = u

vars.B1    = bulk.B1[a,i,j]
vars.B1p   = -u*u * Du_B1[a,i,j]

# vars.B2    = bulk.B2[a,i,j]
vars.B2p   = -u*u * Du_B2[a,i,j]

vars.G     = bulk.G[a,i,j]

vars.G     = bulk.G[a,i,j]
vars.Gp    = -u*u * Du_G[a,i,j]

vars.phip  = -u*u * Du_phi[a,i,j]

vars.S     = bulk.S[a,i,j]

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

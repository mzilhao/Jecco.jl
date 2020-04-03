
using Jecco
using LinearAlgebra
using SparseArrays

sine2D(x, y, kx, ky) = sin(kx * x) * sin(ky * y)

source1(x, y, kx, ky) = -(kx * kx + ky * ky) * sin(kx * x) * sin(ky * y)

source2(x, y, kx, ky) = 1 + 2 * kx * ky * cos(kx * x) * cos(ky * y) +
    ky * sin(kx * x) * cos(ky * y) + kx * cos(kx * x) * sin(ky * y) -
    kx * kx * sin(kx * x) * sin(ky * y) - ky * ky * sin(kx * x) * sin(ky * y)

function getind2D(Nx::Int)
    (i,j) -> i + (j-1)*Nx
end


function build_operators(hx, hy, Nx::Int, Ny::Int)
    ord = 4

    Dx_op  = CenteredDiff{1}(1, ord, hx, Nx)
    Dxx_op = CenteredDiff{1}(2, ord, hx, Nx)

    Dy_op  = CenteredDiff{2}(1, ord, hy, Ny)
    Dyy_op = CenteredDiff{2}(2, ord, hy, Ny)

    Dx  = kron(I(Ny), SparseMatrixCSC(Dx_op))
    Dxx = kron(I(Ny), SparseMatrixCSC(Dxx_op))
    Dy  = kron(SparseMatrixCSC(Dy_op), I(Nx))
    Dyy = kron(SparseMatrixCSC(Dyy_op), I(Nx))

    Dx, Dy, Dxx, Dyy
end


x_min            = -5.0
x_max            =  5.0
x_nodes          =  128
# x_nodes          =  12
y_min            = -5.0
y_max            =  5.0
y_nodes          =  64
# y_nodes          =  4

nx  = 2
ny  = 4

xcoord  = CartCoord{1}("x", x_min, x_max, x_nodes, endpoint=false)
ycoord  = CartCoord{2}("y", y_min, y_max, y_nodes, endpoint=false)

hx = Jecco.delta(xcoord)
hy = Jecco.delta(ycoord)

Nx = xcoord.nodes
Ny = ycoord.nodes

Dx, Dy, Dxx, Dyy = build_operators(hx, hy, Nx, Ny)


kx = 2*pi * nx / (x_max-x_min)
ky = 2*pi * ny / (y_max-y_min)

f_exact = [sine2D(xcoord[i], ycoord[j], kx, ky)
           for i in 1:Nx, j in 1:Ny]

f0 = zeros(Nx,Ny)

# ind2D = getind2D(Nx)
ind2D = LinearIndices(f0)

M = Nx * Ny

A_mat = spzeros(M,M)
b_vec = zeros(M)

for j in 1:Ny, i in 1:Nx
    # idx = ind2D(i,j)
    idx = ind2D[i,j]

    xi  = xcoord[i]
    yi  = ycoord[j]

    b_vec[idx] = source1(xi, yi, kx, ky)
end


for kl in 1:M, ij in 1:M
    A_mat[ij,kl] = Dxx[ij,kl] + Dyy[ij,kl]
end

# f_test = similar(f_exact)
# f_xx   = similar(f_exact)
# f_yy   = similar(f_exact)

# tmp1 = Dxx * f_exact[:]
# tmp2 = Dyy * f_exact[:]
# lap = A_mat * f_exact[:]
# for idx in eachindex(tmp)
#     f_xx[idx] = tmp1[idx]
#     f_yy[idx] = tmp2[idx]
#     f_test[idx] = lap[idx]
# end




sol = A_mat \ b_vec

for idx in eachindex(f0)
    f0[idx] = sol[idx]
end


############### FIM #####################

# Axx = ones(Nx,Ny)
# Axy = ones(Nx,Ny)
# Ayy = ones(Nx,Ny)

# Bx  = ones(Nx,Ny)
# By  = ones(Nx,Ny)

# C   = ones(Nx,Ny)


# for i in 1:Nx, j in 1:Ny
#     idx1= ind2D(i,j)

#     axx = Axx[i,j]
#     axy = Axy[i,j]
#     ayy = Ayy[i,j]
#     bx  =  Bx[i,j]
#     by  =  By[i,j]
#     c   =   C[i,j]

#     xi = xcoord[i]
#     yi = ycoord[j]

#     b_vec[idx1] = source2(xi, yi, kx, ky)

#     for ii in 1:Nx, jj in 1:Ny
#         idx2 = ind2D(ii,jj)
#         A_mat[idx1,idx2] = axx *  Dxx[idx1,idx2] + 2 * axy * Dx[idx1,idx2] * Dy[idx1,idx2] + ayy * Dyy[idx1,idx2] +
#             bx * Dx[idx1,idx2] + by * Dy[idx1,idx2] + c * I[idx1,idx2]
#     end
# end


# sol = A_mat \ b_vec



# C_mat = zeros(Nx, Nx)
# # B_mat = zeros(Nx, Nx)

# α  = ones(Nx)
# # α[1:2:Nx] .= 0.0

# for i in 1:Nx
#     a0 = α[i]
#     for j in 1:Nx
#         C_mat[i,j] = a0 * Dxx[i,j]
#     end
#     C_mat[i,i] += 1.0
# end


# # for j in 1:Nx, i in 1:Nx
# #     a0 = α[i]
# #     B_mat[i,j] = a0 * Dx[i,j] + 1.0*I[i,j]
# # end


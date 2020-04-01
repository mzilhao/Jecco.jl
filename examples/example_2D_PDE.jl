
using Jecco
using LinearAlgebra

sine2D(x, y, kx, ky) = sin(kx * x) * sin(ky * y)

source(x, y, kx, ky) = 1 + 2 * kx * ky * cos(kx * x) * cos(ky * y) +
    ky * sin(kx * x) * cos(ky * y) + kx * cos(kx * x) * sin(ky * y) -
    kx * kx * sin(kx * x) * sin(ky * y) - ky * ky * sin(kx * x) * sin(ky * y)

function getind2D(Nx::Int)
    (i,j) -> i + (j-1)*Nx
end

x_min            = -5.0
x_max            =  5.0
# x_nodes          =  128
x_nodes          =  16
y_min            = -5.0
y_max            =  5.0
# y_nodes          =  64
y_nodes          =  4

ord = 4

nx  = 2
ny  = 4


xcoord  = CartCoord{1}("x", x_min, x_max, x_nodes, endpoint=false)
ycoord  = CartCoord{2}("y", y_min, y_max, y_nodes, endpoint=false)

Dx  = CenteredDiff{1}(1, ord, Jecco.delta(xcoord), xcoord.nodes)
Dxx = CenteredDiff{1}(2, ord, Jecco.delta(xcoord), xcoord.nodes)

Dy  = CenteredDiff{2}(1, ord, Jecco.delta(ycoord), ycoord.nodes)
Dyy = CenteredDiff{2}(2, ord, Jecco.delta(ycoord), ycoord.nodes)


kx = 2*pi * nx / (x_max-x_min)
ky = 2*pi * ny / (y_max-y_min)

Nx = xcoord.nodes
Ny = ycoord.nodes

f_exact = [sine2D(xcoord[i], ycoord[j], kx, ky)
           for i in 1:Nx, j in 1:Ny]

f0 = zeros(Nx,Ny)

ind2D = getind2D(Nx)
# ind2D = LinearIndices(f0)

A_mat = zeros(Nx * Ny, Nx * Ny)
b_vec = zeros(Nx * Ny)

Axx = ones(Nx,Ny)
Axy = ones(Nx,Ny)
Ayy = ones(Nx,Ny)

Bx  = ones(Nx,Ny)
By  = ones(Nx,Ny)

C   = ones(Nx,Ny)



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


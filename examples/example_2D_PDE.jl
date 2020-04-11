
using Jecco
using LinearAlgebra
using SparseArrays
using Plots
source(x,y) = exp(-x^2 - y^2) * (-4 + 3 * (x^2 + y^2) + 4 * x^2 * y^2)

# https://en.wikipedia.org/wiki/Kronecker_product
# https://arxiv.org/pdf/1801.01483.pdf
function deriv_operators(hx, hy, Nx::Int, Ny::Int, ord::Int)
    Dx_op  = CenteredDiff{1}(1, ord, hx, Nx)
    Dxx_op = CenteredDiff{1}(2, ord, hx, Nx)

    Dy_op  = CenteredDiff{2}(1, ord, hy, Ny)
    Dyy_op = CenteredDiff{2}(2, ord, hy, Ny)

    Dx  = kron(I(Ny), SparseMatrixCSC(Dx_op))
    Dxx = kron(I(Ny), SparseMatrixCSC(Dxx_op))
    Dy  = kron(SparseMatrixCSC(Dy_op), I(Nx))
    Dyy = kron(SparseMatrixCSC(Dyy_op), I(Nx))

    Dx, Dy, Dxx, Dyy, Dx * Dy
end


x_min            = -5.0
x_max            =  5.0
x_nodes          =  128
# x_nodes          =  12
y_min            = -5.0
y_max            =  5.0
y_nodes          =  64
# y_nodes          =  4


xcoord  = CartCoord{1}("x", x_min, x_max, x_nodes, endpoint=false)
ycoord  = CartCoord{2}("y", y_min, y_max, y_nodes, endpoint=false)

hx = Jecco.delta(xcoord)
hy = Jecco.delta(ycoord)

Nx = xcoord.nodes
Ny = ycoord.nodes

f_exact = [exp(-xcoord[i]^2 - ycoord[j]^2) for i in 1:Nx, j in 1:Ny]


Dx, Dy, Dxx, Dyy, Dxy = deriv_operators(hx, hy, Nx, Ny, ord)

f0    = zeros(Nx,Ny)
ind2D = LinearIndices(f0)

M = Nx * Ny

# A_mat = spzeros(M,M)
b_vec = zeros(M)

axx     = ones(M)
ayy     = ones(M)
axy     = zeros(M)
bx      = zeros(M)
by      = zeros(M)
cc      = zeros(M)


for j in 1:Ny, i in 1:Nx
    idx = ind2D[i,j]

    xi  = xcoord[i]
    yi  = ycoord[j]

    axy[idx] = xi * yi
    bx[idx]  = xi
    by[idx]  = yi
    cc[idx]  = xi^2 + yi^2

    b_vec[idx] = source(xi, yi)
end

# TODO: improve this construction
A_mat = axx .* Dxx + ayy .* Dyy + axy .* Dxy + bx .* Dx + by .* Dy + Diagonal(cc)

# for kl in 1:M, ij in 1:M
#     A_mat[ij,kl] = Dxx[ij,kl] + Dyy[ij,kl]
# end


# better to do: broadcast!(*, A_mat, Axy, Dxy)
# for kl in 1:M, ij in 1:M
#     axy = Axy[ij]
#     A_mat[ij,kl] = axy * Dxy[ij,kl]
# end


A_fact = factorize(A_mat)
sol    = A_fact \ b_vec

@inbounds for idx in eachindex(f0)
    f0[idx] = sol[idx]
end

j_slice = div(Ny,2)
x = xcoord[:]

plot(x, f_exact[:,j_slice])
scatter!(x, f0[:,j_slice])

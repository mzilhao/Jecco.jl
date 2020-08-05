
using Jecco
using LinearAlgebra
using SparseArrays
using Plots

#= solving equation

 f_xx + f_yy + x y f_xy + x f_x + y f_y + (x^2 + y^2) f
  + (f_x)^2 + (f_y)^2 + x y f_x f_y = Source(x,y)

for the source below, the solution should be given by

 f = exp(-x^2 - y^2)
=#

source(x,y) = exp(-x^2 - y^2) * (-4 + 3 * (x^2 + y^2) + 4 * x^2 * y^2) +
    4 * exp(-2*x^2 - 2*y^2) * (x^2 + y^2 + x^2 * y^2)

function compute_residual!(res::Array, f::Array, x::Vector, y::Vector,
                           Dxx_op::Jecco.FiniteDiffDeriv, Dyy_op::Jecco.FiniteDiffDeriv,
                           Dx_op::Jecco.FiniteDiffDeriv, Dy_op::Jecco.FiniteDiffDeriv)
    Nx = length(x)
    Ny = length(y)

    f_xx = Dxx_op * f
    f_yy = Dyy_op * f
    f_x  = Dx_op  * f
    f_y  = Dy_op  * f
    f_xy = Dx_op  * f_y

    for j in 1:Ny, i in 1:Nx
        x1 = x[i]
        y1 = y[j]

        res[i,j] = f_xx[i,j] + f_yy[i,j] + x1 * y1 * f_xy[i,j] +
            x1 * f_x[i,j] + y1 * f_y[i,j] +
            (x1*x1 + y1*y1) * f[i,j] + f_x[i,j]^2 + f_y[i,j]^2 +
            x1 * y1 * f_x[i,j] * f_y[i,j] - source(x1,y1)
    end

    res
end


function deriv_operators(hx, hy, Nx::Int, Ny::Int, ord::Int)
    Dx_op  = CenteredDiff{1}(1, ord, hx, Nx)
    Dxx_op = CenteredDiff{1}(2, ord, hx, Nx)

    Dy_op  = CenteredDiff{2}(1, ord, hy, Ny)
    Dyy_op = CenteredDiff{2}(2, ord, hy, Ny)

    Dxx_op, Dyy_op, Dx_op, Dy_op
end

#=
use the Kronecker product (kron) to build the 2-dimensional derivation matrices
from the 1-dimensional ones. see for instance:

  https://en.wikipedia.org/wiki/Kronecker_product

  https://arxiv.org/pdf/1801.01483.pdf (section 5)
=#
function deriv_matrices(Dxx_op::Jecco.FiniteDiffDeriv, Dyy_op::Jecco.FiniteDiffDeriv,
                        Dx_op::Jecco.FiniteDiffDeriv, Dy_op::Jecco.FiniteDiffDeriv)
    Dx  = kron(I(Ny), SparseMatrixCSC(Dx_op))
    Dxx = kron(I(Ny), SparseMatrixCSC(Dxx_op))
    Dy  = kron(SparseMatrixCSC(Dy_op), I(Nx))
    Dyy = kron(SparseMatrixCSC(Dyy_op), I(Nx))

    Dxx, Dyy, Dx * Dy, Dx, Dy
end


# returns axx Dxx + ayy Dyy + axy Dxy + bx Dx + by Dy + cc. note that this
# function overwrites the input matrices to save memory
function build_jacobian(Dxx::SparseMatrixCSC, Dyy::SparseMatrixCSC, Dxy::SparseMatrixCSC,
                        Dx::SparseMatrixCSC, Dy::SparseMatrixCSC,
                        axx::Vector, ayy::Vector, axy::Vector,
                        bx::Vector, by::Vector, cc::Vector)
    Jecco.mul_col!(axx, Dxx)
    Jecco.mul_col!(ayy, Dyy)
    Jecco.mul_col!(axy, Dxy)
    Jecco.mul_col!(bx,  Dx)
    Jecco.mul_col!(by,  Dy)
    ccId = Diagonal(cc)

    Dxx + Dyy + Dxy + Dx + Dy + ccId
end



x_min    = -5.0
x_max    =  5.0
x_nodes  =  128
# x_nodes  =  1024
y_min    = -5.0
y_max    =  5.0
y_nodes  =  64
# y_nodes  =  1024

ord = 4
# ord = 2

xcoord  = Cartesian{1}("x", x_min, x_max, x_nodes, endpoint=false)
ycoord  = Cartesian{2}("y", y_min, y_max, y_nodes, endpoint=false)

hx = Jecco.delta(xcoord)
hy = Jecco.delta(ycoord)

Nx = xcoord.nodes
Ny = ycoord.nodes

f_exact = [exp(-xcoord[i]^2 - ycoord[j]^2) for i in 1:Nx, j in 1:Ny]


Dxx_op, Dyy_op, Dx_op, Dy_op = deriv_operators(hx, hy, Nx, Ny, ord)
Dxx_, Dyy_, Dxy_, Dx_, Dy_   = deriv_matrices(Dxx_op, Dyy_op, Dx_op, Dy_op)

# initial guess
# fsol   = zeros(Nx,Ny)
fsol   = copy(f_exact)

res    = zeros(Nx,Ny)
ind2D  = LinearIndices(fsol)

M = Nx * Ny

b_vec   = zeros(M)
axx     = ones(M)
ayy     = ones(M)
axy     = zeros(M)
bx      = zeros(M)
by      = zeros(M)
cc      = zeros(M)
f0      = zeros(M)


@inbounds for idx in eachindex(f0)
    f0[idx] = fsol[idx]
end


# start looping here

compute_residual!(res, fsol, xcoord[:], ycoord[:], Dxx_op, Dyy_op, Dx_op, Dy_op)


# FIXME

# build (linearized) operator A = Dxx + Dyy + x y Dxy + f_y(0) Dx + f_x(0) Dy + (x^2 + y^2)


f0_x = Dx * f0
f0_y = Dy * f0

# FIXME
for j in 1:Ny, i in 1:Nx
    idx = ind2D[i,j]

    xi  = xcoord[i]
    yi  = ycoord[j]

    axy[idx] = xi * yi
    bx[idx]  = f0_x[idx]
    by[idx]  = f0_y[idx]
    cc[idx]  = xi^2 + yi^2

    b_vec[idx] = source(xi, yi)
end


Dxx = copy(Dxx_)
Dyy = copy(Dyy_)
Dxy = copy(Dxy_)
Dx  = copy(Dx_)
Dy  = copy(Dy_)

A_mat = build_jacobian(Dxx, Dyy, Dxy, Dx, Dy, axx, ayy, axy, bx, by, cc)

# since we're using periodic boundary conditions, the operator A_mat (just like
# the Dx and Dxx operators) is strictly speaking not invertible (it has zero
# determinant) since the solution is not unique. indeed, its LU decomposition
# shouldn't even be defined. for some reason, however, the call to "lu" does in
# fact factorize the matrix. in any case, to be safer, let's instead call
# "factorize", which uses fancy algorithms to determine which is the best way to
# factorize (and which performs a QR decomposition if the LU fails). the inverse
# that is performed probably returns the minimum norm least squares solution, or
# something similar. in any case, for our purposes here we mostly care about
# getting a solution (not necessarily the minimum norm least squares one).
A_fact = factorize(A_mat)
ldiv!(f0, A_fact, b_vec)

@inbounds for idx in eachindex(fsol)
    fsol[idx] = f0[idx]
end

# finish looping


j_slice = div(Ny,2) + 1
x = xcoord[:]

plot(x, f_exact[:,j_slice])
scatter!(x, fsol[:,j_slice])

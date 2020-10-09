
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
                        axx::Diagonal, ayy::Diagonal, axy::Diagonal,
                        bx::Diagonal, by::Diagonal, cc::Diagonal)
    Dxx_ = similar(Dxx)
    Dyy_ = similar(Dyy)
    Dxy_ = similar(Dxy)
    Dx_  = similar(Dx)
    Dy_  = similar(Dy)

    mul!(Dxx_, axx, Dxx)
    mul!(Dyy_, ayy, Dyy)
    mul!(Dxy_, axy, Dxy)
    mul!(Dx_, bx, Dx)
    mul!(Dy_, by, Dy)

    Dxx_ + Dyy_ + Dxy_ + Dx_ + Dy_ + cc
end



x_min    = -5.0
x_max    =  5.0
x_nodes  =  128
# x_nodes  =  256
y_min    = -5.0
y_max    =  5.0
y_nodes  =  64
# y_nodes  =  256

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
Dxx, Dyy, Dxy, Dx, Dy = deriv_matrices(Dxx_op, Dyy_op, Dx_op, Dy_op)

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

itmax   = 8
epsilon = 1e-12

for it in 1:itmax

    compute_residual!(res, fsol, xcoord[:], ycoord[:], Dxx_op, Dyy_op, Dx_op, Dy_op)
    max_res = maximum(abs.(res))

    @show it, max_res
    if max_res < epsilon
        break
    end

    # build (linearized) operator
    # A = Dxx + Dyy + x y Dxy + (x + 2 f0_x + x y f0_y) Dx + (y + 2 f0_y + x y f0_x) Dy + (x^2 + y^2)

    @inbounds for idx in eachindex(f0)
        f0[idx] = fsol[idx]
    end
    f0_x = Dx * f0
    f0_y = Dy * f0

    for j in 1:Ny, i in 1:Nx
        idx = ind2D[i,j]

        x1  = xcoord[i]
        y1  = ycoord[j]

        axy[idx] = x1 * y1
        bx[idx]  = x1 + 2 * f0_x[idx] + x1 * y1 * f0_y[idx]
        by[idx]  = y1 + 2 * f0_y[idx] + x1 * y1 * f0_x[idx]
        cc[idx]  = x1^2 + y1^2

        b_vec[idx] = -res[idx]
    end

    A_mat = build_jacobian(Dxx, Dyy, Dxy, Dx, Dy,
                           Diagonal(axx), Diagonal(ayy), Diagonal(axy),
                           Diagonal(bx), Diagonal(by), Diagonal(cc))

    A_fact = factorize(A_mat)
    ldiv!(f0, A_fact, b_vec)

    # update solution
    @inbounds for idx in eachindex(fsol)
        fsol[idx] += f0[idx]
    end

end


j_slice = div(Ny,2) + 1
x = xcoord[:]

plot(x, f_exact[:,j_slice])
scatter!(x, fsol[:,j_slice])

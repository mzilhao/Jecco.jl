
using Jecco
using LinearAlgebra
using SparseArrays

source(x) = exp(-x^2) * (-2 + 3*x^2)

x_min   = -5.0
x_max   = 5.0
x_nodes = 256

xcoord  = CartCoord{1}("x", x_min, x_max, x_nodes, endpoint=false)

hx = Jecco.delta(xcoord)
Nx = xcoord.nodes
x  = xcoord[:]

f_exact = exp.(-x.^2)

ord = 4

Dx_op  = CenteredDiff{1}(1, ord, hx, Nx)
Dxx_op = CenteredDiff{1}(2, ord, hx, Nx)

Dx  = SparseMatrixCSC(Dx_op)
Dxx = SparseMatrixCSC(Dxx_op)

# TODO: building the operator this way takes longer than performing the LU
# factorization below. see if there are better ways to do this.
A_mat = Dxx + x .* Dx + x.^2 .* sparse(I, Nx, Nx)
b_vec = source.(x)

A_fact = lu(A_mat)
sol    = zero(b_vec)
ldiv!(sol, A_fact, b_vec)

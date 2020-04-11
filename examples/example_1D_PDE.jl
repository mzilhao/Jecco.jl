
using Jecco
using LinearAlgebra
using SparseArrays
using Plots

# return M_ij = x_i A_ij  (no sum in i)
function mul_col(x::Vector, A::SparseMatrixCSC)
    @assert size(x)[1] == size(A)[1]

    ii, jj, V = findnz(A)
    @inbounds for idx in eachindex(V)
        V[idx] *= x[ii[idx]]
    end
    sparse(ii,jj,V)
end

function build_operator(aa::Vector, Dxx::SparseMatrixCSC,
                        bb::Vector, Dx::SparseMatrixCSC, cc::Vector)
    aaDxx = mul_col(aa, Dxx)
    bbDx  = mul_col(bb, Dx)
    ccId  = Diagonal(cc)

    aaDxx + bbDx + ccId
end

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

aa = ones(Nx)
bb = copy(x)
cc = x.^2

# build operator A = Dxx + x Dx + x^2
A_mat = build_operator(aa, Dxx, bb, Dx, cc)
b_vec = source.(x)

A_fact = factorize(A_mat)
sol    = A_fact \ b_vec

plot(x, f_exact)
scatter!(x, sol)

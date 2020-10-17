
using Jecco
using LinearAlgebra
using SparseArrays
using Plots

# returns A_ij = x_i A_ij (no sum in i). note that A itself is also changed.
# this is equivalent to the operation A = x .* A, but it's much more efficient
function mul_col!(x::Vector, A::SparseMatrixCSC)
    @assert size(x)[1] == size(A)[1]

    rows = rowvals(A)
    vals = nonzeros(A)
    @inbounds for idx in eachindex(vals)
        row = rows[idx]
        vals[idx] *= x[row]
    end
    A
end

#=
build operator

   aa Dxx + bb Dx + cc

by overwriting Dxx with the returned value. the construction here is more
efficient than directly summing the sparse arrays. note that this construction
only works because we know that the returned array has the *same* sparsity
pattern as Dxx! we couldn't, for instance, mutate Dx instead since the Dx array
has zeros along the diagonal (which would be skipped in the construction below).
=#
function build_operator!(Dxx::SparseMatrixCSC, Dx::SparseMatrixCSC,
                         aa::Vector, bb::Vector, cc::Vector)
    mul_col!(aa, Dxx)
    mul_col!(bb, Dx)
    ccId = Diagonal(cc)

    rows = rowvals(Dxx)
    vals = nonzeros(Dxx)
    m, n = size(Dxx)

    # cf. https://docs.julialang.org/en/v1/stdlib/SparseArrays/#SparseArrays.nzrange
    @inbounds for j = 1:n
        @inbounds for i in nzrange(Dxx, j)
            row = rows[i]
            vals[i] += Dx[row,j] + ccId[row,j]
        end
    end
    Dxx
end

source(x) = exp(-x^2) * (-2 + 3*x^2)


x_min   = -5.0
x_max   = 5.0
x_nodes = 256

xcoord  = Cartesian{1}("x", x_min, x_max, x_nodes, endpoint=false)

hx = Jecco.delta(xcoord)
Nx = xcoord.nodes
x  = xcoord[:]

f_exact = exp.(-x.^2)

ord = 4

Dx_op  = CenteredDiff{1}(1, ord, hx, Nx)
Dxx_op = CenteredDiff{1}(2, ord, hx, Nx)
Dx     = SparseMatrixCSC(Dx_op)
Dxx    = SparseMatrixCSC(Dxx_op)

aa = ones(Nx)
bb = copy(x)
cc = x.^2

A_mat  = copy(Dxx)
Dx_tmp = copy(Dx)

# build operator A = Dxx + x Dx + x^2
build_operator!(A_mat, Dx_tmp, aa, bb, cc)

b_vec = source.(x)

A_fact = factorize(A_mat)
sol    = A_fact \ b_vec

plot(x, f_exact)
scatter!(x, sol)

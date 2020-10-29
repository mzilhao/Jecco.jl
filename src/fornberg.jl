
#= Fornberg algorithm

Implements the Fornberg (1988) algorithm
(https://doi.org/10.1090/S0025-5718-1988-0935077-0), adapted from the
DiffEqOperators module.
=#
"""
    calculate_weights(m, n, s)

Return finite difference weights. Notation follows that of Fornberg, used in
Mathematica as well:
https://reference.wolfram.com/language/tutorial/NDSolveMethodOfLines.html

# Arguments
* `m::Int`: order of the derivative
* `n::Int`: number of grid intervals enclosed in the stencil
* `s::Int`: number of grid intervals between the point at which the derivative
            is approximated and the leftmost edge of the stencil. centred
            formulas always have s=n/2.
"""
function calculate_weights(m::Int, n::Int, s::Int)
    stencil_length = n + 1

    # dummy array with the relative coordinates of the stencil (eg. central
    # differences need coordinates centred at 0)
    x  = -s : (stencil_length - s - 1)
    # point in the array 'x' for which we need the coefficients
    x0 = x[1] + s

    N = length(x)
    @assert m < N "Not enough points for the requested order."

    c1 = 1
    c4 = x[1] - x0
    C = zeros(Rational{Int}, N, m+1)
    C[1,1] = 1
    @inbounds for i in 1:N-1
        i1 = i + 1
        mn = Base.min(i, m)
        c2 = 1
        c5 = c4
        c4 = x[i1] - x0
        for j in 0:i-1
            j1 = j + 1
            c3 = x[i1] - x[j1]
            c2 *= c3
            if j == i-1
                for s in mn:-1:1
                    s1 = s + 1
                    C[i1,s1] = c1*(s*C[i,s] - c5*C[i,s1]) // c2
                end
                C[i1,1] = -c1*c5*C[i,1] // c2
           end
            for s in mn:-1:1
                s1 = s + 1
                C[j1,s1] = (c4*C[j1,s1] - s*C[j1,s]) // c3
            end
            C[j1,1] = c4 * C[j1,1] // c3
        end
        c1 = c2
    end

    #=
    Fix the problem of numerical instability which occurs when the sum of the
    stencil_coefficients is not exactly 0:
https://scicomp.stackexchange.com/questions/11249/numerical-derivative-and-finite-difference-coefficients-any-update-of-the-fornb

    Stack Overflow answer on this issue:
    http://epubs.siam.org/doi/pdf/10.1137/S0036144596322507 - Modified Fornberg Algorithm
    =#

    _C = C[:,end]
    _C[div(N,2)+1] -= sum(_C)
    return _C
end

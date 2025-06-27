
import Statistics: mean

function solve(integrator::ODEIntegrator, tmax::Number)
    while true
        t = integrator.t
        if t >= tmax
            break
        end
        step!(integrator)
    end

    tf  = integrator.t
    sol = integrator.u

    tf, sol
end

function solve(prob::ODEProblem, alg::ODEAlgorithm, dt::Number)
    integrator = ODEIntegrator(prob, alg, dt)
    solve(integrator, prob.tspan[2])
end

function convergence_order(prob::ODEProblem, alg::ODEAlgorithm, dts::AbstractArray,
                          sol_expected::AbstractArray)
    _sols = [solve(prob, alg, dt) for dt in dts]

    tfs  = [sol[1] for sol in _sols]
    sols = [sol[2] for sol in _sols]

    tmax = prob.tspan[2]
    @assert all(tfs .≈ tmax)

    # l2-norm of the difference to the analytical solution. we take the l2 norm to
    # reduce array to a number
    diffs_l2 = [sqrt(sum((sol .- sol_expected).^2)) for sol in sols]

    N = length(diffs_l2)
    S = Vector{eltype(diffs_l2)}(undef, N-1)
    for i in 1:(N-1)
        S[i] = log2(diffs_l2[i+1] / diffs_l2[i])
    end

    conv_ord = mean(S)
end

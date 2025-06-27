using Random
Random.seed!(42)

const ODE = Jecco.ODESolver

function f_analytic(u0, p, t)
    @. u0 * exp(p * t)
end

function rhs!(du,u,p,t)
    @. du = p * u
end


@testset "convergence tests:" begin
    u0 = rand(4,2)
    p  = 1.01
    tspan = (0.0, 1.0)

    prob = ODE.ODEProblem(rhs!, u0, tspan, p)

    tmax = tspan[2]
    sol_analytic = f_analytic(u0, p, tmax)

    # ensure that there is a factor of 2 between each dt so that in the end we can
    # simply take a log2 to check for convergence
    dts = 1 ./ 2 .^ (8:-1:4)

    ord = ODE.convergence_order(prob, ODE.RK4(), dts, sol_analytic)
    @test ord ≈ 4 atol=0.2

    ord = ODE.convergence_order(prob, ODE.AB3(), dts, sol_analytic)
    @test ord ≈ 3 atol=0.2
end

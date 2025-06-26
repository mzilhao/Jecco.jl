module ODESolver

import Base.Threads.@threads

import OrdinaryDiffEq: step!
export step!

abstract type ODEAlgorithm end

struct ODEProblem{F,uType,tType,P}
    "Assuming an ODE of the form `du/dt = f(u,p,t)` where `p` are constant parameters,
the function `f` should have the form `f(du, u, p, t)` where `du` is the RHS vector."
    f  :: F
    "Initial conditions (state vector at t=t0)"
    u0 :: uType
    "Time span of the problem"
    tspan :: tType
    "Parameters to be passed as the third argument of `f`."
    p  :: P

    function ODEProblem(f, u0, tspan, p)
        @assert length(tspan) == 2
        new{typeof(f), typeof(u0), typeof(tspan), typeof(p)}(f, u0, tspan, p)
    end
end

ODEProblem(f::Function, u0::Number, t, p) = ODEProblem(f, [u0], t, p)

ODEProblem(f::Function, u0, t) = ODEProblem(f, u0, t, nothing)

include("algorithms.jl")
include("integrators.jl")
include("utils.jl")

end

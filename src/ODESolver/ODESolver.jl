module ODESolver

import Base.Threads.@threads

import OrdinaryDiffEq: step!
export step!

abstract type ODEAlgorithm end
abstract type Tableau end

struct ODEProblem{F,uType,tType,P}
    "Assuming an ODE of the form `du/dt = f(u,p,t)` where `p` are constant parameters,
the function `f` should have the form `f(du, u, p, t)` where `du` is the RHS vector."
    f  :: F
    "Initial conditions (state vector at t=t0)"
    u0 :: uType
    "Initial time"
    t0 :: tType
    "Parameters to be passed as the third argument of `f`."
    p  :: P
end

ODEProblem(f::Function, u0::Number, t0, p) = ODEProblem(f, [u0], t0, p)

ODEProblem(f::Function, u0, t0) = ODEProblem(f, u0, t0, nothing)

include("tableaus.jl")
include("algorithms.jl")
include("integrators.jl")

end

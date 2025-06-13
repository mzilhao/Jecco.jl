mutable struct ODEIntegrator{F, P, uType, tType, algType, cacheType}
    f     :: F
    p     :: P
    u     :: uType
    t     :: tType
    uprev :: uType
    tprev :: tType
    dt    :: tType
    alg   :: algType
    cache :: cacheType
end


"""
    ODEIntegrator(prob, alg, dt)

# Arguments
* `prob::ODEProblem`
* `alg::ODEAlgorithm`
* `dt::AbstractFloat`
"""
function ODEIntegrator(prob::ODEProblem{F,uType,tType,P}, alg::algType,
                       dt::tType) where {F,uType,tType,P,algType<:ODEAlgorithm}
    f = prob.f
    p = prob.p
    u = deepcopy(prob.u0)
    t = prob.t0

    uprev = deepcopy(u)
    tprev = t

    cache = alg_cache(alg, u)
    ODEIntegrator{F,P,uType,tType,algType,
                  typeof(cache)}(f, p, u, t, uprev, tprev, dt, alg, cache)
end

@inline step!(integrator::ODEIntegrator) = step!(integrator, integrator.alg)


function step!(integrator::ODEIntegrator, alg::RK2)
    rhs!   = integrator.f
    p      = integrator.p
    u      = integrator.u
    t      = integrator.t
    uprev  = integrator.uprev
    dt     = integrator.dt
    cache  = integrator.cache

    tab = alg.tableau
    a21 = tab.a21
    b1  = tab.b1
    b2  = tab.b2
    c2  = tab.c2

    k1  = cache.k1
    k2  = cache.k2
    tmp = cache.tmp

    # store previous timelevel
    @inbounds @threads for i in eachindex(uprev)
        uprev[i] = u[i]
    end

    rhs!(k1, uprev, p, t)

    @inbounds @threads for i in eachindex(tmp)
        tmp[i] = u[i] + a21 * k1[i] * dt
    end

    rhs!(k2, tmp, p, t + c2 * dt)

    # update u
    @inbounds @threads for i in eachindex(u)
        u[i] += (b1 * k1[i] + b2 * k2[i]) * dt
    end

    # update t
    integrator.tprev = t
    integrator.t += dt

    nothing
end

function step!(integrator::ODEIntegrator, alg::RK4)
    rhs!   = integrator.f
    p      = integrator.p
    u      = integrator.u
    t      = integrator.t
    uprev  = integrator.uprev
    dt     = integrator.dt
    cache  = integrator.cache

    tab = alg.tableau
    a21 = tab.a21
    a31 = tab.a31
    a32 = tab.a32
    a41 = tab.a41
    a42 = tab.a42
    a43 = tab.a43
    b1  = tab.b1
    b2  = tab.b2
    b3  = tab.b3
    b4  = tab.b4
    c2  = tab.c2
    c3  = tab.c3
    c4  = tab.c4

    k1  = cache.k1
    k2  = cache.k2
    k3  = cache.k3
    k4  = cache.k4
    tmp = cache.tmp

    # store previous timelevel
    @inbounds @threads for i in eachindex(uprev)
        uprev[i] = u[i]
    end

    rhs!(k1, uprev, p, t)
    @inbounds @threads for i in eachindex(tmp)
        tmp[i] = u[i] + a21 * k1[i] * dt
    end
    rhs!(k2, tmp, p, t + c2 * dt)

    @inbounds @threads for i in eachindex(tmp)
        tmp[i] = u[i] + (a31 * k1[i] + a32 * k2[i]) * dt
    end
    rhs!(k3, tmp, p,  t + c3 * dt)

    @inbounds @threads for i in eachindex(tmp)
        tmp[i] = u[i] + (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]) * dt
    end
    rhs!(k4, tmp, p, t + c4 * dt)

    # update u
    @inbounds @threads for i in eachindex(u)
        u[i] += (b1 * k1[i] + b2 * k2[i] + b3 * k3[i] + b4 * k4[i]) * dt
    end

    # update t
    integrator.tprev = t
    integrator.t += dt

    nothing
end

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

    k1  = cache.k1
    k2  = cache.k2
    tmp = cache.tmp

    # store previous timelevel
    @inbounds @threads for i in eachindex(uprev)
        uprev[i] = u[i]
    end

    rhs!(k1, uprev, p, t)

    @inbounds @threads for i in eachindex(tmp)
        tmp[i] = u[i] + k1[i] * (dt/2)
    end

    rhs!(k2, tmp, p, t + dt/2)

    # update u
    @inbounds @threads for i in eachindex(u)
        u[i] += k2[i] * dt
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
        tmp[i] = u[i] + k1[i] * (dt/2)
    end
    rhs!(k2, tmp, p, t + dt/2)

    @inbounds @threads for i in eachindex(tmp)
        tmp[i] = u[i] + k2[i] * (dt/2)
    end
    rhs!(k3, tmp, p,  t + dt/2)

    @inbounds @threads for i in eachindex(tmp)
        tmp[i] = u[i] + k3[i] * dt
    end
    rhs!(k4, tmp, p, t + dt)

    # update u
    @inbounds @threads for i in eachindex(u)
        u[i] += (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) * (dt/6)
    end

    # update t
    integrator.tprev = t
    integrator.t += dt

    nothing
end

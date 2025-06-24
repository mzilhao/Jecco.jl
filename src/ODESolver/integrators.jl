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

@inline step!(integrator::ODEIntegrator) = step!(integrator, integrator.cache)


function step!(integrator::ODEIntegrator, cache::RK2Cache)
    rhs!   = integrator.f
    p      = integrator.p
    u      = integrator.u
    t      = integrator.t
    uprev  = integrator.uprev
    dt     = integrator.dt

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

function step!(integrator::ODEIntegrator, cache::RK4Cache)
    rhs!   = integrator.f
    p      = integrator.p
    u      = integrator.u
    t      = integrator.t
    uprev  = integrator.uprev
    dt     = integrator.dt

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

function step!(integrator::ODEIntegrator, cache::AB3Cache)
    rhs! = integrator.f
    p    = integrator.p

    k1  = cache.k1
    k2  = cache.k2
    k3  = cache.k3

    # when calling this stepping function for the first time, we don't yet have
    # the values for the previous time-levels to perform the update. therefore,
    # we step backwards in time using RK4
    if cache.init
        integrator.dt *= -1
        u0 = deepcopy(integrator.u)
        t0 = integrator.t
        rk4_cache = alg_cache(RK4(), u0)

        rhs!(k1, integrator.u, p, integrator.t)
        step!(integrator, rk4_cache)
        rhs!(k2, integrator.u, p, integrator.t)
        step!(integrator, rk4_cache)
        rhs!(k3, integrator.u, p, integrator.t)

        # now reset starting values
        integrator.u  .= u0
        integrator.t   = t0
        integrator.dt *= -1

        cache.init = false
        rk4_cache  = nothing
        u0         = nothing
    end

    u      = integrator.u
    uprev  = integrator.uprev
    dt     = integrator.dt

    # rotate timelevels
    u, uprev = uprev, u

    # update u. we need to have the multithreaded kernel in a separate function
    # to avoid type-instabilities. see discussion here:
    # https://discourse.julialang.org/t/mysterious-type-instability-performance-hit-with-simple-threads/117404
    _AB3_update_u!(u, uprev, k1, k2, k3, dt)

    # rotate timelevels
    cache.k3, cache.k2, cache.k1 = k2, k1, k3
    rhs!(cache.k1, u, p, integrator.t + dt)

    integrator.u = u
    integrator.uprev = uprev

    # update t
    integrator.tprev = integrator.t
    integrator.t += dt

    nothing
end

function _AB3_update_u!(u,uprev,k1,k2,k3,dt)
    @inbounds @threads for i in eachindex(u)
        u[i] = uprev[i] + (dt / 12) * (23 * k1[i] - 16 * k2[i] + 5 * k3[i])
    end
    nothing
end

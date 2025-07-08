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
function ODEIntegrator(prob::ODEProblem, alg::ODEAlgorithm, dt::AbstractFloat)
    f = prob.f
    p = prob.p
    u = deepcopy(prob.u0)
    t = prob.tspan[1]

    uprev = deepcopy(u)
    tprev = t

    cache = alg_cache(alg, u)
    ODEIntegrator(f, p, u, t, uprev, tprev, dt, alg, cache)
end

@inline step!(integrator::ODEIntegrator) = step!(integrator, integrator.cache)


# RK2

function step!(integrator::ODEIntegrator, cache::RK2Cache)
    rhs!   = integrator.f
    p      = integrator.p
    u      = integrator.u
    t      = integrator.t
    uprev  = integrator.uprev
    dt     = integrator.dt

    tmp = cache.tmp

    # rotate timelevels
    u, uprev = uprev, u

    # k1
    rhs!(tmp, uprev, p, t)
    _RK2_update_k1!(u, uprev, tmp, dt)

    # k2
    rhs!(tmp, u, p, t + dt/2)

    # update u
    _RK2_update_u!(u, uprev, tmp, dt)

    integrator.u = u
    integrator.uprev = uprev

    # update t
    integrator.tprev = t
    integrator.t += dt

    nothing
end
function _RK2_update_k1!(u, uprev, rhs, dt)
    @inbounds @threads for i in eachindex(u)
        u[i] = uprev[i] + rhs[i] * (dt/2)
    end
    nothing
end
function _RK2_update_u!(u, uprev, rhs, dt)
    @inbounds @threads for i in eachindex(u)
        u[i] = uprev[i] + rhs[i] * dt
    end
    nothing
end


# RK 4

function step!(integrator::ODEIntegrator, cache::RK4Cache)
    rhs!   = integrator.f
    p      = integrator.p
    u      = integrator.u
    t      = integrator.t
    uprev  = integrator.uprev
    dt     = integrator.dt

    k   = cache.k
    tmp = cache.tmp

    # rotate timelevels
    u, uprev = uprev, u

    # k1
    rhs!(tmp, uprev, p, t)
    _RK4_update_k1!(u,k,uprev,tmp,dt)

    # k2
    rhs!(tmp, u, p, t + dt/2)
    _RK4_update_k2!(u,k,uprev,tmp,dt)

    # k3
    rhs!(tmp, u, p, t + dt/2)
    _RK4_update_k3!(u,k,uprev,tmp,dt)

    # k4, final update
    rhs!(tmp, u, p, t + dt)
    _RK4_update_k4!(u,k,uprev,tmp,dt)

    integrator.u = u
    integrator.uprev = uprev

    # update t
    integrator.tprev = t
    integrator.t += dt

    nothing
end
function _RK4_update_k1!(u, k, uprev, rhs, dt)
    @inbounds @threads for i in eachindex(u)
        u[i] = uprev[i] + rhs[i] * (dt/2)
        k[i] = rhs[i]
    end
    nothing
end
function _RK4_update_k2!(u, k, uprev, rhs, dt)
    @inbounds @threads for i in eachindex(u)
        u[i]  = uprev[i] + rhs[i] * (dt/2)
        k[i] += 2 * rhs[i]
    end
    nothing
end
function _RK4_update_k3!(u, k, uprev, rhs, dt)
    @inbounds @threads for i in eachindex(u)
        u[i]  = uprev[i] + rhs[i] * dt
        k[i] += 2 * rhs[i]
    end
    nothing
end
function _RK4_update_k4!(u, k, uprev, rhs, dt)
    @inbounds @threads for i in eachindex(u)
        k[i] += rhs[i]
        u[i]  = uprev[i] + k[i] * (dt/6)
    end
    nothing
end


# AB3

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

"""
    AlgorithmCache

Abstract type for the cache of each algorithm. Subtypes of this Abstract type
need to implement the following method:

```alg_cache(alg::ODEAlgorithm, u)```

"""
abstract type AlgorithmCache end

function alg_cache(alg::ODEAlgorithm, u) end

@inline alg_cache(alg::ODEAlgorithm, u::Number) = alg_cache(alg, [u])


struct RK2 <: ODEAlgorithm end

mutable struct RK2Cache{A} <: AlgorithmCache
    k1  :: A
    k2  :: A
    tmp :: A
end

function alg_cache(alg::RK2, u::AbstractArray)
    k1  = zero(u)
    k2  = zero(u)
    tmp = zero(u)
    RK2Cache(k1, k2, tmp)
end


struct RK4 <: ODEAlgorithm end

mutable struct RK4Cache{A} <: AlgorithmCache
    k1  :: A
    k2  :: A
    k3  :: A
    k4  :: A
    tmp :: A
end

function alg_cache(alg::RK4, u::AbstractArray)
    k1  = zero(u)
    k2  = zero(u)
    k3  = zero(u)
    k4  = zero(u)
    tmp = zero(u)
    RK4Cache(k1, k2, k3, k4, tmp)
end

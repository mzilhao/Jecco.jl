
"""
    AbstractPartition{N,A}

Container to store a set of objects across different partitions. Needs a field
`x`.
"""
abstract type AbstractPartition{N,A} end

@inline Base.iterate(ff::AbstractPartition)         = iterate(ff.x)
@inline Base.iterate(ff::AbstractPartition, i::Int) = iterate(ff.x, i)

@inline Base.length(ff::AbstractPartition{N}) where{N} = N
@inline Base.firstindex(ff::AbstractPartition) = 1
@inline Base.lastindex(ff::AbstractPartition)  = length(ff)

@inline Base.getindex(ff::AbstractPartition, i::Int) = ff.x[i]


"""
    FlattenedVector{T,N,A} <: AbstractVector{T}

Container to store a set of `Arrays` as elements of an `NTuple`, or similar. The
idea is to treat them as a single column vector for all intents and purposes.
Inspired in `ArrayPartition` from `RecursiveArrayTools`. Needs a field `x`.
"""
abstract type FlattenedVector{T,N,A} <: AbstractVector{T} end

@inline Base.length(ff::FlattenedVector) = sum((length(x) for x in ff.x))
@inline Base.size(ff::FlattenedVector)   = (length(ff),)

# indexing. this is just a linear indexing through all the arrays. adapted from
# RecursiveArrayTools
@inline Base.firstindex(ff::FlattenedVector) = 1
@inline Base.lastindex(ff::FlattenedVector)  = length(ff)

@inline function Base.getindex(ff::FlattenedVector, i::Int)
    @inbounds for j in 1:length(ff.x)
        f  = ff.x[j]
        i -= length(f)
        if i <= 0
            return f[length(f)+i]
        end
    end
end
@inline function Base.setindex!(ff::FlattenedVector, v, i::Int)
    @inbounds for j in 1:length(ff.x)
        f  = ff.x[j]
        i -= length(f)
        if i <= 0
            f[length(f)+i] = v
            break
        end
    end
end

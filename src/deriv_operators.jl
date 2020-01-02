
import Base: *

abstract type AbstractDerivOperator{T,N} end

D1_42_weights() = [1.0, -8.0, 0.0, 8.0, -1.0] ./ 12.0

D2_42_weights() = [-1.0, 16.0, -30.0, 16.0, -1.0] ./ 12.0


struct FiniteDiffDeriv{T<:Real,N,T2,S} <: AbstractDerivOperator{T,N}
    derivative_order        :: Int
    approximation_order     :: Int
    dx                      :: T2
    len                     :: Int
    stencil_length          :: Int
    stencil_coefs           :: S
end

struct SpectralDeriv{T<:Real,N,S} <: AbstractDerivOperator{T,N}
    derivative_order        :: Int
    len                     :: Int
    D                       :: S
end

struct CenteredDiff{N} end

function CenteredDiff{N}(derivative_order::Int,
                         approximation_order::Int, dx::T,
                         len::Int) where {T<:Real,N}
    @assert approximation_order > 1 "approximation_order must be greater than 1."

    stencil_length = derivative_order + approximation_order - 1 +
        (derivative_order+approximation_order)%2


    # TODO: this can all be improved...

    if approximation_order != 4
        error("approximation_order not implemented yet")
    end

    if derivative_order == 1
        weights = D1_42_weights()
    elseif derivative_order == 2
        weights = D2_42_weights()
    else
        error("derivative_order not implemented yet")
    end

    stencil_coefs = (1/dx^derivative_order) .* weights

    FiniteDiffDeriv{T,N,T,typeof(stencil_coefs)}(derivative_order, approximation_order,
                                                 dx, len, stencil_length, stencil_coefs)
end

CenteredDiff(args...) = CenteredDiff{1}(args...)


struct ChebDeriv{N} end

# TODO
function ChebDeriv{N}(derivative_order::Int, delta::T, len::Int) where {T<:Real,N}

end

ChebDeriv(args...) = ChebDeriv{1}(args...)


# mul! done by convolution
function LinearAlgebra.mul!(x_temp::AbstractVector, A::FiniteDiffDeriv, x::AbstractVector)
    convolve!(x_temp, x, A)
end


function LinearAlgebra.mul!(fout::AbstractArray{T}, A::FiniteDiffDeriv{T,N},
                            f::AbstractArray{T}) where {T,N}
    # dimension of f
    ndim   = ndims(f)
    dims   = [axes(f)...]

    otherdims = setdiff(1:ndim, N)

    # an array of type "Any" with the first index of each dimension
    idx = Any[first(ind) for ind in axes(f)]
    # and now set the Nth dimension entry to a Colon, ie, if ndim = 4 and N = 2,
    # idx = [1,:,1,1]
    setindex!(idx, :, N)

    itershape = tuple(dims[otherdims]...)
    # this generates an iterator along all the otherdims, ie, without touching the
    # N-dimension
    indices = Iterators.drop(CartesianIndices(itershape), 0)

    nidx = length(otherdims)
    # index I will loop along all indices in otherdims, without touching the
    # N-dimension, and idx is updated through the call to replace_tuples!
    @fastmath @inbounds for I in indices
        Base.replace_tuples!(nidx, idx, idx, otherdims, I)
        mul!(view(fout, idx...), A, view(f, idx...))
    end

    nothing
end

# convolution operation to act with derivatives on vectors. this currently
# assumes periodic BCs. adapted from
# DiffEqOperators.jl/src/derivative_operators/convolutions.jl
function convolve!(xout::AbstractVector{T}, x::AbstractVector{T},
                   A::FiniteDiffDeriv) where {T<:Real}
    N = length(x)
    coeffs = A.stencil_coefs
    mid = div(A.stencil_length, 2) + 1

    @fastmath @inbounds for i in 1:N
        sum_i = zero(T)
        @inbounds for idx in 1:A.stencil_length
            # imposing periodicity
            i_circ = 1 + mod(i - (mid-idx) - 1, N)
            sum_i += coeffs[idx] * x[i_circ]
        end
        xout[i] = sum_i
    end
    nothing
end

function *(A::AbstractDerivOperator, x::AbstractArray)
    y = similar(x)
    LinearAlgebra.mul!(y, A, x)
    y
end

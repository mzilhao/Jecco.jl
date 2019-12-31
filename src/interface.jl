
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

struct SpectralDeriv{T<:Real,N} <: AbstractDerivOperator{T,N}
    derivative_order        :: Int
    len                     :: Int
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

    stencil_coefs = (1/dx^derivative_order) * weights

    FiniteDiffDeriv{T,N,T,typeof(stencil_coefs)}(derivative_order, approximation_order,
                                                 dx, len, stencil_length, stencil_coefs)
end

CenteredDiff(args...) = CenteredDiff{1}(args...)

diff_axis(A::AbstractDerivOperator{T,N}) where {T,N} = N

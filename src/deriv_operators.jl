
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

function LinearAlgebra.mul!(df::AbstractArray{T}, A::FiniteDiffDeriv{T,N},
                            f::AbstractArray{T}) where {T,N}
    Rpre  = CartesianIndices(axes(f)[1:N-1])
    Rpost = CartesianIndices(axes(f)[N+1:end])

    _mul_loop!(df, A, f, Rpre, Rpost)
end

# this works as follows. suppose we are taking the derivative of a 4-dimensional
# array g, with coordinates w,x,y,z, and size
#
#   size(g) = (nw, nx, ny, nz)
#
# let's further suppose that we want the derivative along the x-direction (here,
# the second entry). this fixes
#
#   N = 2
#
# and then
#
#   Rpre  = CartesianIndices((nw,))
#   Rpost = CartesianIndices((ny,nz))
#
# now, the entry [Ipre,:,Ipost]
#
# slices *only* along the x-direction, for fixed Ipre and Ipost. Ipre and Ipost
# loop along Rpre and Rpost (respectively) without touching the x-direction. we
# further take care to loop in memory-order, since julia's arrays are faster in
# the first dimension.
#
# adapted from http://julialang.org/blog/2016/02/iteration
#
@noinline function _mul_loop!(df, A, f, Rpre, Rpost)
    @fastmath @inbounds for Ipost in Rpost
        @inbounds for Ipre in Rpre
            @views mul!(df[Ipre,:,Ipost], A, f[Ipre,:,Ipost])
        end
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

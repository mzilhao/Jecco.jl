
import Base: *, copyto!
import LinearAlgebra: mul!
using SparseArrays

abstract type AbstractDerivOperator{T,N} end

abstract type AbstractFiniteDiff{T,N} <: AbstractDerivOperator{T,N} end

struct PeriodicFD{T<:AbstractFloat,N,S} <: AbstractFiniteDiff{T,N}
    derivative_order        :: Int
    approximation_order     :: Int
    dx                      :: T
    len                     :: Int
    stencil_length          :: Int
    stencil_offset          :: Int
    stencil_coefs           :: S
end

struct FiniteDiffDeriv{T<:AbstractFloat,N,S1,S2,S3} <: AbstractFiniteDiff{T,N}
    derivative_order        :: Int
    approximation_order     :: Int
    dx                      :: T
    len                     :: Int
    stencil_length          :: Int
    stencil_offset          :: Int
    stencil_coefs           :: S1
    boundary_stencil_length :: Int
    # these are assumed to be ordered from the boundary inwards
    low_boundary_offsets    :: Vector{Int}
    low_boundary_coefs      :: S2
    high_boundary_offsets   :: Vector{Int}
    high_boundary_coefs     :: S3
end

struct SpectralDeriv{T<:Real,N,S} <: AbstractDerivOperator{T,N}
    derivative_order        :: Int
    len                     :: Int
    D                       :: S
end

struct CenteredDiff{N} end

function CenteredDiff{N}(derivative_order::Int,
                         approximation_order::Int, dx::T,
                         len::Int) where {T<:AbstractFloat,N}
    @assert approximation_order > 1 "approximation_order must be greater than 1."

    stencil_length = derivative_order + approximation_order - 1 +
        (derivative_order+approximation_order)%2

    stencil_offset = div(stencil_length-1,2)
    weights = calculate_weights(derivative_order, stencil_length-1, stencil_offset)

    stencil_coefs = (1/dx^derivative_order) .* weights

    PeriodicFD{T,N,typeof(stencil_coefs)}(derivative_order, approximation_order,
                                          dx, len, stencil_length,
                                          stencil_offset, stencil_coefs)
end

CenteredDiff(args...) = CenteredDiff{1}(args...)


struct EqualSizeStencilFD{N} end

"""
    EqualSizeStencilFD{N}(derivative_order, approximation_order, dx, len)

Finite difference operator where interior and boundary stencils have the same
length and same approximation order. Low (high) boundary operators will use as
many points to the left (right) as possible. For approximation_order > 2 there
is more than one boundary point, and each boundary point uses a different
stencil. Note that these are *not* truncation matched.
"""
function EqualSizeStencilFD{N}(derivative_order::Int,
                               approximation_order::Int, dx::T,
                               len::Int) where {T<:AbstractFloat,N}
    @assert approximation_order > 1 "approximation_order must be greater than 1."

    stencil_length = derivative_order + approximation_order - 1 +
        (derivative_order+approximation_order)%2
    stencil_offset = div(stencil_length-1,2)
    weights = calculate_weights(derivative_order, stencil_length-1, stencil_offset)

    stencil_coefs = (1/dx^derivative_order) .* weights

    boundary_stencil_length = stencil_length

    low_boundary_size  = div(approximation_order,2)
    high_boundary_size = div(approximation_order,2)

    # remember to order from the boundary inwards, which means the offset array
    # is increasing in the low boundary and decreasing in the high boundary
    low_boundary_offsets  = collect(0:low_boundary_size-1)
    high_boundary_offsets = collect(boundary_stencil_length-1:-1:high_boundary_size+1)

    low_weights = [calculate_weights(derivative_order, stencil_length-1, offset)
                   for offset in low_boundary_offsets]
    low_boundary_coefs = (1/dx^derivative_order) .* low_weights

    high_weights = [calculate_weights(derivative_order, stencil_length-1, offset)
                    for offset in high_boundary_offsets]
    high_boundary_coefs = (1/dx^derivative_order) .* high_weights

    FiniteDiffDeriv{T,N,typeof(stencil_coefs),typeof(low_boundary_coefs),
                    typeof(high_boundary_coefs)}(
                        derivative_order, approximation_order,
                        dx, len, stencil_length,
                        stencil_offset, stencil_coefs,
                        boundary_stencil_length, low_boundary_offsets,
                        low_boundary_coefs, high_boundary_offsets, high_boundary_coefs
                    )
end

EqualSizeStencilFD(args...) = EqualSizeStencilFD{1}(args...)


struct ChebDeriv{N} end

function ChebDeriv{N}(derivative_order::Int, xmin::T, xmax::T, len::Int) where {T<:Real,N}
    x, Dx, Dxx = cheb(xmin, xmax, len)

    if derivative_order == 1
        D = Dx
    elseif derivative_order == 2
        D = Dxx
    else
        error("derivative_order not implemented")
    end

    SpectralDeriv{T,N,typeof(D)}(derivative_order, len, D)
end

ChebDeriv(args...) = ChebDeriv{1}(args...)


struct FourierDeriv{N} end

function FourierDeriv{N}(derivative_order::Int, xsize::T, len::Int) where {T<:Real,N}
    Dx, Dxx = fourier(xsize, len)

    if derivative_order == 1
        D = Dx
    elseif derivative_order == 2
        D = Dxx
    else
        error("derivative_order not implemented")
    end

    SpectralDeriv{T,N,typeof(D)}(derivative_order, len, D)
end

FourierDeriv(args...) = FourierDeriv{1}(args...)


function Base.getindex(A::FiniteDiffDeriv{T}, i::Int, j::Int) where {T}
    s = A.stencil_offset + 1

    # interior points
    if s <= i <= (A.len-s+1)
        idx  = j - i + s
        coeffs = A.stencil_coefs
        ss     = s
        stencil_length = A.stencil_length
    elseif i < s
        coeffs = A.low_boundary_coefs[i]
        ss     = A.low_boundary_offsets[i] + 1
        idx    = j - i + ss
        stencil_length = A.boundary_stencil_length
    else # i > (A.len-s+1)
        ii     = A.len-i+1
        coeffs = A.high_boundary_coefs[ii]
        ss     = A.high_boundary_offsets[ii] + 1
        idx    = j - i + ss
        stencil_length = A.boundary_stencil_length
    end

    if idx < 1 || idx > stencil_length
        return zero(T)
    else
        return coeffs[idx]
    end
end

function Base.getindex(A::PeriodicFD{T}, i::Int, j::Int) where {T}
    N      = A.len
    coeffs = A.stencil_coefs
    s      = A.stencil_offset + 1

    if 1 <= j - i + s <= N
        idx  = j - i + s
    else
        # imposing periodicity
        idx  = mod1(j - i + s, N)
    end

    if idx < 1 || idx > A.stencil_length
        return zero(T)
    else
        return coeffs[idx]
    end
end


@inline Base.getindex(A::AbstractFiniteDiff, i::Int, ::Colon) =
    [A[i,j] for j in 1:A.len]

@inline Base.getindex(A::SpectralDeriv, i, j) = A.D[i,j]


# make it a callable struct, to compute derivatives at a given point
function (A::AbstractFiniteDiff{T,N})(f::AbstractArray{T,M},
                                      idx::Vararg{Int,M}) where {T<:Real,N,M}
    # make sure axis of differentiation is contained in the dimensions of f
    @assert N <= M

    s = A.stencil_offset + 1
    i = idx[N] # point where derivative will be taken (with respect to the N-axis)

    if s <= i <= (A.len-s+1)
        return _D_interior(A, f, idx)
    elseif i < s
        return _D_low_bdr(A, f, idx)
    else # i > (A.len-s+1)
        return _D_high_bdr(A, f, idx)
    end
end

# interior points
function _D_interior(A::AbstractFiniteDiff{T,N}, f::AbstractArray, idx) where {T<:Real,N}
    coeffs = A.stencil_coefs
    s      = A.stencil_offset + 1
    i      = idx[N] # point where derivative will be taken (with respect to the N-axis)

    sum_i  = zero(T)
    @fastmath @inbounds for aa in 1:A.stencil_length
        i_circ = i - (s - aa)
        I      = Base.setindex(idx, i_circ, N)
        sum_i += coeffs[aa] * f[I...]
    end

    sum_i
end

# boundary points
function _D_low_bdr(A::PeriodicFD{T,N}, f::AbstractArray, idx) where {T<:Real,N}
    coeffs = A.stencil_coefs
    s      = A.stencil_offset + 1
    i      = idx[N]

    sum_i  = zero(T)
    @fastmath @inbounds for aa in 1:A.stencil_length
        # imposing periodicity
        i_circ = mod1(i - (s-aa), A.len)
        I      = Base.setindex(idx, i_circ, N)
        sum_i += coeffs[aa] * f[I...]
    end

    sum_i
end
function _D_high_bdr(A::PeriodicFD{T,N}, f::AbstractArray, idx) where {T<:Real,N}
    coeffs = A.stencil_coefs
    s      = A.stencil_offset + 1
    i      = idx[N]

    sum_i  = zero(T)
    @fastmath @inbounds for aa in 1:A.stencil_length
        # imposing periodicity
        i_circ = mod1(i - (s-aa), A.len)
        I      = Base.setindex(idx, i_circ, N)
        sum_i += coeffs[aa] * f[I...]
    end

    sum_i
end

# low boundary points
function _D_low_bdr(A::FiniteDiffDeriv{T,N}, f::AbstractArray, idx) where {T<:Real,N}
    i      = idx[N]

    sum_i  = zero(T)
    @fastmath @inbounds for aa in 1:A.boundary_stencil_length
        coeffs = A.low_boundary_coefs[i]
        s      = A.low_boundary_offsets[i] + 1
        i_circ = i - (s - aa)
        I      = Base.setindex(idx, i_circ, N)
        sum_i += coeffs[aa] * f[I...]
    end

    sum_i
end

# high boundary points
function _D_high_bdr(A::FiniteDiffDeriv{T,N}, f::AbstractArray, idx) where {T<:Real,N}
    i      = idx[N]
    ii     = A.len-i+1

    sum_i  = zero(T)
    @fastmath @inbounds for aa in 1:A.boundary_stencil_length
        coeffs = A.high_boundary_coefs[ii]
        s      = A.high_boundary_offsets[ii] + 1
        i_circ = i - (s - aa)
        I      = Base.setindex(idx, i_circ, N)
        sum_i += coeffs[aa] * f[I...]
    end

    sum_i
end


function LinearAlgebra.mul!(df::AbstractVector, A::AbstractFiniteDiff, f::AbstractVector)
    @fastmath @inbounds for idx in eachindex(f)
        df[idx] = A(f,idx)
    end
    nothing
end

function LinearAlgebra.mul!(df::AbstractArray, A::AbstractFiniteDiff, f::AbstractArray)
    @fastmath @inbounds for idx in CartesianIndices(f)
        df[idx] = A(f,idx.I...)
    end
    nothing
end



# mul! done by standard matrix multiplication for Chebyshev differentiation
# matrices. This can also be done through FFT, but after some testing the FFT
# route actually seemed to be performing slower, so let's stick with this
function LinearAlgebra.mul!(df::AbstractVector, A::SpectralDeriv, f::AbstractVector)
    mul!(df, A.D, f)
end

# and now for Arrays
function LinearAlgebra.mul!(df::AbstractArray{T}, A::SpectralDeriv{T,N,S},
                            f::AbstractArray{T}) where {T,N,S}
    Rpre  = CartesianIndices(axes(f)[1:N-1])
    Rpost = CartesianIndices(axes(f)[N+1:end])

    _mul_loop!(df, A, f, Rpre, size(f,N), Rpost)
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
@noinline function _mul_loop!(df::AbstractArray{T}, A, f, Rpre, n, Rpost) where {T}
    @fastmath @inbounds for Ipost in Rpost
        @inbounds for i in 1:n
            @inbounds for Ipre in Rpre
                sum_i = zero(T)
                @inbounds @simd for j in 1:n
                    sum_i += A.D[i,j] * f[Ipre,j,Ipost]
                end
                df[Ipre,i,Ipost] = sum_i
            end
        end
    end
    nothing
end

function *(A::AbstractDerivOperator, x::AbstractArray)
    y = similar(x)
    mul!(y, A, x)
    y
end


# make SpectralDeriv a callable struct, to compute derivatives only at a given point.
function (A::SpectralDeriv{T,N,S})(f::AbstractArray{T,M},
                                   idx::Vararg{Int,M}) where {T<:Real,N,M,S}
    # make sure axis of differentiation is contained in the dimensions of f
    @assert N <= M

    i  = idx[N] # point where derivative will be taken (with respect to the N-axis)

    sum_i = zero(T)
    @fastmath @inbounds for ii in 1:A.len
        I = Base.setindex(idx, ii, N)
        sum_i += A.D[i,ii] * f[I...]
    end
    sum_i
end


# now for cross-derivatives. we assume that A acts on the first and B on the
# second axis of the x Matrix.
function (A::AbstractFiniteDiff{T,N1})(B::AbstractFiniteDiff{T,N2},
                                       f::AbstractArray{T,M},
                                       idx::Vararg{Int,M}) where {T<:Real,N1,N2,M}
    NA   = A.len
    NB   = B.len
    sA   = A.stencil_offset + 1
    sB   = B.stencil_offset + 1

    # make sure axes of differentiation are contained in the dimensions of f
    @assert N1 < N2 <= M

    # points where derivative will be taken (along their respective axes)
    i  = idx[N1]
    j  = idx[N2]

    if sA <= i <= NA-sA+1 && sB <= j <= NB-sB+1
        return _D_interior(A, B, f, idx)

    elseif sA <= i <= NA-sA+1 && j < sB
        return _D_intA_lowB(A, B, f, idx)

    elseif sA <= i <= NA-sA+1 && j > NB-sB+1
        return _D_intA_highB(A, B, f, idx)

    elseif i < sA && sB <= j <= NB-sB+1
        return _D_lowA_intB(A, B, f, idx)

    elseif i > NA-sA+1 && sB <= j <= NB-sB+1
        return _D_highA_intB(A, B, f, idx)

    elseif i < sA && j < sB
        return _D_lowA_lowB(A, B, f, idx)

    elseif i < sA && j > NB-sB+1
        return _D_lowA_highB(A, B, f, idx)

    elseif i > NA-sA+1 && j < sB
        return _D_highA_lowB(A, B, f, idx)

    elseif i > NA-sA+1 && j > NB-sB+1
        return _D_highA_highB(A, B, f, idx)

    else
        error("I shouldn't be here...")
    end

end

# interior points
function _D_interior(A::AbstractFiniteDiff{T,N1}, B::AbstractFiniteDiff{T,N2},
                     f::AbstractArray, idx) where {T<:Real,N1,N2}
    NA   = A.len
    NB   = B.len
    qA   = A.stencil_coefs
    qB   = B.stencil_coefs
    sA   = A.stencil_offset + 1
    sB   = B.stencil_offset + 1

    # points where derivative will be taken (along their respective axes)
    i  = idx[N1]
    j  = idx[N2]

    sum_ij = zero(T)
    @fastmath @inbounds for jj in 1:B.stencil_length
        j_circ = j - (sB-jj)
        Itmp   = Base.setindex(idx, j_circ, N2)
        sum_i  = zero(T)
        @inbounds for ii in 1:A.stencil_length
            i_circ = i - (sA-ii)
            I      = Base.setindex(Itmp, i_circ, N1)
            sum_i += qA[ii] * qB[jj] * f[I...]
        end
        sum_ij += sum_i
    end

    sum_ij
end

# boundary points
function _D_intA_lowB(A::AbstractFiniteDiff{T,N1}, B::PeriodicFD{T,N2},
                      f::AbstractArray, idx) where {T<:Real,N1,N2}
    NA   = A.len
    NB   = B.len
    qA   = A.stencil_coefs
    qB   = B.stencil_coefs
    sA   = A.stencil_offset + 1
    sB   = B.stencil_offset + 1

    # points where derivative will be taken (along their respective axes)
    i  = idx[N1]
    j  = idx[N2]

    sum_ij = zero(T)
    @fastmath @inbounds for jj in 1:B.stencil_length
        # imposing periodicity
        j_circ = mod1(j - (sB-jj), NB)
        Itmp   = Base.setindex(idx, j_circ, N2)
        sum_i  = zero(T)
        @inbounds for ii in 1:A.stencil_length
            i_circ = i - (sA-ii)
            I      = Base.setindex(Itmp, i_circ, N1)
            sum_i += qA[ii] * qB[jj] * f[I...]
        end
        sum_ij += sum_i
    end

    sum_ij
end
function _D_intA_highB(A::AbstractFiniteDiff{T,N1}, B::PeriodicFD{T,N2},
                       f::AbstractArray, idx) where {T<:Real,N1,N2}
    NA   = A.len
    NB   = B.len
    qA   = A.stencil_coefs
    qB   = B.stencil_coefs
    sA   = A.stencil_offset + 1
    sB   = B.stencil_offset + 1

    # points where derivative will be taken (along their respective axes)
    i  = idx[N1]
    j  = idx[N2]

    sum_ij = zero(T)
    @fastmath @inbounds for jj in 1:B.stencil_length
        # imposing periodicity
        j_circ = mod1(j - (sB-jj), NB)
        Itmp   = Base.setindex(idx, j_circ, N2)
        sum_i  = zero(T)
        @inbounds for ii in 1:A.stencil_length
            i_circ = i - (sA-ii)
            I      = Base.setindex(Itmp, i_circ, N1)
            sum_i += qA[ii] * qB[jj] * f[I...]
        end
        sum_ij += sum_i
    end

    sum_ij
end
function _D_lowA_intB(A::PeriodicFD{T,N1}, B::AbstractFiniteDiff{T,N2},
                      f::AbstractArray, idx) where {T<:Real,N1,N2}
    NA   = A.len
    NB   = B.len
    qA   = A.stencil_coefs
    qB   = B.stencil_coefs
    sA   = A.stencil_offset + 1
    sB   = B.stencil_offset + 1

    # points where derivative will be taken (along their respective axes)
    i  = idx[N1]
    j  = idx[N2]

    sum_ij = zero(T)
    @fastmath @inbounds for jj in 1:B.stencil_length
        j_circ = j - (sB-jj)
        Itmp   = Base.setindex(idx, j_circ, N2)
        sum_i  = zero(T)
        @inbounds for ii in 1:A.stencil_length
            # imposing periodicity
            i_circ = mod1(i - (sA-ii), NA)
            I      = Base.setindex(Itmp, i_circ, N1)
            sum_i += qA[ii] * qB[jj] * f[I...]
        end
        sum_ij += sum_i
    end

    sum_ij
end
function _D_highA_intB(A::PeriodicFD{T,N1}, B::AbstractFiniteDiff{T,N2},
                       f::AbstractArray, idx) where {T<:Real,N1,N2}
    NA   = A.len
    NB   = B.len
    qA   = A.stencil_coefs
    qB   = B.stencil_coefs
    sA   = A.stencil_offset + 1
    sB   = B.stencil_offset + 1

    # points where derivative will be taken (along their respective axes)
    i  = idx[N1]
    j  = idx[N2]

    sum_ij = zero(T)
    @fastmath @inbounds for jj in 1:B.stencil_length
        j_circ = j - (sB-jj)
        Itmp   = Base.setindex(idx, j_circ, N2)
        sum_i  = zero(T)
        @inbounds for ii in 1:A.stencil_length
            # imposing periodicity
            i_circ = mod1(i - (sA-ii), NA)
            I      = Base.setindex(Itmp, i_circ, N1)
            sum_i += qA[ii] * qB[jj] * f[I...]
        end
        sum_ij += sum_i
    end

    sum_ij
end
function _D_lowA_lowB(A::PeriodicFD{T,N1}, B::PeriodicFD{T,N2},
                      f::AbstractArray, idx) where {T<:Real,N1,N2}
    NA   = A.len
    NB   = B.len
    qA   = A.stencil_coefs
    qB   = B.stencil_coefs
    sA   = A.stencil_offset + 1
    sB   = B.stencil_offset + 1

    # points where derivative will be taken (along their respective axes)
    i  = idx[N1]
    j  = idx[N2]

    sum_ij = zero(T)
    @fastmath @inbounds for jj in 1:B.stencil_length
        # imposing periodicity
        j_circ = mod1(j - (sB-jj), NB)
        Itmp   = Base.setindex(idx, j_circ, N2)
        sum_i  = zero(T)
        @inbounds for ii in 1:A.stencil_length
            # imposing periodicity
            i_circ = mod1(i - (sA-ii), NA)
            I      = Base.setindex(Itmp, i_circ, N1)
            sum_i += qA[ii] * qB[jj] * f[I...]
        end
        sum_ij += sum_i
    end

    sum_ij
end
function _D_lowA_highB(A::PeriodicFD{T,N1}, B::PeriodicFD{T,N2},
                       f::AbstractArray, idx) where {T<:Real,N1,N2}
    NA   = A.len
    NB   = B.len
    qA   = A.stencil_coefs
    qB   = B.stencil_coefs
    sA   = A.stencil_offset + 1
    sB   = B.stencil_offset + 1

    # points where derivative will be taken (along their respective axes)
    i  = idx[N1]
    j  = idx[N2]

    sum_ij = zero(T)
    @fastmath @inbounds for jj in 1:B.stencil_length
        # imposing periodicity
        j_circ = mod1(j - (sB-jj), NB)
        Itmp   = Base.setindex(idx, j_circ, N2)
        sum_i  = zero(T)
        @inbounds for ii in 1:A.stencil_length
            # imposing periodicity
            i_circ = mod1(i - (sA-ii), NA)
            I      = Base.setindex(Itmp, i_circ, N1)
            sum_i += qA[ii] * qB[jj] * f[I...]
        end
        sum_ij += sum_i
    end

    sum_ij
end
function _D_highA_lowB(A::PeriodicFD{T,N1}, B::PeriodicFD{T,N2},
                       f::AbstractArray, idx) where {T<:Real,N1,N2}
    NA   = A.len
    NB   = B.len
    qA   = A.stencil_coefs
    qB   = B.stencil_coefs
    sA   = A.stencil_offset + 1
    sB   = B.stencil_offset + 1

    # points where derivative will be taken (along their respective axes)
    i  = idx[N1]
    j  = idx[N2]

    sum_ij = zero(T)
    @fastmath @inbounds for jj in 1:B.stencil_length
        # imposing periodicity
        j_circ = mod1(j - (sB-jj), NB)
        Itmp   = Base.setindex(idx, j_circ, N2)
        sum_i  = zero(T)
        @inbounds for ii in 1:A.stencil_length
            # imposing periodicity
            i_circ = mod1(i - (sA-ii), NA)
            I      = Base.setindex(Itmp, i_circ, N1)
            sum_i += qA[ii] * qB[jj] * f[I...]
        end
        sum_ij += sum_i
    end

    sum_ij
end
function _D_highA_highB(A::PeriodicFD{T,N1}, B::PeriodicFD{T,N2},
                        f::AbstractArray, idx) where {T<:Real,N1,N2}
    NA   = A.len
    NB   = B.len
    qA   = A.stencil_coefs
    qB   = B.stencil_coefs
    sA   = A.stencil_offset + 1
    sB   = B.stencil_offset + 1

    # points where derivative will be taken (along their respective axes)
    i  = idx[N1]
    j  = idx[N2]

    sum_ij = zero(T)
    @fastmath @inbounds for jj in 1:B.stencil_length
        # imposing periodicity
        j_circ = mod1(j - (sB-jj), NB)
        Itmp   = Base.setindex(idx, j_circ, N2)
        sum_i  = zero(T)
        @inbounds for ii in 1:A.stencil_length
            # imposing periodicity
            i_circ = mod1(i - (sA-ii), NA)
            I      = Base.setindex(Itmp, i_circ, N1)
            sum_i += qA[ii] * qB[jj] * f[I...]
        end
        sum_ij += sum_i
    end

    sum_ij
end


# now for FiniteDiffDeriv structs

function _D_intA_lowB(A::AbstractFiniteDiff{T,N1}, B::FiniteDiffDeriv{T,N2},
                      f::AbstractArray, idx) where {T<:Real,N1,N2}
    NA   = A.len
    NB   = B.len
    qA   = A.stencil_coefs
    sA   = A.stencil_offset + 1

    # points where derivative will be taken (along their respective axes)
    i  = idx[N1]
    j  = idx[N2]

    qB  = B.low_boundary_coefs[j]
    sB  = B.low_boundary_offsets[j] + 1

    sum_ij = zero(T)
    @fastmath @inbounds for jj in 1:B.boundary_stencil_length
        j_circ = j - (sB-jj)
        Itmp   = Base.setindex(idx, j_circ, N2)
        sum_i  = zero(T)
        @inbounds for ii in 1:A.stencil_length
            i_circ = i - (sA-ii)
            I      = Base.setindex(Itmp, i_circ, N1)
            sum_i += qA[ii] * qB[jj] * f[I...]
        end
        sum_ij += sum_i
    end

    sum_ij
end

function _D_intA_highB(A::AbstractFiniteDiff{T,N1}, B::FiniteDiffDeriv{T,N2},
                       f::AbstractArray, idx) where {T<:Real,N1,N2}
    NB   = B.len
    qA   = A.stencil_coefs
    sA   = A.stencil_offset + 1

    i     = idx[N1]
    j     = idx[N2]
    jhigh = NB-j+1

    qB    = B.high_boundary_coefs[jhigh]
    sB    = B.high_boundary_offsets[jhigh] + 1

    sum_ij = zero(T)
    @fastmath @inbounds for jj in 1:B.boundary_stencil_length
        j_circ = j - (sB-jj)
        Itmp   = Base.setindex(idx, j_circ, N2)
        sum_i  = zero(T)
        @inbounds for ii in 1:A.stencil_length
            i_circ = i - (sA-ii)
            I      = Base.setindex(Itmp, i_circ, N1)
            sum_i += qA[ii] * qB[jj] * f[I...]
        end
        sum_ij += sum_i
    end

    sum_ij
end

function _D_lowA_intB(A::FiniteDiffDeriv{T,N1}, B::AbstractFiniteDiff{T,N2},
                      f::AbstractArray, idx) where {T<:Real,N1,N2}
    qB = B.stencil_coefs
    sB = B.stencil_offset + 1

    i  = idx[N1]
    j  = idx[N2]

    qA = A.low_boundary_coefs[i]
    sA = A.low_boundary_offsets[i] + 1

    sum_ij = zero(T)
    @fastmath @inbounds for jj in 1:B.stencil_length
        j_circ = j - (sB-jj)
        Itmp   = Base.setindex(idx, j_circ, N2)
        sum_i  = zero(T)
        @inbounds for ii in 1:A.boundary_stencil_length
            i_circ = i - (sA-ii)
            I      = Base.setindex(Itmp, i_circ, N1)
            sum_i += qA[ii] * qB[jj] * f[I...]
        end
        sum_ij += sum_i
    end

    sum_ij
end

function _D_highA_intB(A::FiniteDiffDeriv{T,N1}, B::AbstractFiniteDiff{T,N2},
                       f::AbstractArray, idx) where {T<:Real,N1,N2}
    NA   = A.len
    qB   = B.stencil_coefs
    sB   = B.stencil_offset + 1

    i     = idx[N1]
    j     = idx[N2]
    ihigh = NA-i+1

    qA    = A.high_boundary_coefs[ihigh]
    sA    = A.high_boundary_offsets[ihigh] + 1

    sum_ij = zero(T)
    @fastmath @inbounds for jj in 1:B.stencil_length
        j_circ = j - (sB-jj)
        Itmp   = Base.setindex(idx, j_circ, N2)
        sum_i  = zero(T)
        @inbounds for ii in 1:A.boundary_stencil_length
            i_circ = i - (sA-ii)
            I      = Base.setindex(Itmp, i_circ, N1)
            sum_i += qA[ii] * qB[jj] * f[I...]
        end
        sum_ij += sum_i
    end

    sum_ij
end

function _D_lowA_lowB(A::FiniteDiffDeriv{T,N1}, B::FiniteDiffDeriv{T,N2},
                      f::AbstractArray, idx) where {T<:Real,N1,N2}
    i  = idx[N1]
    j  = idx[N2]

    qA = A.low_boundary_coefs[i]
    sA = A.low_boundary_offsets[i] + 1
    qB = B.low_boundary_coefs[j]
    sB = B.low_boundary_offsets[j] + 1

    sum_ij = zero(T)
    @fastmath @inbounds for jj in 1:B.boundary_stencil_length
        j_circ = j - (sB-jj)
        Itmp   = Base.setindex(idx, j_circ, N2)
        sum_i  = zero(T)
        @inbounds for ii in 1:A.boundary_stencil_length
            i_circ = i - (sA-ii)
            I      = Base.setindex(Itmp, i_circ, N1)
            sum_i += qA[ii] * qB[jj] * f[I...]
        end
        sum_ij += sum_i
    end

    sum_ij
end

function _D_lowA_highB(A::FiniteDiffDeriv{T,N1}, B::FiniteDiffDeriv{T,N2},
                       f::AbstractArray, idx) where {T<:Real,N1,N2}
    NB = B.len

    i     = idx[N1]
    j     = idx[N2]
    jhigh = NB-j+1

    qA    = A.low_boundary_coefs[i]
    sA    = A.low_boundary_offsets[i] + 1
    qB    = B.high_boundary_coefs[jhigh]
    sB    = B.high_boundary_offsets[jhigh] + 1

    sum_ij = zero(T)
    @fastmath @inbounds for jj in 1:B.boundary_stencil_length
        j_circ = j - (sB-jj)
        Itmp   = Base.setindex(idx, j_circ, N2)
        sum_i  = zero(T)
        @inbounds for ii in 1:A.boundary_stencil_length
            i_circ = i - (sA-ii)
            I      = Base.setindex(Itmp, i_circ, N1)
            sum_i += qA[ii] * qB[jj] * f[I...]
        end
        sum_ij += sum_i
    end

    sum_ij
end

function _D_highA_lowB(A::FiniteDiffDeriv{T,N1}, B::FiniteDiffDeriv{T,N2},
                       f::AbstractArray, idx) where {T<:Real,N1,N2}
    NA = A.len

    i     = idx[N1]
    j     = idx[N2]
    ihigh = NA-i+1

    qA    = A.high_boundary_coefs[ihigh]
    sA    = A.high_boundary_offsets[ihigh] + 1
    qB    = B.low_boundary_coefs[j]
    sB    = B.low_boundary_offsets[j] + 1

    sum_ij = zero(T)
    @fastmath @inbounds for jj in 1:B.boundary_stencil_length
        j_circ = j - (sB-jj)
        Itmp   = Base.setindex(idx, j_circ, N2)
        sum_i  = zero(T)
        @inbounds for ii in 1:A.boundary_stencil_length
            i_circ = i - (sA-ii)
            I      = Base.setindex(Itmp, i_circ, N1)
            sum_i += qA[ii] * qB[jj] * f[I...]
        end
        sum_ij += sum_i
    end

    sum_ij
end

function _D_highA_highB(A::FiniteDiffDeriv{T,N1}, B::FiniteDiffDeriv{T,N2},
                        f::AbstractArray, idx) where {T<:Real,N1,N2}
    NA = A.len
    NB = B.len

    i     = idx[N1]
    j     = idx[N2]
    ihigh = NA-i+1
    jhigh = NB-j+1

    qA    = A.high_boundary_coefs[ihigh]
    sA    = A.high_boundary_offsets[ihigh] + 1
    qB    = B.high_boundary_coefs[jhigh]
    sB    = B.high_boundary_offsets[jhigh] + 1

    sum_ij = zero(T)
    @fastmath @inbounds for jj in 1:B.boundary_stencil_length
        j_circ = j - (sB-jj)
        Itmp   = Base.setindex(idx, j_circ, N2)
        sum_i  = zero(T)
        @inbounds for ii in 1:A.boundary_stencil_length
            i_circ = i - (sA-ii)
            I      = Base.setindex(Itmp, i_circ, N1)
            sum_i += qA[ii] * qB[jj] * f[I...]
        end
        sum_ij += sum_i
    end

    sum_ij
end


# cross derivatives with FiniteDiffDeriv and PeriodicFD

function _D_lowA_lowB(A::PeriodicFD{T,N1}, B::FiniteDiffDeriv{T,N2},
                      f::AbstractArray, idx) where {T<:Real,N1,N2}
    NA = A.len
    qA = A.stencil_coefs
    sA = A.stencil_offset + 1

    i  = idx[N1]
    j  = idx[N2]

    qB = B.low_boundary_coefs[j]
    sB = B.low_boundary_offsets[j] + 1

    sum_ij = zero(T)
    @fastmath @inbounds for jj in 1:B.boundary_stencil_length
        j_circ = j - (sB-jj)
        Itmp   = Base.setindex(idx, j_circ, N2)
        sum_i  = zero(T)
        @inbounds for ii in 1:A.stencil_length
            # imposing periodicity
            i_circ = mod1(i - (sA-ii), NA)
            I      = Base.setindex(Itmp, i_circ, N1)
            sum_i += qA[ii] * qB[jj] * f[I...]
        end
        sum_ij += sum_i
    end

    sum_ij
end
function _D_lowA_lowB(A::FiniteDiffDeriv{T,N1}, B::PeriodicFD{T,N2},
                      f::AbstractArray, idx) where {T<:Real,N1,N2}
    NB = B.len
    qB = B.stencil_coefs
    sB = B.stencil_offset + 1

    i  = idx[N1]
    j  = idx[N2]

    qA = A.low_boundary_coefs[i]
    sA = A.low_boundary_offsets[i] + 1

    sum_ij = zero(T)
    @fastmath @inbounds for jj in 1:B.stencil_length
        # imposing periodicity
        j_circ = mod1(j - (sB-jj), NB)
        Itmp   = Base.setindex(idx, j_circ, N2)
        sum_i  = zero(T)
        @inbounds for ii in 1:A.boundary_stencil_length
            i_circ = i - (sA-ii)
            I      = Base.setindex(Itmp, i_circ, N1)
            sum_i += qA[ii] * qB[jj] * f[I...]
        end
        sum_ij += sum_i
    end

    sum_ij
end

function _D_lowA_highB(A::PeriodicFD{T,N1}, B::FiniteDiffDeriv{T,N2},
                       f::AbstractArray, idx) where {T<:Real,N1,N2}
    NA    = A.len
    qA    = A.stencil_coefs
    sA    = A.stencil_offset + 1

    i     = idx[N1]
    j     = idx[N2]
    jhigh = NB-j+1

    qB    = B.high_boundary_coefs[jhigh]
    sB    = B.high_boundary_offsets[jhigh] + 1

    sum_ij = zero(T)
    @fastmath @inbounds for jj in 1:B.boundary_stencil_length
        j_circ = j - (sB-jj)
        Itmp   = Base.setindex(idx, j_circ, N2)
        sum_i  = zero(T)
        @inbounds for ii in 1:A.stencil_length
            # imposing periodicity
            i_circ = mod1(i - (sA-ii), NA)
            I      = Base.setindex(Itmp, i_circ, N1)
            sum_i += qA[ii] * qB[jj] * f[I...]
        end
        sum_ij += sum_i
    end

    sum_ij
end
function _D_lowA_highB(A::FiniteDiffDeriv{T,N1}, B::PeriodicFD{T,N2},
                       f::AbstractArray, idx) where {T<:Real,N1,N2}
    NB    = B.len
    qB    = B.stencil_coefs
    sB    = B.stencil_offset + 1

    i     = idx[N1]
    j     = idx[N2]
    jhigh = NB-j+1

    qA    = A.low_boundary_coefs[i]
    sA    = A.low_boundary_offsets[i] + 1

    sum_ij = zero(T)
    @fastmath @inbounds for jj in 1:B.stencil_length
        # imposing periodicity
        j_circ = mod1(j - (sB-jj), NB)
        Itmp   = Base.setindex(idx, j_circ, N2)
        sum_i  = zero(T)
        @inbounds for ii in 1:A.boundary_stencil_length
            i_circ = i - (sA-ii)
            I      = Base.setindex(Itmp, i_circ, N1)
            sum_i += qA[ii] * qB[jj] * f[I...]
        end
        sum_ij += sum_i
    end

    sum_ij
end

function _D_highA_lowB(A::PeriodicFD{T,N1}, B::FiniteDiffDeriv{T,N2},
                       f::AbstractArray, idx) where {T<:Real,N1,N2}
    NA    = A.len
    qA    = A.stencil_coefs
    sA    = A.stencil_offset + 1

    i     = idx[N1]
    j     = idx[N2]
    ihigh = NA-i+1

    qB    = B.low_boundary_coefs[j]
    sB    = B.low_boundary_offsets[j] + 1

    sum_ij = zero(T)
    @fastmath @inbounds for jj in 1:B.boundary_stencil_length
        j_circ = j - (sB-jj)
        Itmp   = Base.setindex(idx, j_circ, N2)
        sum_i  = zero(T)
        @inbounds for ii in 1:A.stencil_length
            # imposing periodicity
            i_circ = mod1(i - (sA-ii), NA)
            I      = Base.setindex(Itmp, i_circ, N1)
            sum_i += qA[ii] * qB[jj] * f[I...]
        end
        sum_ij += sum_i
    end

    sum_ij
end
function _D_highA_lowB(A::FiniteDiffDeriv{T,N1}, B::PeriodicFD{T,N2},
                       f::AbstractArray, idx) where {T<:Real,N1,N2}
    NA    = A.len
    NB    = B.len
    qB    = B.stencil_coefs
    sB    = B.stencil_offset + 1

    i     = idx[N1]
    j     = idx[N2]
    ihigh = NA-i+1

    qA    = A.high_boundary_coefs[ihigh]
    sA    = A.high_boundary_offsets[ihigh] + 1

    sum_ij = zero(T)
    @fastmath @inbounds for jj in 1:B.stencil_length
        # imposing periodicity
        j_circ = mod1(j - (sB-jj), NB)
        Itmp   = Base.setindex(idx, j_circ, N2)
        sum_i  = zero(T)
        @inbounds for ii in 1:A.boundary_stencil_length
            i_circ = i - (sA-ii)
            I      = Base.setindex(Itmp, i_circ, N1)
            sum_i += qA[ii] * qB[jj] * f[I...]
        end
        sum_ij += sum_i
    end

    sum_ij
end

function _D_highA_highB(A::PeriodicFD{T,N1}, B::FiniteDiffDeriv{T,N2},
                        f::AbstractArray, idx) where {T<:Real,N1,N2}
    NA    = A.len
    NB    = B.len

    qA    = A.stencil_coefs
    sA    = A.stencil_offset + 1

    i     = idx[N1]
    j     = idx[N2]
    ihigh = NA-i+1
    jhigh = NB-j+1

    qB    = B.high_boundary_coefs[jhigh]
    sB    = B.high_boundary_offsets[jhigh] + 1

    sum_ij = zero(T)
    @fastmath @inbounds for jj in 1:B.boundary_stencil_length
        j_circ = j - (sB-jj)
        Itmp   = Base.setindex(idx, j_circ, N2)
        sum_i  = zero(T)
        @inbounds for ii in 1:A.stencil_length
            # imposing periodicity
            i_circ = mod1(i - (sA-ii), NA)
            I      = Base.setindex(Itmp, i_circ, N1)
            sum_i += qA[ii] * qB[jj] * f[I...]
        end
        sum_ij += sum_i
    end

    sum_ij
end
function _D_highA_highB(A::FiniteDiffDeriv{T,N1}, B::PeriodicFD{T,N2},
                        f::AbstractArray, idx) where {T<:Real,N1,N2}
    NA    = A.len
    NB    = B.len
    qB    = B.stencil_coefs
    sB    = B.stencil_offset + 1

    i     = idx[N1]
    j     = idx[N2]
    ihigh = NA-i+1
    jhigh = NB-j+1

    qA    = A.high_boundary_coefs[ihigh]
    sA    = A.high_boundary_offsets[ihigh] + 1

    sum_ij = zero(T)
    @fastmath @inbounds for jj in 1:B.stencil_length
        # imposing periodicity
        j_circ = mod1(j - (sB-jj), NB)
        Itmp   = Base.setindex(idx, j_circ, N2)
        sum_i  = zero(T)
        @inbounds for ii in 1:A.boundary_stencil_length
            i_circ = i - (sA-ii)
            I      = Base.setindex(Itmp, i_circ, N1)
            sum_i += qA[ii] * qB[jj] * f[I...]
        end
        sum_ij += sum_i
    end

    sum_ij
end


# Casting to matrix types

copyto!(M::AbstractMatrix{T}, A::SpectralDeriv) where {T<:Real} =
    copyto!(M, A.D)

function copyto!(M::AbstractMatrix{T}, A::AbstractFiniteDiff) where {T<:Real}
    for idx in CartesianIndices(M)
        M[idx] = A[idx.I...]
    end
    M
end

LinearAlgebra.Array(A::AbstractDerivOperator{T}) where {T} =
    copyto!(zeros(T, A.len, A.len), A)

SparseArrays.SparseMatrixCSC(A::AbstractFiniteDiff{T}) where {T} =
    copyto!(spzeros(T, A.len, A.len), A)

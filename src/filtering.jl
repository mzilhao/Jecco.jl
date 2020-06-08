
abstract type Filter{T,N} end

function KO_kernel(N::Int, sigma_diss::T) where {T<:Real}
    if N == 3
        kernel = -T[1, -4, 6, -4, 1]
        fac = sigma_diss / 16
    elseif N == 5
        kernel = T[1, -6, 15, -20, 15, -6, 1]
        fac = sigma_diss / 64
    else
        error("N = $N Kreiss-Oliger filter not yet implemented.")
    end

    rmul!(kernel, fac)

    # adding one to the middle entry so that it becomes a low-pass filter
    kernel[N+1] += 1

    kernel
end

struct ConvFilter{T<:Real,M} <: Filter{T,M}
    kernel :: Vector{T}
end

struct KO_Filter{M} end

function KO_Filter{M}(N::Int, sigma_diss::T) where {T<:Real,M}
    kernel = KO_kernel(N, sigma_diss)
    ConvFilter{T,M}(kernel)
end

KO_Filter(args...) = KO_Filter{1}(args...)

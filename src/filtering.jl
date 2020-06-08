
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
    mid = div(length(kernel) + 1, 2)
    kernel[mid] += 1

    kernel
end

struct ConvFilter{T<:Real,M} <: Filter{T,M}
    # 'len' is not the length of the kernel, but rather the length of the target vectors
    len    :: Int
    kernel :: Vector{T}
    _cache :: Vector{T}
end

struct KO_Filter{M} end

function KO_Filter{M}(N::Int, sigma_diss::T, len::Int) where {T<:Real,M}
    kernel = KO_kernel(N, sigma_diss)
    _cache = Vector{T}(undef, len)
    ConvFilter{T,M}(len, kernel, _cache)
end

KO_Filter(args...) = KO_Filter{1}(args...)

function convolution!(fout::AbstractVector{T}, f::AbstractVector{T},
                      g::AbstractVector{T}) where {T}
    f_len = length(f)
    g_len = length(g)
    mid   = div(g_len + 1, 2)

    @fastmath @inbounds for i in eachindex(f)
        sum_i = zero(T)

        if mid <= i <= (f_len - mid + 1)
            @inbounds for aa in 1:g_len
                i_circ = i - (mid - aa)
                sum_i += g[aa] * f[i_circ]
            end
        else
            @inbounds for aa in 1:g_len
                # imposing periodicity
                i_circ = 1 + mod(i - (mid-aa) - 1, f_len)
                sum_i += g[aa] * f[i_circ]
            end
        end
        fout[i] = sum_i
    end

    fout
end

function (filter::ConvFilter)(f::AbstractVector)
    @assert filter.len == length(f)
    convolution!(filter._cache, f, filter.kernel)
    copyto!(f, filter._cache)
end

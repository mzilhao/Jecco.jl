
"""
    cheb(xmin, xmax, N)

Return a Chebyshev-Lobatto grid, together with its first and second
differentiation matrices.

# Arguments
* `xmin::Real`: rightmost grid point
* `xmax::Real`: leftmost grid point
* `N::Integer`: total number of grid points
"""
function cheb(xmin::T, xmax::T, N::Integer) where {T<:Real}
    x, D, D2 = cheb(N)
    x = 0.5 * (xmax + xmin .+ (xmax - xmin) * x)
    D  ./= 0.5 * (xmax - xmin)
    D2 ./= 0.25 * (xmax - xmin)^2
    x, D, D2
end


# taken from Trefethen (2000), "Spectral Methods in MatLab"

"When omitting grid limits, default to [-1,1] interval"
function cheb(N::Integer)
    @assert(N > 1, "number of points should be greater than 1...")
    M = N - 1
    x = -cos.(pi*(0:M)/M)

    c  = [2; ones(M-1, 1); 2] .* (-1).^(0:M)
    X  = repeat(x, 1, N)
    dX = X - X'

    D = (c * (1 ./ c)') ./ (dX + Matrix(I, N, N))  # off-diagonal entries
    D = D - diagm(0 => sum(D', dims=1)[:])         # diagonal entries

    x, D, D*D
end


"""
    fourier(xsize, N)

Return the first and second derivative matrices on a Fourier grid

# Arguments
* `xsize::Real`: size of (periodic) domain
* `N::Integer`: total number of grid points
"""
function fourier(xsize::T, N::Integer) where {T<:Real}
    D, D2 = fourier(N)
    pi_ = T(pi)
    D  ./= xsize / (2*pi_)
    D2 ./= xsize^2 / (4*pi_^2)
    D, D2
end

# taken from Trefethen (2000), "Spectral Methods in MatLab"

"When omitting `xsize`, default to 2π"
function fourier(N::Integer)
    @assert(mod(N,2)==0, "number of points needs to be even")
    h = 2*pi/N

    column = [0; 0.5*(-1).^(1:N-1) .* cot.((1:N-1)*h/2)]
    tmp = [circshift(column, i) for i in 0:N-1]
    D   = hcat(tmp...)

    column2 = [-pi^2/(3*h^2) - 1/6; -0.5*(-1).^(1:N-1) ./ sin.((1:N-1)*h/2).^2]
    tmp = [circshift(column2, i) for i in 0:N-1]
    D2  = hcat(tmp...)

    D, D2
end


struct ChebInterpolator{T,A,FT<:FFTW.r2rFFTWPlan}
    xmin     :: T
    xmax     :: T
    c        :: A
    fft_plan :: FT
end
"""
    ChebInterpolator(xmin, xmax, N)

Build an interpolator to act on functions defined on a Chebyshev-Lobatto grid.
Uses `FFTW`.

# Arguments
* `xmin::Real`: rightmost grid point
* `xmax::Real`: leftmost grid point
* `N::Integer`: total number of grid points
"""
function ChebInterpolator(xmin::T, xmax::T, N::Int) where {T<:Real}
    M  = N - 1
    x  = Vector{T}(undef, N)

    # use only one thread for the FFTW! by default it allocates a bunch of them,
    # and it hurts performance (by a huge amount!) when using this inside loops.
    # it's much better to thread the loops themselves, which is what we do.
    FFTW.set_num_threads(1)

    # Create the FFT plan for the DCT-I
    fft_plan = FFTW.plan_r2r(x, FFTW.REDFT00)

    c  = M * [2; ones(Int,M-1); 2] .* (-1).^(0:M)

    ChebInterpolator{T,typeof(c),typeof(fft_plan)}(xmin, xmax, c, fft_plan)
end
"""
    ChebInterpolator(xp::Vector)

Build it directly from a `Vector` with the grid points
"""
ChebInterpolator(xp::Vector) = ChebInterpolator(xp[1], xp[end], length(xp))

"""
# Examples

```
julia> xx, = Jecco.cheb(0.0, 2.0, 16);

julia> f = xx.^2;

julia> interp = Jecco.ChebInterpolator(xx);

julia> f_interp = interp(f);

julia> f_interp(0.2)
0.03999999999999987
```
"""
function (interp::ChebInterpolator)(fp)
    # compute the DCT-I of the coefficients
    fft_fp = interp.fft_plan * fp

    # compute the spectral coefficients
    spec_coeff = fft_fp ./ interp.c

    function (x0::T) where {T<:Real}
        @assert interp.xmin <= x0 <= interp.xmax
        X = (2 * x0 - (interp.xmin + interp.xmax)) / (interp.xmax - interp.xmin)
        sum_l = zero(T)
        @fastmath @inbounds for i in LinearIndices(spec_coeff)
            sum_l += spec_coeff[i] * cos( (i-1)*acos(X) )
        end
        sum_l
    end
end

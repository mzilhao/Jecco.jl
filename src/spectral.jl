
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
    D ./= 0.5 * (xmax - xmin)
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
    fourier(xmin, xmax, N)

Return a (periodic) Fourier grid, together with its first and second derivative matrices

# Arguments
* `xmin::Real`: rightmost grid point
* `xmax::Real`: leftmost grid point
* `N::Integer`: total number of grid points
"""
function fourier(xmin::T, xmax::T, N::Integer) where {T<:Real}
    x, D, D2 = fourier(N)
    x = (xmin .+ (xmax - xmin) * x) / (2*pi)
    D  ./= 0.5 * (xmax - xmin) / pi
    D2 ./= 0.25 * (xmax - xmin)^2 / pi^2
    x, D, D2
end

# taken from Trefethen (2000), "Spectral Methods in MatLab"

"When omitting grid limits, default to ]0, 2pi] interval"
function fourier(N::Integer)
    @assert(mod(N,2)==0, "number of points needs to be even")
    h = 2*pi/N
    x = h*(1:N)

    column = [0; 0.5*(-1).^(1:N-1) .* cot.((1:N-1)*h/2)]
    tmp = [circshift(column, i) for i in 0:N-1]
    D   = hcat(tmp...)

    column2 = [-pi^2/(3*h^2) - 1/6; -0.5*(-1).^(1:N-1) ./ sin.((1:N-1)*h/2).^2]
    tmp = [circshift(column2, i) for i in 0:N-1]
    D2  = hcat(tmp...)

    x, D, D2
end

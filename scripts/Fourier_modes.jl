using NumericalIntegration
using LinearAlgebra
using FFTW

#=
This function computes the Fourier modes either integrating "by hand" or using FFTW. We split the comlex modes, associated to the complex exponentials into a and b, associated to the cos and sin functions by a(kx,ky) = 2*Re(f(kx,ky)) and b(kx,ky) = -2*Im(f(kx,ky)).
=#


function Fourier_modes_2D(e :: Array{T,2}, x :: Array{T,1}, y :: Array{T,1}, order_x :: Integer, order_y :: Integer, integration_method :: String) where {T<:Real}
    Lx, Ly = x[end]-x[1], y[end]-y[1]
    xmid, ymid = (x[end]+x[1])/2, (y[end]+y[1])/2
    nx, ny = length(x),length(y)

    if order_x > Integer(floor(nx)) || order_y > Integer(floor(ny))
        print("\n Error, requested order is higher than the number of grid points\n")
        nothing
    end
    
    if integration_method == "trapezoidal"

        a, b = zeros(order_x+1,order_y+1), zeros(order_x+1,order_y+1)
        xx, yy = zeros(nx,ny), zeros(nx,ny)
        #We create the 2D grid
        for i in 1:ny
            xx[:,i] = x
        end
        for i in 1:nx
            yy[i,:] = y
        end

        for my in 0:order_y
            for mx in 0:order_x
                integral_cos, integral_sin = zeros(nx), zeros(nx)
                integrand_cos = 2/(Lx*Ly).*e.*cos.(2*π*(mx/Lx.*(xx.-xmid).+my/Ly.*(yy.-ymid)))
                integrand_sin = 2/(Lx*Ly).*e.*sin.(2*π*(mx/Lx.*(xx.-xmid).+my/Ly.*(yy.-ymid)))
                for i in 1:nx
                    integral_cos[i] = integrate(y,integrand_cos[i,:], Trapezoidal())
                    integral_sin[i] = integrate(y,integrand_sin[i,:], Trapezoidal())
                end
                a[mx+1,my+1] = integrate(x,integral_cos, Trapezoidal())
                b[mx+1,my+1] = integrate(x,integral_sin, Trapezoidal())
            end
        end

    elseif integration_method == "FFTW"

        plan = plan_rfft(e)
        fft_z = (1/(nx*ny).*(plan*e)[1:order_x+1,1:order_y+1])
        a = 2 .*real(fft_z) 
        b = -2 .*imag(fft_z)

    else
        print("\n No such integration method\n")
        nothing
    end
    return a, b
end

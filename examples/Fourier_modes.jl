using NumericalIntegration
using QuadGK
using Interpolations
using FastGaussQuadrature
using LinearAlgebra
#Function to compute Fourier modes in 2D for real funcitons.
#All methods seem to give really bad results for such small amount number of points besides the trapezoidal that is fast
#and has a discrepancy of 10^-6 with respect to mathematica.
#Maybe try our own method? I don't know why these other methods work so badly.
function Fourier_modes_2D(e :: Array{T,2}, x :: Array{T,1}, y :: Array{T,1}, order_x :: Integer, order_y :: Integer, integration_method :: String) where {T<:Real}
    Lx, Ly = x[end]-x[1], y[end]-y[1]
    xmid, ymid = (x[end]+x[1])/2, (y[end]+y[1])/2
    nx, ny = length(x),length(y)
    a, b = zeros(order_x+1,order_y+1), zeros(order_x+1,order_y+1)
    if integration_method == "trapezoidal" || integration_method == "simps"
        if integration_method == "trapezoidal"
            method = Trapezoidal()
        else
            method = SimpsonEven()
        end
        xx, yy = zeros(nx,ny), zeros(nx,ny)
        for i in 1:ny
            xx[:,i] = x
        end
        for i in 1:nx
            yy[i,:] = y
        end

        for mx in 0:order_x
            for my in 0:order_y
                integral_cos, integral_sin = zeros(nx), zeros(nx)
                integrand_cos = 2/(Lx*Ly).*e.*cos.(2*π*(mx/Lx.*(xx.-xmid).+my/Ly.*(yy.-ymid)))
                integrand_sin = 2/(Lx*Ly).*e.*sin.(2*π*(mx/Lx.*(xx.-xmid).+my/Ly.*(yy.-ymid)))
                for i in 1:nx
                    integral_cos[i] = integrate(y,integrand_cos[i,:], method)
                    integral_sin[i] = integrate(y,integrand_sin[i,:], method)
                end
                a[mx+1,my+1] = integrate(x,integral_cos, method)
                b[mx+1,my+1] = integrate(x,integral_sin, method)
            end
        end

    elseif integration_method == "quad"
        for mx in 0:order_x
            for my in 0:order_y
                integral_cos, integral_sin = zeros(nx), zeros(nx)
                for i in 1:nx
                    integrand_cos = 2/(Lx*Ly).*e[i,:].*cos.(2*π*(mx/Lx*(x[i]-xmid).+my/Ly.*(y.-ymid)))
                    integrand_sin = 2/(Lx*Ly).*e[i,:].*sin.(2*π*(mx/Lx*(x[i]-xmid).+my/Ly.*(y.-ymid)))
                    interpolated_cos = interpolate((y,),integrand_cos,Gridded(Linear()))
                    interpolated_sin = interpolate((y,),integrand_sin,Gridded(Linear()))
                    integral_cos[i] = quadgk(interpolated_cos,y[1],y[end])[1]
                    integral_sin[i] = quadgk(interpolated_sin,y[1],y[end])[1]
                end
                interpolated_cos = interpolate((x,),integral_cos,Gridded(Linear()))
                interpolated_sin = interpolate((x,),integral_sin,Gridded(Linear()))
                a[mx+1,my+1] = quadgk(interpolated_cos,x[1],x[end])[1]
                b[mx+1,my+1] = quadgk(interpolated_sin,x[1],x[end])[1]
            end
        end
    elseif integration_method == "gauss"
        num_nodes_x, num_nodes_y = 1000, 1000
        xx, yy, aux = zeros(nx,ny), zeros(nx,ny), zeros(num_nodes_x)
        for i in 1:ny
            xx[:,i] = x
        end
        for i in 1:nx
            yy[i,:] = y
        end
        for mx in 0:order_x
            for my in 0:order_y
                integrand_cos = 2/(Lx*Ly).*e.*cos.(2*π*(mx/Lx.*(xx.-xmid).+my/Ly.*(yy.-ymid)))
                integrand_sin = 2/(Lx*Ly).*e.*sin.(2*π*(mx/Lx.*(xx.-xmid).+my/Ly.*(yy.-ymid)))

                eta, xi = 2 .*(x.-x[1])./Lx.-1, 2 .*(y.-y[1])./Ly.-1
                nodes_eta, w_eta = gausslegendre(num_nodes_x)
                nodes_xi, w_xi = gausslegendre(num_nodes_y)
                interpolated_cos = interpolate((eta,xi,),integrand_cos,Gridded(Linear()))(nodes_eta,nodes_xi)
                interpolated_sin = interpolate((eta,xi,),integrand_sin,Gridded(Linear()))(nodes_eta,nodes_xi)

                mul!(aux, interpolated_cos, w_xi)
                a[mx+1,my+1] = sum(w_eta.*aux)
                mul!(aux, interpolated_sin, w_xi)
                b[mx+1,my+1] = sum(w_eta.*aux)
            end
        end

    else
        print("\n No such integration method\n")
        nothing
    end
    return a, b
end

using Jecco, Jecco.AdS5_3_1
import Base.Threads.@threads
using Plots
gr()


dirname = "/home/mikel/Documents/Jecco.jl/data/bubbles/end_data/"

a4_series = BoundaryTimeSeries(dirname, :a4)
t, x, y   = get_coords(a4_series,:,:,:)
a4        = a4_series[:,:,:]
a4_new_1  = similar(a4)
a4_new_2  = similar(a4)
a4_new_3  = similar(a4)
a4_new_4  = similar(a4)
interp    = AdS5_3_1.xy_interpolator(x, y, tolerance=10^-5)
a4_inter  = interp(a4)
_, Nx, Ny = size(a4)
dx, dy    = (x[2]-x[1], y[2]-y[1])

#xx        = zeros(Nx, Ny)
#yy        = zeros(Nx, Ny)
#=
@fastmath @inbounds @threads for j in 1:Ny
        xx[:,j] = x
        yy[:,j].= y[j]
end
=#
#=
@time @fastmath @inbounds @threads for j in 1:Ny
        for i in 1:Nx
                a4_new_1[:,i,j] = a4[:,i,j]
        end
end

@time a4_new_2[:,:,:] = @view a4[:,:,:]



@time @fastmath @inbounds @threads for j in 1:Ny
        for i in 1:Nx
                a4_new_3[:,i,j] = a4_inter(x[i], y[j])
        end
end

@time a4_new_4[:,:,:] = @view a4_inter(xx, yy)[:,:,:]
=#

Dx = CenteredDiff(1,4,dx,Nx)
Dy = CenteredDiff(1,4,dy,Ny)

x_new = x
y_new = y

Nx_new, Ny_new = (length(x_new), length(y_new))
a4_new_3       = zeros(1,Nx_new, Ny_new)

@time @fastmath @inbounds @threads for j in 1:Ny_new
        for i in 1:Nx_new
                xx = x_new[i]
                yy = y_new[j]
                i0 = findlast(x .<= x_new[i])
                if i0 < Nx i1 = i0 + 1 else i1 = i0 end
                j0 = findlast(y .<= y_new[j])
                if j0 < Ny j1 = j0 + 1 else j1 = j0 end
                if abs(x[i0]-xx) <= abs(x[i1]-xx)
                        x_old = x[i0]
                        ii    = i0
                else
                        x_old = x[i1]
                        ii  = i1
                end
                if abs(y[j0]-yy) <= abs(y[j1]-yy)
                        y_old = y[j0]
                        jj  = j0
                else
                        y_old = y[j1]
                        jj  = j1
                end

                a4_new_3[1,i,j] = a4[1,ii,jj]+Dx(a4[1,:,jj],ii)*(xx-x_old)+
                                                Dy(a4[1,ii,:],jj)*(yy-y_old)
        end
end

plot(x_new, a4_new_3[1,:,1], lw=3)
plot!(x, a4[1,:,1], lw=3)

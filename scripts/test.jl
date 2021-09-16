using Jecco, Jecco.AdS5_3_1
import Base.Threads.@threads
import Base.Threads.@spawn
import Base.Threads.@sync
using Plots
gr()


dirname = "/Users/apple/Documents/Jecco.jl/data/bubbles/end_data/"

a4_series  = BoundaryTimeSeries(dirname, :a4)
B1_series  = BulkTimeSeries(dirname, :B1, 2)
t, x, y    = get_coords(a4_series,:,:,:)
a4         = a4_series[:,:,:]
B1         = B1_series[1,:,:,:]
_, Nx, Ny  = size(a4)
Nu, Nx, Ny = size(B1)
dx, dy     = (x[2]-x[1], y[2]-y[1])

_, chart2D = get_field(a4_series.ts, it = a4_series.ts.iterations[end],field="a4")
interp   = AdS5_3_1.Linear_Interpolator(chart2D)

Dx = CenteredDiff(1,4,dx,Nx)
Dy = CenteredDiff(1,4,dy,Ny)

x_new = x
y_new = y

Nx_new, Ny_new = (length(x_new), length(y_new))
a4_new_1       = zeros(1,Nx_new, Ny_new)
a4_new_2       = zeros(1,Nx_new, Ny_new)
a4_new_3       = zeros(1,Nx_new, Ny_new)
B1_new         = zeros(Nu, Nx_new, Ny_new)


a4_inter = interp(a4[1,:,:])
@time @fastmath @inbounds @threads for j in 1:Ny_new
        for i in 1:Nx_new
                xx = x_new[i]
                yy = y_new[j]

                a4_new_1[1,i,j] = a4_inter(xx, yy)
        end
end

xm = zeros(Nx_new, Ny_new)
ym = zeros(Nx_new, Ny_new)

@threads for j in 1:Ny_new
        for i in 1:Nx_new
                xm[:,j]  = x_new[:]
                ym[:,j] .= y_new[:]
        end
end



plot(x_new, a4_new_1[1,:,1], lw=3)
plot!(x_new, a4_new_2[1,:,1], lw=3)
plot!(x, a4[1,:,1], lw=3)

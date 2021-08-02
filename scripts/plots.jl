using Jecco, Jecco.AdS5_3_1
using  PyPlot, PyCall
import PyPlot.view_init
#pyplot()

dirname = "/home/mikel/Dropbox/CollisionNewPotential/bubble_expansion/2+1/expansions/data/2D_eA_1.318_eB_ecold"

en         = VEVTimeSeries(dirname, :energy)
t, x, y    = get_coords(en, :, :, :)
Nt, Nx, Ny = size(en)

xx = zeros(Nx, Ny)
yy = zeros(Nx, Ny)

@fastmath @inbounds Threads.@threads for j in 1:Ny
    for i in 1:Nx
        xx[i,j] = x[i]
        yy[i,j] = y[j]
    end
end

xcord = reshape(xx, Nx*Ny)
ycord = reshape(yy, Nx*Ny)
een   = reshape(en[Nt,:,:], Nx*Ny)



figure(figsize=(8,8))
grid(false)
view_init(50,50)
xlabel("xΛ")
ylabel("yΛ")
zlabel("ℰ/Λ⁴")
surf(xcord[1:20:end], ycord[1:20:end], een[1:20:end], cmap=ColorMap(:jet)) 


#=
function plot(a::Real, b::Real)
    xlabel("xΛ")
    ylabel("yΛ")
    zlabel("ℰ/Λ⁴")
    grid(false)
    view_init(a,b)
    surf(xcord[1:20:end], ycord[1:20:end], een[1:20:end], cmap=ColorMap("jet"))
end
=#

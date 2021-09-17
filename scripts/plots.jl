using Jecco, Jecco.AdS5_3_1
using Plots
#using  PyPlot, PyCall
#import PyPlot.view_init
pyplot()

dirname   = "/home/mikel/Documents/Jecco.jl/data/e_1.8_L_20_AH_0.95/"
video_file = "/home/mikel/Documents/Jecco.jl/data/e_1.8_L_20_AH_0.95/e_T2.gif"

en         = VEVTimeSeries(dirname, :energy)
T2         = AdS5_3_1.TTTimeSeries(dirname, :T2)
#hd2        = GWTimeSeries(dirname, :hd2)
t, x, y    = get_coords(T2, :, :, :)
Nt, Nx, Ny = size(T2)



#=
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
emax   = maximum(en[:,:,:])
emin   = minimum(en[:,:,:])
hd2min = minimum(hd2[:,:,:])
hd2max = maximum(hd2[:,:,:])





anim = @animate for n in 1:Nt
    @time begin
        println("progress = $(Int(floor(n/Nt*100)))%")
        een = reshape(en[n,:,:], Nx*Ny)
        p1  = surface(xcord[1:20:end], ycord[1:20:end], een[1:20:end], xlabel="xΛ", ylabel="yΛ", zlabel="ℰ/Λ⁴",color=:jet, zlim=(emin,emax), clims=(emin,emax), camera=(50,50), title="tΛ = $(t[n])")
        nothing
        een  = reshape(hd2[n,:,:], Nx*Ny)
        p2   = surface(xcord[1:20:end], ycord[1:20:end], een[1:20:end], xlabel="xΛ", ylabel="yΛ", zlabel="hd²/Λ⁸",color=:jet, zlim=(hd2min,hd2max), clims=(hd2min,hd2max), camera=(50,50))
        nothing
        een   = reshape(T2[n,:,:], Nx*Ny)
        T2min = minimum(een[:])
        T2max = maximum(een[:])
        p3 = surface(xcord[1:20:end], ycord[1:20:end], een[1:20:end], xlabel="xΛ", ylabel="yΛ", zlabel="T²/Λ⁸",color=:jet, zlim=(T2min,T2max), clims=(T2min,T2max), camera=(50,50))
        nothing
        plot(p1,p2,p3,size=(1500,1500))
        nothing
    end
end
=#


#If x and y are of same length
emax   = maximum(en[:,:,:])
emin   = minimum(en[:,:,:])

anim = @animate for n in 1:Nt
    @time begin
        println("progress = $(Int(floor(n/Nt*100)))%")
        p1  = surface(x, y, en[n,:,:], xlabel="xΛ", ylabel="yΛ", zlabel="ℰ/Λ⁴",color=:jet, zlim=(emin,emax), clims=(emin,emax), camera=(50,50), title="tΛ = $(t[n])")
        nothing
        T2min = minimum(T2[n,:,:])
        T2max = maximum(T2[n,:,:])
        p2 = surface(x, y, T2[n,:,:], xlabel="xΛ", ylabel="yΛ", zlabel="T²/Λ⁸",color=:jet, zlim=(T2min,T2max), clims=(T2min,T2max), camera=(50,50))
        nothing
        plot(p1, p2, size=(2000,2000))
        nothing
    end
end

gif(anim, video_file, fps=1)
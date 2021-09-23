using Jecco, Jecco.AdS5_3_1
using Plots
pyplot()

dirname    = "/Users/apple/Documents/Jecco.jl/data/spinodal_2/"
out_file = "/Users/apple/Dropbox/CollisionNewPotential/GW/spinodal/fourier_modes/"

#=
en         = VEVTimeSeries(dirname, :energy)
T2         = AdS5_3_1.TTTimeSeries(dirname, :T2)
#hd2        = GWTimeSeries(dirname, :hd2)
t, x, y    = get_coords(en, :, :, :)
Nt, Nx, Ny = size(en)
Ntt, _, _ = size(T2)

nt = minimum((Ntt, Nt))



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
#T2min  = minimum(T2[:,:,:])
#T2max  = maximum(T2[:,:,:])

anim = @animate for n in 1:nt
    @time begin
        println("progress = $(n/nt*100)%")
        p1  = surface(x, y, en[n,:,:], xlabel="xΛ", ylabel="yΛ", zlabel="ℰ/Λ⁴",color=:jet, zlim=(emin,emax), clims=(emin,emax), camera=(50,50), title="tΛ = $(t[n])")
        #plot(p1, size=(2000,2000))
        nothing
        T2min = minimum(T2[n,:,:])
        T2max = maximum(T2[n,:,:])
        p2 = surface(x, y, T2[n,:,:], xlabel="xΛ", ylabel="yΛ", zlabel="T²/Λ⁸",color=:jet, zlim=(T2min,T2max), clims=(T2min,T2max), camera=(50,50), title="tΛ = $(t[n])")
        nothing
        plot(p1, p2, size=(2000,1000))
        nothing
    end
end

gif(anim, out_file*"spinodal_a4_e_T2.gif", fps=12)

=#

#Fourier mode plots
nxmax, nymax = (4, 4)
@time a, b, c, d  = AdS5_3_1.Fourier_cos_sin(dirname, :energy)
println("Max cos_cos: $(findmax(abs.(a))[2])")
println("Max cos_sin: $(findmax(abs.(b))[2])")
println("Max sin_cos: $(findmax(abs.(c))[2])")
println("Max sin_sin: $(findmax(abs.(d))[2])")

_, Nkx, Nky = size(a)
a4          = BoundaryTimeSeries(dirname, :a4)
t, _, _     = get_coords(a4, :, 1, 1)
#idx         = findfirst(t .> 100)
#println("Max cos_cos: $(findmax(abs.(a[1:idx,:,:]))[2])")
#println("Max cos_sin: $(findmax(abs.(b[1:idx,:,:]))[2])")
#println("Max sin_cos: $(findmax(abs.(c[1:idx,:,:]))[2])")
#println("Max sin_sin: $(findmax(abs.(d[1:idx,:,:]))[2])")
idx = length(t)
p1 = plot()
p2 = plot()
p3 = plot()
p4 = plot()
for j in 1:nxmax
    for i in 1:nymax
        plot!(p1, t[1:idx], abs.(a[1:idx,i,j]), lw = 1, label="nx=$i ny=$j")
        plot!(p2, t[1:idx], abs.(b[1:idx,i,j]), lw = 1, label="nx=$i ny=$j")
        plot!(p3, t[1:idx], abs.(c[1:idx,i,j]), lw = 1, label="nx=$i ny=$j")
        plot!(p4, t[1:idx], abs.(d[1:idx,i,j]), lw = 1, label="nx=$i ny=$j")
    end
end
xlabel!(p1, "tΛ")
xlabel!(p2, "tΛ")
xlabel!(p3, "tΛ")
xlabel!(p4, "tΛ")
title!(p1, "Cos*Cos")
title!(p2, "Cos*Sin")
title!(p3, "Sin*Cos")
title!(p4, "Sin*Sin")

pfinal = plot(p1, p2, p3, p4, size=(1000, 1000))
nothing
savefig(pfinal, out_file*"spinodal_a4.pdf")
#=
savefig(p1, out_file*"spinodal_2_cos_cos.pdf")
savefig(p2, out_file*"spinodal_2_cos_sin.pdf")
savefig(p3, out_file*"spinodal_2_sin_cos.pdf")
savefig(p4, out_file*"spinodal_2_sin_sin.pdf")
=#

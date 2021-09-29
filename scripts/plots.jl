using Jecco, Jecco.AdS5_3_1
using Plots
pyplot()

dirname    = "/Users/apple/Documents/Jecco.jl/data/new_data/"
out_file   = "/Users/apple/Dropbox/CollisionNewPotential/bubble_expansion/2+1/expansions/eA_1.318_eB_elow/angular_modes/"
name       = "Phi_modes_circ_symm_100.pdf"


en         = VEVTimeSeries(dirname, :energy)
#T2         = AdS5_3_1.TTTimeSeries(dirname, :T2)
#hd2        = GWTimeSeries(dirname, :hd2)
t, x, y    = get_coords(en, :, :, :)
#Nt, Nx, Ny = size(en)
#Nt2, _, _  = size(T2)
#Nt3, _, _  = size(hd2)

#nt = minimum((Nt, Nt2, Nt3))


#=
#ANIMATIONS
#Rectangle boxes
xx = zeros(Nx, Ny)
yy = zeros(Nx, Ny)

@fastmath @inbounds Threads.@threads for j in 1:Ny
    for i in 1:Nx
        xx[i,j] = x[i]
        yy[i,j] = y[j]
    end
end

xcord  = reshape(xx, Nx*Ny)
ycord  = reshape(yy, Nx*Ny)

emax   = maximum(en[:,:,:])
emin   = minimum(en[:,:,:])
hd2min = minimum(hd2[:,:,:])
hd2max = maximum(hd2[:,:,:])

skip   = 20

anim = @animate for n in 1:nt
    @time begin
        println("progress = $(Int(floor(n/nt*100)))%")
        een = reshape(en[n,:,:], Nx*Ny)
        p1  = surface(xcord[1:skip:end], ycord[1:skip:end], een[1:skip:end], xlabel="xΛ", ylabel="yΛ", zlabel="ℰ/Λ⁴",color=:jet, zlim=(emin,emax), clims=(emin,emax), camera=(50,50), title="tΛ = $(t[n])")
        nothing
        een  = reshape(hd2[n,:,:], Nx*Ny)
        p2   = surface(xcord[1:skip:end], ycord[1:skip:end], een[1:skip:end], xlabel="xΛ", ylabel="yΛ", zlabel="hd²/Λ⁸",color=:jet, zlim=(hd2min,hd2max), clims=(hd2min,hd2max), camera=(50,50))
        nothing
        een   = reshape(T2[n,:,:], Nx*Ny)
        T2min = minimum(een[:])
        T2max = maximum(een[:])
        p3 = surface(xcord[1:skip:end], ycord[1:skip:end], een[1:skip:end], xlabel="xΛ", ylabel="yΛ", zlabel="T²/Λ⁸",color=:jet, zlim=(T2min,T2max), clims=(T2min,T2max), camera=(50,50))
        nothing
        plot(p1,p2,p3,size=(3000,2000))
        nothing
    end
end



#=
#SQUARE BOXES
emax   = maximum(en[:,:,:])
emin   = maximum(en[:,:,:])
hd2min = minimum(hd2[:,:,:])
hd2max = maximum(hd2[:,:,:])

anim = @animate for n in 1:nt
    @time begin
        println("progress = $(n/nt*100)%")
        p1  = surface(x, y, en[n,:,:], xlabel="xΛ", ylabel="yΛ", zlabel="ℰ/Λ⁴",color=:jet, zlim=(emin,emax), clims=(emin,emax), camera=(50,50), title="tΛ = $(t[n])")
        nothing
        T2min = minimum(T2[n,:,:])
        T2max = maximum(T2[n,:,:])
        p2 = surface(x, y, T2[n,:,:], xlabel="xΛ", ylabel="yΛ", zlabel="T²/Λ⁸",color=:jet, zlim=(T2min,T2max), clims=(T2min,T2max), camera=(50,50), title="tΛ = $(t[n])")
        nothing
        #hd2min = minimum(hd2[n,:,:])
        #hd2max = maximum(hd2[n,:,:])
        p3 = surface(x, y, hd2[n,:,:], xlabel="xΛ", ylabel="yΛ", zlabel="hd²/Λ⁸",color=:jet, zlim=(hd2min,hd2max), clims=(hd2min,hd2max), camera=(50,50), title="tΛ = $(t[n])")
        nothing
        plot(p1, p2, p3, size=(3000,2000))
        nothing
    end
end
=#

gif(anim, out_file*video_name*".gif", fps=12)
try
    run(`rm $(out_file*video_name*".mp4")`)
catch
end
run(`ffmpeg -i $(out_file*video_name*".gif") -pix_fmt yuv420p $(out_file*video_name*".mp4")`)
run(`rm $(out_file*video_name*".gif")`)
=#



#=
#FOURIER MODES PLOTS
nxmax, nymax = (3, 3)
#=
@time a, b, c, d  = AdS5_3_1.Fourier_cos_sin(dirname, :energy)
println("Max cos_cos: $(findmax(abs.(a))[2])")
println("Max cos_sin: $(findmax(abs.(b))[2])")
println("Max sin_cos: $(findmax(abs.(c))[2])")
println("Max sin_sin: $(findmax(abs.(d))[2])")
_, Nkx, Nky = size(a)
a4          = BoundaryTimeSeries(dirname, :a4)
t, _, _     = get_coords(a4, :, 1, 1)
=#
idx         = findfirst(t .> 100)
println("Max cos_cos: $(findmax(abs.(a[1:idx,:,:]))[2])")
println("Max cos_sin: $(findmax(abs.(b[1:idx,:,:]))[2])")
println("Max sin_cos: $(findmax(abs.(c[1:idx,:,:]))[2])")
println("Max sin_sin: $(findmax(abs.(d[1:idx,:,:]))[2])")
p1 = plot()
p2 = plot()
p3 = plot()
p4 = plot()
legend_pos = :bottomright
for j in 1:nxmax
    for i in 1:nymax
        plot!(p1, t[1:idx], log.(abs.(a[1:idx,i,j])), lw = 1, label="nx=$i ny=$j", legend = legend_pos)
        plot!(p2, t[1:idx], log.(abs.(b[1:idx,i,j])), lw = 1, label="nx=$i ny=$j", legend = legend_pos)
        plot!(p3, t[1:idx], log.(abs.(c[1:idx,i,j])), lw = 1, label="nx=$i ny=$j", legend = legend_pos)
        plot!(p4, t[1:idx], log.(abs.(d[1:idx,i,j])), lw = 1, label="nx=$i ny=$j", legend = legend_pos)
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
savefig(pfinal, out_file)
#=
savefig(p1, out_file*"spinodal_2_cos_cos.pdf")
savefig(p2, out_file*"spinodal_2_cos_sin.pdf")
savefig(p3, out_file*"spinodal_2_sin_cos.pdf")
savefig(p4, out_file*"spinodal_2_sin_sin.pdf")
=#
=#


#FOURIER MODES IN ANGLE

nmax       = 6
legend_pos = :topright

r, (a, b) = AdS5_3_1.Fourier_ϕ(dirname, :energy)

p1 = plot()
p2 = plot()
for n in 2:nmax
    plot!(p1, r, a[:,n], lw=3, label="n=$(n-1)", legend=legend_pos)
    plot!(p2, r, b[:,n], lw=3, label="n=$(n-1)", legend=legend_pos)
end
xlabel!(p1, "rΛ")
xlabel!(p2, "rΛ")
title!(p1, "Cos Modes")
title!(p2, "Sin Modes")

pfinal = plot(p1, p2, size=(1000, 500))
savefig(pfinal, out_file*name)

#=
#SINGLE RECTANGULAR PLOT
xx = zeros(Nx, Ny)
yy = zeros(Nx, Ny)

@inbounds Threads.@threads for j in 1:Ny
    for i in 1:Nx
        xx[i,j] = x[i]
        yy[i,j] = y[j]
    end
end

xcord  = reshape(xx, Nx*Ny)
ycord  = reshape(yy, Nx*Ny)

emax   = maximum(en[:,:,:])
emin   = minimum(en[:,:,:])
T2min  = minimum(T2[:,:,:])
T2max  = maximum(T2[:,:,:])

skip   = 10

een    = reshape(en[1,:,:], Nx*Ny)
p1     = surface(xcord[1:skip:end], ycord[1:skip:end], een[1:skip:end], xlabel="xΛ", ylabel="yΛ", zlabel="ℰ/Λ⁴",color=:jet, zlim=(emin,emax), clims=(emin,emax), camera=(50,50))
nothing
een    = reshape(T2[1,:,:], Nx*Ny)
p2     = surface(xcord[1:skip:end], ycord[1:skip:end], een[1:skip:end], xlabel="xΛ", ylabel="yΛ", zlabel="T2/Λ⁴",color=:jet, zlim=(T2min,T2max), clims=(T2min,T2max), camera=(50,50))
nothing
plot(p1, p2, size=(2000,1000))
=#

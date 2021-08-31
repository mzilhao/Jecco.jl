using Jecco, Jecco.AdS5_3_1
using PyPlot, PyCall
#gr()



main_directory = "/home/mikel/Dropbox/CollisionNewPotential/bubble_expansion/2+1/expansions/data/"

#=
folders        = ["2D_bubble_collision_Lx_120_Ly_100_Nx_300_Ny_250/",
                    "2D_bubble_collision_Lx_150_Ly_100_Nx_375_Ny_250/",
                    "2D_bubble_collision_Lx_200_Ly_140_Nx_500_Ny_350/",
                    "2D_bubble_collision_Lx_210_Ly_160_Nx_525_Ny_400/",
                    "2D_bubble_collision_Lx_210_Ly_180_Nx_525_Ny_450/",
                    "2D_bubble_collision_Lx_230_Ly_180_Nx_575_Ny_450/",
                    "2D_bubble_collision_Lx_230_Ly_200_Nx_575_Ny_500/"
]
=#
#=
folders = ["2D_eA_1.318_eB_ecold_L_100_N_400",
                        "2D_eA_1.318_eB_ecold_L_120_N_480",
                        "2D_eA_1.318_eB_ecold_L_140_N_560",
                        "2D_eA_1.318_eB_ecold_L_152_N_608",
                        "2D_eA_1.318_eB_ecold_L_172_N_688"
]
=#
#directories = main_directory.*folders

#outdir = main_directory.*"2D_eA_1.318_eB_ecold"

#@time AdS5_3_1.same_box_all_runs(outdir, directories)


#=
dirname1 = "/home/mikel/Dropbox/CollisionNewPotential/GW/bubble_collisions/data/2D_bubble_collision"
dirname2 = "/home/mikel/Dropbox/CollisionNewPotential/GW/bubble_collisions/data/2D_bubble_collision_Lx_120_Ly_100_Nx_300_Ny_250"

e1 = BoundaryTimeSeries(dirname1,:a4)
e2 = BoundaryTimeSeries(dirname2,:a4)

_, x1, y1     = get_coords(e1,1,:,:)
Nt1, Nx1, Ny1 = size(e1)
_, x2, y2     = get_coords(e2,1,:,:)
Nt2, Nx2, Ny2 = size(e2)

#println(maximum(abs.(e1[1,:,:]-e2[1,:,:])))

plot(x1, e1[1,:,Int(floor(Ny1/2))], lw=3)
plot!(x2, e2[1,:,Int(floor(Ny2/2))], lw=3)


=#
dirname = "/home/mikel/Dropbox/CollisionNewPotential/GW/bubble_collisions/data/2D_bubble_collision/"
#dirname = "/home/mikel/Documents/Jecco.jl/data/spinodal2D_e_1.0_L_20_N_80_AH_0.95/"
#AdS5_3_1.Energy_to_mathematica(outdir, dit=2)
#AdS5_3_1.GW_to_mathematica(outdir, dit=2)

#e          = VEVTimeSeries(dirname, :energy)
#t, x, y    = get_coords(e, :, :, :)
hxx        = GWTimeSeries(dirname, :hxx)
t, x, y    = get_coords(hxx, :, :, :)
Nt, Nx, Ny = size(hxx)
#Nt, Nx, Ny  = size(e)

xx = zeros(Nx*Ny)
yy = zeros(Nx*Ny)
zz = zeros(Nx*Ny)


xx = zeros(Nx, Ny)
yy = zeros(Nx, Ny)
h  = hxx[end,:,:]

for j in 1:Ny
    for i in 1:Nx
        xx[i,j] = x[i]
        yy[i,j] = y[j]
    end
end

xxx = reshape(xx, Nx*Ny)
yyy = reshape(yy, Nx*Ny)
zzz = reshape(h, Nx*Ny)

plot1 = surf(xxx, yyy, zzz)
#plot2 = plot_surface(xxx, yyy, zzz)
nothing



#plot(x, y, e[end,:,:], st = :surface)
#=
hplot = wireframe(
        x, y, h, transpose=true,
        xlabel = "x",
        ylabel = "y",
        zlabel = "energy",
        reuse  = false,
)
=#

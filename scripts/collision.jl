using Jecco, Jecco.AdS5_3_1
using Plots
gr()


#=
main_directory = "/home/mikel/Dropbox/CollisionNewPotential/GW/bubble_collisions/data/"

folders        = ["2D_bubble_collision_Lx_120_Ly_100_Nx_300_Ny_250/",
                    "2D_bubble_collision_Lx_150_Ly_100_Nx_375_Ny_250/",
                    "2D_bubble_collision_Lx_200_Ly_140_Nx_500_Ny_350/",
                    "2D_bubble_collision_Lx_210_Ly_160_Nx_525_Ny_400/",
                    "2D_bubble_collision_Lx_210_Ly_180_Nx_525_Ny_450/",
                    "2D_bubble_collision_Lx_230_Ly_180_Nx_575_Ny_450/",
                    "2D_bubble_collision_Lx_230_Ly_200_Nx_575_Ny_500/"
]

directories = main_directory.*folders

outdir = main_directory.*"2D_bubble_collision"

@time AdS5_3_1.same_box_all_runs(outdir, directories)




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


outdir = "/home/mikel/Dropbox/CollisionNewPotential/GW/bubble_collisions/data/2D_bubble_collision/"
#AdS5_3_1.Energy_to_mathematica(outdir, dit=2)
AdS5_3_1.GW_to_mathematica(outdir, dit=2)
nothing

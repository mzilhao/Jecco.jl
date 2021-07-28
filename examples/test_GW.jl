using Jecco
using Jecco.AdS5_3_1
using Plots
gr()
include("../src/AdS5_3_1/GW_3.jl")
#dirname = "/home/mikel/Dropbox/CollisionNewPotential/GW/bubble_collisions/data/2D_bubble_collision_Lx_120_Ly_100_Nx_300_Ny_250/"
dirname = "/home/mikel/Dropbox/CollisionNewPotential/GW/bubble_collisions/data/2D_bubble_collision/"
outdir  = "/home/mikel/Documents/Jecco.jl/data/test/"
dt = 0.01
alg = AdS5_3_1.AB3()
@time solve_GW(outdir, dirname, dt=dt, dt_output=0.0, alg=alg, tol=10^-12)
nothing
#=
plot(t_2, real.(h[:,2,2,1]), lw=3)
plot!(t_2, real.(pxk[:,2,2]), lw=3)
xlabel!("t")
=#

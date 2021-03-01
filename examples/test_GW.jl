using Jecco
using Jecco.AdS5_3_1
using Plots
gr()
include("../src/AdS5_3_1/GW_2.jl")
#/Users/apple/Dropbox/CollisionNewPotential/3+1 project/runs/spinodal2D_model1_en0.9_run01/#
dirname = "/Users/apple/Documents/Jecco.jl/data/spinodal2D_e_0.817_Lxy_30_adaptive"
field = :energy
dt = 0.1
alg = AB3()
h, h_t, Tk, kx, ky = solve_GW(dirname, field, dt=dt, alg=alg)


#=
tmax = (length(h[:,1,1])-1)*dt
t = 0:dt:tmax
imax = 1:length(t);
j=15;
plot(t[imax], real.(h[imax,j,1]), label="h", lw=3)
plot!(t[imax], real.(h_t[imax,j,1]), label="h_t", lw=3)
plot!(t[imax], real.(Tk[imax,j,1]), label="Tk", lw=3)
=#

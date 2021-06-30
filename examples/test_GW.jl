using Jecco
using Jecco.AdS5_3_1
using Plots
gr()
include("../src/AdS5_3_1/GW.jl")
#/Users/apple/Dropbox/CollisionNewPotential/3+1 project/runs/spinodal2D_model1_en0.9_run01/#
dirname = "/home/mikel/Documents/Jecco.jl/data/spinodal2D_e_1.0_L_20_N_80_AH_0.95"
dt = 0.08
alg = AB3()
@time h, h_t, t_2, kxx, kyy, pxk, pyk, pzk = solve_GW(dirname, dt=dt, alg=alg, tol=10^-7)
nothing

plot(t_2, real.(h[:,2,2,1]), lw=3)
plot!(t_2, real.(pxk[:,2,2]), lw=3)
xlabel!("t")

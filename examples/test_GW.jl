using Jecco, Jecco.AdS5_3_1


include("/home/mikel/Documents/Jecco.jl/src/AdS5_3_1/GW_2.jl")
dirname = "/home/mikel/Documents/Jecco.jl/data/bubbles/phiM_0.85_phiQ_10/phase_separated/e_1.8_L_20_AH_0.95"
#outdir  = "/home/mikel/Documents/Jecco.jl/data/test/"
outdir  = dirname
#dt = 0.01
#alg = AdS5_3_1.VCABM3()
#alg = AdS5_3_1.AB3()
#adaptive_step = false
#@time solve_GW(outdir, dirname, dt=dt, alg=AdS5_3_1.AB3())
@time get_TT_stress_tensor(outdir, dirname)
nothing
#=
plot(t_2, real.(h[:,2,2,1]), lw=3)
plot!(t_2, real.(pxk[:,2,2]), lw=3)
xlabel!("t")
=#

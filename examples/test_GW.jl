using Jecco, Jecco.AdS5_3_1


dirname = "/Users/apple/Documents/Jecco.jl/data/spinodal_1/"
outdir  = dirname
dt = 0.01
solve_GW(outdir, dirname, dt=dt, alg=AdS5_3_1.AB3())
nothing

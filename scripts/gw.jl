using Jecco, Jecco.AdS5_3_1
using Plots
gr()

dirname = "/home/mikel/Documents/Jecco.jl/data/spinodal_e_1.0_L_20_N_80_AH_1.0/"

dt=0.1
@time h, h_t, t, kx, ky, pxk, pyk, pzk = AdS5_3_1.solve_GW(dirname, dt=dt, alg = AdS5_3_1.AB3())

first_i = findfirst(t .== 0)
m = 1
r = 1:first_i-1
maxk = findfirst(kx .>= 1/dt)
println("Max nx = $maxk")


#=
plot(t[r], real.(h[r,m,1,3]), lw=3)
plot!(t[r], real.(h_t[r,m,1,3]), lw=3)

plot!(t[r], real.(h[r,m,1,4]), lw=3)
plot!(t[r], real.(h_t[r,m,1,4]), lw=3)

plot!(t[r], real.(pyk[r,m,1]-1/3*pxk[r,m,1]-1/3*pyk[r,m,1]-1/3*pzk[r,m,1]), lw=3)
=#

function plots(rr, i, j, k)
    plot(t[rr], real.(h[rr,i,j,k]), lw=3, label = "h")
    plot!(t[rr], real.(h_t[rr,i,j,k]), lw=3, label = "h_t")
    plot!(t[rr], real.(pyk[rr,i,j]-pzk[rr,i,j]), lw=3, label = "source")
    xlabel!("time")
end



plots(r,7,1,3)

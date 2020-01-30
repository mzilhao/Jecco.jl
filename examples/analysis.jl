using Jecco
using Plots
gr()

ts_multi = OpenPMDTimeSeries("./data")

it = 48

phi_1, grid_1 = get_field(ts_multi, it=it, field="phi c=1")
phi_2, grid_2 = get_field(ts_multi, it=it, field="phi c=2")
phi_3, grid_3 = get_field(ts_multi, it=it, field="phi c=3")
phi_4, grid_4 = get_field(ts_multi, it=it, field="phi c=4")

u_1, x, y = grid_1[:]
u_2, x, y = grid_2[:]
u_3, x, y = grid_3[:]
u_4, x, y = grid_4[:]

# heatmap(x, y, phi_1[1,:,:])
heatmap(x, y, phi_2[end,:,:])


# plot!(u_1, phi_1[:,10,10], marker=:circle)
# plot!(u_2, phi_2[:,10,10], marker=:circle)
# plot!(u_3, phi_3[:,10,10], marker=:circle)
# plot!(u_4, phi_4[:,10,10], marker=:circle)

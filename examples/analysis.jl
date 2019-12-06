using Vivi
using Plots
gr()

ts_multi = OpenPMDTimeSeries("./data00_multi")

it = 48

phi_1, info_1 = get_field(ts_multi, it=it, field="phi c=1")
phi_2, info_2 = get_field(ts_multi, it=it, field="phi c=2")
phi_3, info_3 = get_field(ts_multi, it=it, field="phi c=3")
phi_4, info_4 = get_field(ts_multi, it=it, field="phi c=4")

u_1, x, y = info_1.xx
u_2, x, y = info_2.xx
u_3, x, y = info_3.xx
u_4, x, y = info_4.xx

# heatmap(x, y, phi_1[1,:,:])
# heatmap(x, y, phi_2[end,:,:])



ts = OpenPMDTimeSeries("./data00")

it = 48

phi, info = get_field(ts, it=it, field="phi")

u, x, y = info.xx

heatmap(x, y, phi[1,:,:])
# heatmap(x, y, phi_1[1,:,:])

plot(u, phi[:,10,10], marker=:hex)

plot!(u_1, phi_1[:,10,10], marker=:circle)
plot!(u_2, phi_2[:,10,10], marker=:circle)
plot!(u_3, phi_3[:,10,10], marker=:circle)
plot!(u_4, phi_4[:,10,10], marker=:circle)

using Vivi
using Plots
gr()

ts = OpenPMDTimeSeries("./data")

it = 25

phi_1, info_1 = get_field(ts, it=it, field="phi c=1")
phi_2, info_2 = get_field(ts, it=it, field="phi c=2")

u_1, x, y = info_1.xx
u_2, x, y = info_2.xx

heatmap(x, y, phi_1[1,:,:])
# heatmap(x, y, phi_2[end,:,:])

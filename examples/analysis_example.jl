using Jecco, Jecco.AdS5_3_1
using Plots
pyplot()

en_ts    = VEVTimeSeries("spinodal1D_x_model1_en0.9_run02", :energy)
xi_ts    = XiTimeSeries("spinodal1D_x_model1_en0.9_run02")
B1_c2_ts = BulkTimeSeries("spinodal1D_x_model1_en0.9_run02", :B1, 2)

xi = xi_ts[1500, :,3]
t, x, y = get_coords(xi_ts, 1500,:,3)

xi_plot = plot(
    x, xi,
    xlabel = "x",
    ylabel = "ξ",
    title  = "t = $t",
    marker = :hex,
    reuse=false,
)


en = en_ts[1:10:1500,:,3]
t, x, y = get_coords(en_ts, 1:10:1500,:,3)

en_plot = wireframe(
    t, x, en, transpose=true,
    xlabel = "t",
    ylabel = "x",
    title  = "energy",
    reuse=false,
)


B1_c2 = B1_c2_ts[300, :, :, 3]
t, u, x, y = get_coords(B1_c2_ts, 300, :, :, 3)

B1_plot = wireframe(
    u, x, B1_c2, transpose=true,
    xlabel = "u",
    ylabel = "x",
    title  = "B1 c=2",
    reuse=false,
)

display(xi_plot)
display(en_plot)
display(B1_plot)

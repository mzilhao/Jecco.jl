using Vivi
using Plots
using Printf
pyplot()

ts = OpenPMDTimeSeries("./data")

for it in ts.iterations
    phi, info=get_field(ts, it=it, field="phi");
    ti = ts.current_t
    heatmap(phi[1,:,:],
            aspect_ratio=:equal,
            clims=(-0.8, 0.8),
            # title="t =  $ti",
            title="t =  $(@sprintf("%.2f", ti))",
            )

    padit = lpad(it, 4, "0")
    savefig("figs/phi_$padit.png")
end

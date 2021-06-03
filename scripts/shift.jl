using Jecco, Jecco.AdS5_3_1
using HDF5
using Plots


in_dir   = "/home/mikel/Documents/Jecco.jl/data/bubbles/phase_separated"
out_dir  = "/home/mikel/Documents/Jecco.jl/data/bubbles/phase_separated_2"

io = InOut(recover_dir = in_dir, out_dir = out_dir, checkpoint_dir = out_dir,
                 out_boundary_every=1,out_gauge_every=1,out_bulk_every=1,remove_existing = true,)

AdS5_3_1.shift(io)

a4_new = BoundaryTimeSeries(out_dir, :a4)
t_new,x_new,y_new = get_coords(a4_new,:,:,:)

a4_old = BoundaryTimeSeries(in_dir, :a4)
t_old, x_old, y_old = get_coords(a4_old,:,:,:)

plot(x_old, a4_old[end,:,1], lw=3, label="old")
plot!(x_old, a4_new[end,:,31], lw=3, label="new")

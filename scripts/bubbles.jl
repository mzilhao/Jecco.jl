using Jecco
using Jecco.AdS5_3_1
using HDF5
using Plots
gr()

out_dir  = "/home/mikel/Documents/Jecco.jl/data/bubbles/phiM_0.85_phiQ_10/bubble"
A_dir    = "/home/mikel/Documents/Jecco.jl/data/bubbles/phiM_0.85_phiQ_10/state_A_e_1.318"
B_dir    = "/home/mikel/Documents/Jecco.jl/data/bubbles/phiM_0.85_phiQ_10/state_B_e_0.133"
PS_dir   = "/home/mikel/Documents/Jecco.jl/data/bubbles/phiM_0.85_phiQ_10/phase_separated"

grid = SpecCartGrid3D(
    x_min            = -10.0,
    x_max            =  10.0,
    x_nodes          =  80,
    y_min            = -0.5,
    y_max            =  0.5,
    y_nodes          =  5,
    u_outer_min      =  0.1,
    u_outer_max      =  1.005,
    u_outer_domains  =  1,
    u_outer_nodes    =  48,
    u_inner_nodes    =  12,
)

#=
io = InOut(recover_dir = B_dir, out_dir = PS_dir, checkpoint_dir = PS_dir,
                 out_boundary_every=1,out_gauge_every=1,out_bulk_every=1,remove_existing = true,)

AdS5_3_1.shift(io, new_center=(10,0))
=#

io = InOut(recover_dir = out_dir, out_dir = out_dir, checkpoint_dir = out_dir,
                 out_boundary_every=1,out_gauge_every=1,out_bulk_every=1,remove_existing = true,)

AdS5_3_1.bubble_expansion(grid, io, A_dir, B_dir, PS_dir)

#=
a4 = BoundaryTimeSeries(out_dir, :a4)
_,x,y = get_coords(a4,1,:,:)


plot(x, a4[1,:,51], lw=3)
xlabel!("x")
ylabel!("a4")
=#

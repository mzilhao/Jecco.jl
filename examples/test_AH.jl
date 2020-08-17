using Jecco
using Jecco.AdS5_3_1

grid = SpecCartGrid3D(
    x_min            = -5.0,
    x_max            =  5.0,
    #x_nodes          =  32,
    x_nodes          =  12,
    y_min            = -5.0,
    y_max            =  5.0,
    # y_nodes          =  32,
    y_nodes          =  12,
    u_outer_min      =  0.1,
    # u_outer_max      =  2.005,
    u_outer_max      =  1.005,
    # u_outer_domains  =  3,
    u_outer_domains  =  5,
    # u_outer_nodes    =  32,
    u_outer_nodes    =  48,
    u_inner_nodes    =  12,
)

id = BlackBranePert(
    # B1_amp  = 0.005,
    # B1_amp = 0.0,
    # B1_amp = 0.0001,
    B1_nx  = 1,
    B1_ny  = 1,
    #a4_amp = 0.00001,
    a4_amp = 0.01,
    a4_k   = 1,
    xmax   = grid.x_max,
    xmin   = grid.x_min,
    ymax   = grid.y_max,
    ymin   = grid.y_min,
    # AH_pos = 2.0,
    AH_pos = 1.0,
)

evoleq = AffineNull(
    phi0           = 0.0,
    potential      = ZeroPotential(),
    gaugecondition = ConstantAH(u_AH = id.AH_pos),
)

# atlas of grid configuration and respective SystemPartition
atlas     = Atlas(grid)
systems   = SystemPartition(grid)

# allocate variables
boundary  = Boundary(grid)
gauge     = Gauge(grid)
gauge_t   = similar(gauge)

bulkevols      = BulkEvolvedPartition(grid)
bulkconstrains = BulkConstrainedPartition(grid)
bulkderivs     = BulkDerivPartition(grid)

AH_fd_order  = 4

horizoncache   = AdS5_3_1.HorizonCache(systems[end], AH_fd_order)

# initial conditions
init_data!(bulkconstrains, bulkevols, boundary, gauge, systems, evoleq, id)

# guess
uAH = id.AH_pos
# uAH = 0.9
# uAH = 1.9
# uAH = 2.0

sigma = similar(gauge.xi)
fill!(sigma, 1/uAH)

res = 0 * sigma


AdS5_3_1.compute_coeffs_AH!(sigma, gauge, horizoncache, systems[end])

AdS5_3_1.find_AH!(sigma, bulkconstrains[end], bulkevols[end], bulkderivs[end], gauge,
                  horizoncache, systems[end])

Dx = systems[end].Dx
bulk = Bulk(bulkevols[end], bulkconstrains[end])

#This two are the same as for this scenario S = r+xi, so S_AH = sigma+xi. We have to derivate and then evaluate at the horizon as the expressions I use already make a chain rule. Or we can evaluate at AH and then derivate, modifying the expressions we are using.
sigma_x, S_x = similar(sigma[1,:,:]), similar(sigma[1,:,:])
for i in 1:12
    for j in 1:12 
        sigma_x = Dx(sigma, 1,i,j)
        S_x = Dx(horizoncache.bulkhorizon.S_uAH, 1,i,j)
    end
end
sigma_x-S_x

# AdS5_3_1.compute_residual_AH!(res, sigma, gauge, horizoncache, systems[end])

#bla = reshape(horizoncache.bx, (12,12));

#maximum(abs.(bla))

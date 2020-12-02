using Jecco
using Jecco.AdS5_3_1
using FFTW
using Test
include("../scripts/test_utils_2.jl")

#We are going to manipulate an analytic configuration and test the output of create_new_initial. Perturbed conformal BB.
#Miguel and Thanasis: if you change the arguments of the function create_new_initial
#you'll have to change line 75 of this script.
dirname = "/home/mikel/Documentos/Jecco.jl/data/test/"
outdir  = "/home/mikel/Documentos/Jecco.jl/data/test_new/"

grid = SpecCartGrid3D(
    x_min            = -15.0,
    x_max            =  15.0,
    x_nodes          =  32,
    y_min            = -15.0,
    y_max            =  15.0,
    y_nodes          =  32,
    u_outer_min      =  0.1,
    u_outer_max      =  1.005,
    u_outer_domains  =  1,
    u_outer_nodes    =  64,
    u_inner_nodes    =  12,
    fd_order         =  4,
    sigma_diss       =  0.2,
)

new_grid = SpecCartGrid3D(
    x_min            = -30.0,
    x_max            =  30.0,
    x_nodes          =  32,
    y_min            = -30.0,
    y_max            =  30.0,
    y_nodes          =  32,
    u_outer_min      =  0.1,
    u_outer_max      =  1.005,
    u_outer_domains  =  1,
    u_outer_nodes    =  64,
    u_inner_nodes    =  12,
    fd_order         =  4,
    sigma_diss       =  0.2,
)

parameters_collision =AdS5_3_1.new_parameters_coll(
    dirname1  = dirname,
    dirname2  = dirname,
    x1_center = -15.0,
    y1_center = 0.0,
    x2_center = 15.0,
    y2_center = 0.0,
    fx21      = -0.1,
    fy21      = 0.0,
    fx22      = 0.1,
    fy22      = 0.0,
    u_AH      = 1.0,
)


analytic = analytic_parameters(
    energy   = 1.0,
    phi0     = 1.0,
    u_AH     = 1.0,
    a4_amp   = 0.01,
    phi2_amp = 0.01,
    b14_amp  = 0.01,
    b24_amp  = 0.01,
)

tolerance = 5.e-15 #tolerance in the error of our routine.

potential = Phi8Potential(oophiM2=-1.0, oophiQ=0.1,)#It doesn't matter for the purpose of changing a config. it will matter for the gauge change
io = InOut(out_dir = dirname, checkpoint_dir = dirname, out_boundary_every=1,
                  out_gauge_every=1,out_bulk_every=1,remove_existing = true,)

#We generate the data, and the file according to analytic expressions, and modify the configuration analytically
boundary_new, gauge_new, bulkevols_new = test_design_collision(grid,new_grid,analytic,parameters_collision,io,potential)

# We call the function here
io = InOut(recover_dir = dirname, out_dir = outdir, checkpoint_dir = outdir,
               out_boundary_every=1,out_gauge_every=1,out_bulk_every=1,remove_existing = true,)

AdS5_3_1.design_collision(new_grid, io, parameters_collision)


#We compare the analytical result with the routine
Nsys = grid.u_outer_domains+1
for i in 1:Nsys
    B1_aux  = BulkTimeSeries(outdir,:B1,i)
    B2_aux  = BulkTimeSeries(outdir,:B2,i)
    G_aux   = BulkTimeSeries(outdir,:G,i)
    phi_aux = BulkTimeSeries(outdir,:phi,i)

    difB1  = maximum(abs.(B1_aux[end,:,:,:]-bulkevols_new[i].B1))
    difB2  = maximum(abs.(B2_aux[end,:,:,:]-bulkevols_new[i].B2))
    difG   = maximum(abs.(G_aux[end,:,:,:]-bulkevols_new[i].G))
    difphi = maximum(abs.(phi_aux[end,:,:,:]-bulkevols_new[i].phi))

    @testset "Bulk grid $i" begin
        @test difB1  < tolerance
        @test difB2  < tolerance
        @test difG   < tolerance
        @test difphi < tolerance
    end
end

a4_aux  = BoundaryTimeSeries(outdir,:a4)
fx2_aux = BoundaryTimeSeries(outdir,:fx2)
fy2_aux = BoundaryTimeSeries(outdir,:fy2)

difa   = maximum(abs.(a4_aux[end,:,:]-boundary_new.a4[1,:,:]))
diffx2 = maximum(abs.(fx2_aux[end,:,:]-boundary_new.fx2[1,:,:]))
diffy2 = maximum(abs.(fy2_aux[end,:,:]-boundary_new.fy2[1,:,:]))

@testset "Boundary" begin
    @test difa   < tolerance
    @test diffx2 < tolerance
    @test diffy2 < tolerance
end

xi    = XiTimeSeries(outdir)
difxi = maximum(abs.(xi[end,:,:]-gauge_new.xi[1,:,:]))

@testset "Gauge" begin
    @test difxi < tolerance
end
nothing

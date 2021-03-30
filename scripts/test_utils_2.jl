##Functions to test create new initial, design collision and change_gauge. Run the corresponding test files.
Base.@kwdef struct analytic_parameters
    energy   = 1.0
    phi0     = 1.0
    u_AH     = 1.0
    a4_amp   = 0.01
    phi2_amp = 0.01
    b14_amp  = 0.01
    b24_amp  = 0.01
end

a4(x,y,a40,a4_amp,xmid,ymid,Lx,Ly) = a40*(1 + a4_amp*cos(2*pi*(x-xmid)/Lx)*cos(2*pi*(y-ymid)/Ly))
phi2(x,y,phi2_amp,xmid,ymid,Lx,Ly) = phi2_amp*cos(2*pi*(x-xmid)/Lx)*cos(2*pi*(y-ymid)/Ly)
b14(x,y,b14_amp,xmid,ymid,Lx,Ly)   = b14_amp*cos(2*pi*(x-xmid)/Lx)*cos(2*pi*(y-ymid)/Ly)
b24(x,y,b24_amp,xmid,ymid,Lx,Ly)   = b24_amp*cos(2*pi*(x-xmid)/Lx)*cos(2*pi*(y-ymid)/Ly)
B1(u,x,y,b14_amp,xmid,ymid,Lx,Ly)  = u^4*b14(x,y,b14_amp,xmid,ymid,Lx,Ly)
B2(u,x,y,b24_amp,xmid,ymid,Lx,Ly)  = u^4*b24(x,y,b24_amp,xmid,ymid,Lx,Ly)
phi(u,x,y,phi0,xi0,phi2_amp,xmid,ymid,Lx,Ly) = phi0*u-xi0*phi0*u^2+u^3*(phi2(x,y,phi2_amp,xmid,ymid,Lx,Ly)+xi0^2*phi0)

function test_output_analytic(grid::SpecCartGrid3D,prm::analytic_parameters, io::InOut, potential::Potential)
    dirname   = io.out_dir
    energy    = prm.energy
    phi0      = prm.phi0
    u_AH      = prm.u_AH
    a4_amp    = prm.a4_amp
    phi2_amp  = prm.phi2_amp
    b14_amp   = prm.b14_amp
    b24_amp   = prm.b24_amp

    a40 = -4/3*energy
    xi0 = (-a40)^0.25-1/u_AH

    Nx   = grid.x_nodes
    Ny   = grid.y_nodes
    dx   = (grid.x_max-grid.x_min)/Nx
    dy   = (grid.y_max-grid.y_min)/Ny
    Lx   = grid.x_max-grid.x_min
    Ly   = grid.y_max-grid.y_min
    xmid = (grid.x_max+grid.x_min)/2
    ymid = (grid.y_max+grid.y_min)/2

    aa4(x,y)    = a4(x,y,a40,a4_amp,xmid,ymid,Lx,Ly)
    bb14(x,y)   = b14(x,y,b14_amp,xmid,ymid,Lx,Ly)
    bb24(x,y)   = b24(x,y,b24_amp,xmid,ymid,Lx,Ly)
    pphi2(x,y)  = phi2(x,y,phi2_amp,xmid,ymid,Lx,Ly)
    BB1(u,x,y)  = B1(u,x,y,b14_amp,xmid,ymid,Lx,Ly)
    BB2(u,x,y)  = B2(u,x,y,b24_amp,xmid,ymid,Lx,Ly)
    pphi(u,x,y) = phi(u,x,y,phi0,xi0,phi2_amp,xmid,ymid,Lx,Ly)

    atlas     = Atlas(grid)
    systems   = SystemPartition(grid)
    Nsys      = length(systems)
    boundary  = Boundary(grid)
    gauge     = Gauge(grid)
    bulkevols = BulkEvolvedPartition(grid)
    tinfo     = Jecco.TimeInfo(0, 0.0, 0.0, 0.0)
    empty     = Cartesian{1}("u", 0.0, 0.0, 1)
    x         = systems[1].xcoord[:]
    y         = systems[1].ycoord[:]
    chart2D   = Chart(empty, systems[1].xcoord, systems[1].ycoord)

    fill!(boundary.fx2,0.0)
    fill!(boundary.fy2,0.0)
    fill!(gauge.xi,xi0)

    for k in 1:Ny
        for j in 1:Nx
            boundary.a4[1,j,k] = aa4(x[j],y[k])
        end
    end

    for l in 1:Nsys
        u  = systems[l].ucoord[:]
        Nu = length(u)
        for k in 1:Ny
            for j in 1:Nx
                if l == 1
                    bulkevols[1].B1[:,j,k] .= bb14(x[j],y[k])
                    bulkevols[1].B2[:,j,k] .= bb24(x[j],y[k])
                    bulkevols[1].G[:,j,k]  .= 0.0
                    bulkevols[1].phi[:,j,k].= pphi2(x[j],y[k])/phi0^3+xi0^2/phi0^2
                else
                    for i in 1:Nu
                        bulkevols[l].B1[i,j,k]  = BB1(u[i],x[j],y[k])
                        bulkevols[l].B2[i,j,k]  = BB2(u[i],x[j],y[k])
                        bulkevols[l].G[i,j,k]   = 0.0
                        bulkevols[l].phi[i,j,k] = pphi(u[i],x[j],y[k])
                    end
                end
            end
        end
    end

    evolvars = AdS5_3_1.EvolVars(boundary, gauge, bulkevols)
    AdS5_3_1.create_outputs(dirname, evolvars, chart2D, atlas.charts, io, potential, phi0)
end

function test_create_new_data(grid::SpecCartGrid3D,new_grid::SpecCartGrid3D,prm::analytic_parameters,new_parameters::AdS5_3_1.NewParameters,
                                 io::InOut, potential::Potential)

    test_output_analytic(grid,prm,io,potential)
    #We modify the configuration now

    energy    = prm.energy
    phi0      = prm.phi0
    u_AH      = prm.u_AH
    a4_amp    = prm.a4_amp
    phi2_amp  = prm.phi2_amp
    b14_amp   = prm.b14_amp
    b24_amp   = prm.b24_amp

    a40 = -4/3*energy
    xi0 = (-a40)^0.25-1/u_AH

    Nx   = grid.x_nodes
    Ny   = grid.y_nodes
    dx   = (grid.x_max-grid.x_min)/Nx
    dy   = (grid.y_max-grid.y_min)/Ny
    Lx   = grid.x_max-grid.x_min
    Ly   = grid.y_max-grid.y_min
    xmid = (grid.x_max+grid.x_min)/2
    ymid = (grid.y_max+grid.y_min)/2

    aa4(x,y)    = a4(x,y,a40,a4_amp,xmid,ymid,Lx,Ly)
    bb14(x,y)   = b14(x,y,b14_amp,xmid,ymid,Lx,Ly)
    bb24(x,y)   = b24(x,y,b24_amp,xmid,ymid,Lx,Ly)
    pphi2(x,y)  = phi2(x,y,phi2_amp,xmid,ymid,Lx,Ly)
    BB1(u,x,y)  = B1(u,x,y,b14_amp,xmid,ymid,Lx,Ly)
    BB2(u,x,y)  = B2(u,x,y,b24_amp,xmid,ymid,Lx,Ly)
    pphi(u,x,y) = phi(u,x,y,phi0,xi0,phi2_amp,xmid,ymid,Lx,Ly)

    atlas         = Atlas(new_grid)
    systems       = SystemPartition(new_grid)
    Nsys          = length(systems)
    boundary_new  = Boundary(new_grid)
    gauge_new     = Gauge(new_grid)
    bulkevols_new = BulkEvolvedPartition(new_grid)
    empty         = Cartesian{1}("u", 0.0, 0.0, 1)
    x             = systems[1].xcoord[:]
    y             = systems[1].ycoord[:]
    chart2D       = Chart(empty, systems[1].xcoord, systems[1].ycoord)

    Nx = length(x)
    Ny = length(y)

    fill!(boundary_new.fx2,0.0)
    fill!(boundary_new.fy2,0.0)
    fill!(gauge_new.xi,xi0)

    phi2_new = zeros(Nx,Ny)

    for i in 1:Nsys
        u  = systems[i].ucoord[:]
        Nu = length(u)
        for l in 1:Ny
            for k in 1:Nx
                if grid.x_min<=x[k]<=grid.x_max && grid.y_min<=y[l]<=grid.y_max
                    if i == 1
                        boundary_new.a4[1,k,l]      = aa4(x[k],y[l])
                        phi2_new[k,l]               = pphi2(x[k],y[l])
                        bulkevols_new[1].B1[:,k,l] .= bb14(x[k],y[l])
                        bulkevols_new[1].B2[:,k,l] .= bb24(x[k],y[l])
                        bulkevols_new[1].G[:,k,l]  .= 0.0
                        bulkevols_new[1].phi[:,k,l].= pphi2(x[k],y[l])/phi0^3+xi0^2/phi0^2
                    else
                        for j in 1:Nu
                            bulkevols_new[i].B1[j,k,l]  = BB1(u[j],x[k],y[l])
                            bulkevols_new[i].B2[j,k,l]  = BB2(u[j],x[k],y[l])
                            bulkevols_new[i].G[j,k,l]   = 0.0
                            bulkevols_new[i].phi[j,k,l] = pphi(u[j],x[k],y[l])
                        end
                    end
                elseif grid.x_min<=x[k]<=grid.x_max
                    if i == 1
                        boundary_new.a4[1,k,l]      = aa4(x[k],grid.y_min)
                        phi2_new[k,l]               = pphi2(x[k],grid.y_min)
                        bulkevols_new[1].B1[:,k,l] .= bb14(x[k],grid.y_min)
                        bulkevols_new[1].B2[:,k,l] .= bb24(x[k],grid.y_min)
                        bulkevols_new[1].G[:,k,l]  .= 0.0
                        bulkevols_new[1].phi[:,k,l].= pphi2(x[k],grid.y_min)/phi0^3+xi0^2/phi0^2
                    else
                        for j in 1:Nu
                            bulkevols_new[i].B1[j,k,l]  = BB1(u[j],x[k],grid.y_min)
                            bulkevols_new[i].B2[j,k,l]  = BB2(u[j],x[k],grid.y_min)
                            bulkevols_new[i].G[j,k,l]   = 0.0
                            bulkevols_new[i].phi[j,k,l] = pphi(u[j],x[k],grid.y_min)
                        end
                    end
                elseif grid.y_min<=y[l]<=grid.y_max
                    if i == 1
                        boundary_new.a4[1,k,l]      = aa4(grid.x_min,y[l])
                        phi2_new[k,l]               = pphi2(grid.x_min,y[l])
                        bulkevols_new[1].B1[:,k,l] .= bb14(grid.x_min,y[l])
                        bulkevols_new[1].B2[:,k,l] .= bb24(grid.x_min,y[l])
                        bulkevols_new[1].G[:,k,l]  .= 0.0
                        bulkevols_new[1].phi[:,k,l].= pphi2(grid.x_min,y[l])/phi0^3+xi0^2/phi0^2
                    else
                        for j in 1:Nu
                            bulkevols_new[i].B1[j,k,l]  = BB1(u[j],grid.x_min,y[l])
                            bulkevols_new[i].B2[j,k,l]  = BB2(u[j],grid.x_min,y[l])
                            bulkevols_new[i].G[j,k,l]   = 0.0
                            bulkevols_new[i].phi[j,k,l] = pphi(u[j],grid.x_min,y[l])
                        end
                    end
                else
                    if i == 1
                        boundary_new.a4[1,k,l]      = aa4(grid.x_min,grid.y_min)
                        phi2_new[k,l]               = pphi2(grid.x_min,grid.y_min)
                        bulkevols_new[1].B1[:,k,l] .= bb14(grid.x_min,grid.y_min)
                        bulkevols_new[1].B2[:,k,l] .= bb24(grid.x_min,grid.y_min)
                        bulkevols_new[1].G[:,k,l]  .= 0.0
                        bulkevols_new[1].phi[:,k,l].= pphi2(grid.x_min,grid.y_min)/phi0^3+xi0^2/phi0^2
                    else
                        for j in 1:Nu
                            bulkevols_new[i].B1[j,k,l]  = BB1(u[j],grid.x_min,grid.y_min)
                            bulkevols_new[i].B2[j,k,l]  = BB2(u[j],grid.x_min,grid.y_min)
                            bulkevols_new[i].G[j,k,l]   = 0.0
                            bulkevols_new[i].phi[j,k,l] = pphi(u[j],grid.x_min,grid.y_min)
                        end
                    end
                end
            end
        end
    end

    e       = AdS5_3_1.compute_energy.(boundary_new.a4[1,:,:], phi2_new, phi0, potential.oophiM2)
    e_new   = new_parameters.e_new
    fx20    = new_parameters.fx20
    fy20    = new_parameters.fy20
    #change of energy
    if e_new != -1.0
        plan    = plan_rfft(e)
        e0      = real(1/(Nx*Ny)*(plan*e)[1])
        a40_new = real(1/(Nx*Ny)*(plan*boundary_new.a4[1,:,:])[1])
        k       = (-4/3*(e_new-e0)*phi0^4-maximum(boundary_new.a4[1,:,:])+a40_new)/(a40_new-maximum(boundary_new.a4[1,:,:]))
        boundary_new.a4[1,:,:] = k*(boundary_new.a4[1,:,:].-maximum(boundary_new.a4[1,:,:])).+maximum(boundary_new.a4[1,:,:])
        #a4_new  = a4_new.+4/3*(e0-e_new)
    end
    #boosting, with static cold phase (edge).
    if new_parameters.boostx
        boundary_new.fx2[1,:,:] .= fx20.*(e.-minimum(e))
    end
    if new_parameters.boosty
        boundary_new.fy2[1,:,:] .= fy20.*(e.-minimum(e))
    end
    #We add the specified a4 cos perturbation
    ampx    = new_parameters.a4_ampx
    ampy    = new_parameters.a4_ampy
    kx      = new_parameters.a4_kx
    ky      = new_parameters.a4_ky
    xmin    = new_grid.x_min
    xmax    = new_grid.x_max
    ymin    = new_grid.y_min
    ymax    = new_grid.y_max
    xmid    = (xmax+xmin)/2
    ymid    = (ymax+ymin)/2

    plan    = plan_rfft(boundary_new.a4[1,:,:])
    a40_new = real(1/(Nx*Ny)*(plan*boundary_new.a4[1,:,:])[1])

    for j in 1:Ny
        for i in 1:Nx
            boundary_new.a4[1,i,j] -= a40_new*(ampx*cos(2*π*kx*(x[i]-xmid)/(xmax-xmin))
                                                  +ampy*cos(2*π*ky*(y[j]-ymid)/(ymax-ymin)))
        end
    end

    return boundary_new, gauge_new, bulkevols_new
end



function test_design_collision(grid::SpecCartGrid3D, new_grid::SpecCartGrid3D, prm::analytic_parameters,
                                           new_parameters_coll::AdS5_3_1.NewParameters, io::InOut, potential::Potential)

    test_output_analytic(grid,prm,io,potential)

    energy    = prm.energy
    phi0      = prm.phi0
    u_AH      = prm.u_AH
    a4_amp    = prm.a4_amp
    phi2_amp  = prm.phi2_amp
    b14_amp   = prm.b14_amp
    b24_amp   = prm.b24_amp

    a40 = -4/3*energy
    xi0 = (-a40)^0.25-1/u_AH

    Nx   = grid.x_nodes
    Ny   = grid.y_nodes
    dx1  = (grid.x_max-grid.x_min)/Nx
    dy1  = (grid.y_max-grid.y_min)/Ny
    dx2  = dx1
    dy2  = dy1
    Lx   = grid.x_max-grid.x_min
    Ly   = grid.y_max-grid.y_min
    xmid = (grid.x_max+grid.x_min)/2
    ymid = (grid.y_max+grid.y_min)/2
    x1mid = xmid
    y1mid = ymid
    x2mid = xmid
    y2mid = ymid

    aa4(x,y)    = a4(x,y,a40,a4_amp,xmid,ymid,Lx,Ly)
    bb14(x,y)   = b14(x,y,b14_amp,xmid,ymid,Lx,Ly)
    bb24(x,y)   = b24(x,y,b24_amp,xmid,ymid,Lx,Ly)
    pphi2(x,y)  = phi2(x,y,phi2_amp,xmid,ymid,Lx,Ly)
    BB1(u,x,y)  = B1(u,x,y,b14_amp,xmid,ymid,Lx,Ly)
    BB2(u,x,y)  = B2(u,x,y,b24_amp,xmid,ymid,Lx,Ly)
    pphi(u,x,y) = phi(u,x,y,phi0,xi0,phi2_amp,xmid,ymid,Lx,Ly)

    systems = SystemPartition(grid)
    x1      = systems[1].xcoord[:]
    y1      = systems[1].ycoord[:]
    x2      = x1
    y2      = y1

    atlas         = Atlas(new_grid)
    systems       = SystemPartition(new_grid)
    Nsys          = length(systems)
    boundary_new  = Boundary(new_grid)
    gauge_new     = Gauge(new_grid)
    bulkevols_new = BulkEvolvedPartition(new_grid)
    empty         = Cartesian{1}("u", 0.0, 0.0, 1)
    x_new         = systems[1].xcoord[:]
    y_new         = systems[1].ycoord[:]
    chart2D       = Chart(empty, systems[1].xcoord, systems[1].ycoord)

    Nx = length(x_new)
    Ny = length(y_new)

#Position of the center of the boxes inside the new box
    x1mid_new = new_parameters_coll.x1_center
    y1mid_new = new_parameters_coll.y1_center
    x2mid_new = new_parameters_coll.x2_center
    y2mid_new = new_parameters_coll.y2_center
#New bounds of the old boxwes inside the new box
    x1max_new = x1[end]-x1mid+x1mid_new
    y1max_new = y1[end]-y1mid+y1mid_new
    x2max_new = x2[end]-x2mid+x2mid_new
    y2max_new = y2[end]-y2mid+y2mid_new
    x1min_new = x1[1]-x1mid+x1mid_new
    y1min_new = y1[1]-y1mid+y1mid_new
    x2min_new = x2[1]-x2mid+x2mid_new
    y2min_new = y2[1]-y2mid+y2mid_new
    #Check that we are within bounds.
    if x1max_new+dx1>new_grid.x_max || x1min_new<x_new[1] || y1max_new+dy1>new_grid.y_max || y1min_new<y_new[1]
        println("Box 1 out of bounds.")
        return nothing
    elseif x2max_new+dx2>new_grid.x_max || x2min_new<x_new[1] || y2max_new+dy2>new_grid.y_max || y2min_new<y_new[1]
        println("Box 2 out of bounds.")
        return nothing
    end

    e = zeros(length(x1),length(y1))
    for j in 1:length(y1)
        for i in 1:length(x1)
            e[i,j] = AdS5_3_1.compute_energy(aa4(x1[i],y1[j]),pphi2(x1[i],y1[j]),phi0,potential.oophiM2)
        end
    end
    e_min = minimum(e)

    #Read the desired new momenta
    fx21 = new_parameters_coll.fx21
    fy21 = new_parameters_coll.fy21
    fx22 = new_parameters_coll.fx22
    fy22 = new_parameters_coll.fy22
    #Embedding in the new box
    phi2_new = zeros(Nx,Ny)
    for i in 1:Nsys
        u  = systems[i].ucoord[:]
        Nu = length(u)
        for l in 1:Ny
            for k in 1:Nx
                if x1min_new<=x_new[k]<=x1max_new && y1min_new<=y_new[l]<=y1max_new #we displace box 1
                    if i == 1
                        boundary_new.a4[1,k,l]      = aa4(x_new[k]-x1mid_new+x1mid,y_new[l]-y1mid_new+y1mid)
                        phi2_new[k,l]               = pphi2(x_new[k]-x1mid_new+x1mid,y_new[l]-y1mid_new+y1mid)
                        e_aux                       = AdS5_3_1.compute_energy(boundary_new.a4[1,k,l],phi2_new[k,l],phi0,potential.oophiM2)-e_min
                        boundary_new.fx2[1,k,l]     = fx21*e_aux
                        boundary_new.fy2[1,k,l]     = fy21*e_aux
                        bulkevols_new[1].B1[:,k,l] .= bb14(x_new[k]-x1mid_new+x1mid,y_new[l]-y1mid_new+y1mid)
                        bulkevols_new[1].B2[:,k,l] .= bb24(x_new[k]-x1mid_new+x1mid,y_new[l]-y1mid_new+y1mid)
                        bulkevols_new[1].G[:,k,l]  .= 0.0
                        bulkevols_new[1].phi[:,k,l].= pphi2(x_new[k]-x1mid_new+x1mid,y_new[l]-y1mid_new+y1mid)/phi0^3+xi0^2/phi0^2
                    else
                        for j in 1:Nu
                            bulkevols_new[i].B1[j,k,l]  = BB1(u[j],x_new[k]-x1mid_new+x1mid,y_new[l]-y1mid_new+y1mid)
                            bulkevols_new[i].B2[j,k,l]  = BB2(u[j],x_new[k]-x1mid_new+x1mid,y_new[l]-y1mid_new+y1mid)
                            bulkevols_new[i].G[j,k,l]   = 0.0
                            bulkevols_new[i].phi[j,k,l] = pphi(u[j],x_new[k]-x1mid_new+x1mid,y_new[l]-y1mid_new+y1mid)
                        end
                    end
                elseif x2min_new<=x_new[k]<=x2max_new && y2min_new<=y_new[l]<=y2max_new #we displace box 2
                    if i == 1
                        boundary_new.a4[1,k,l]      = aa4(x_new[k]-x2mid_new+x2mid,y_new[l]-y2mid_new+y2mid)
                        phi2_new[k,l]               = pphi2(x_new[k]-x2mid_new+x2mid,y_new[l]-y2mid_new+y2mid)
                        e_aux                       = AdS5_3_1.compute_energy(boundary_new.a4[1,k,l],phi2_new[k,l],phi0,potential.oophiM2)-e_min
                        boundary_new.fx2[1,k,l]     = fx22*e_aux
                        boundary_new.fy2[1,k,l]     = fy22*e_aux
                        bulkevols_new[1].B1[:,k,l] .= bb14(x_new[k]-x2mid_new+x2mid,y_new[l]-y2mid_new+y2mid)
                        bulkevols_new[1].B2[:,k,l] .= bb24(x_new[k]-x2mid_new+x2mid,y_new[l]-y2mid_new+y2mid)
                        bulkevols_new[1].G[:,k,l]  .= 0.0
                        bulkevols_new[1].phi[:,k,l].= pphi2(x_new[k]-x2mid_new+x2mid,y_new[l]-y2mid_new+y2mid)/phi0^3+xi0^2/phi0^2
                    else
                        for j in 1:Nu
                            bulkevols_new[i].B1[j,k,l]  = BB1(u[j],x_new[k]-x2mid_new+x2mid,y_new[l]-y2mid_new+y2mid)
                            bulkevols_new[i].B2[j,k,l]  = BB2(u[j],x_new[k]-x2mid_new+x2mid,y_new[l]-y2mid_new+y2mid)
                            bulkevols_new[i].G[j,k,l]   = 0.0
                            bulkevols_new[i].phi[j,k,l] = pphi(u[j],x_new[k]-x2mid_new+x2mid,y_new[l]-y2mid_new+y2mid)
                        end
                    end
                else
                    if i == 1
                        boundary_new.a4[1,k,l]      = aa4(grid.x_min,grid.y_min)
                        phi2_new[k,l]               = pphi2(grid.x_min,grid.y_min)
                        boundary_new.fx2[1,k,l]     = 0.0
                        boundary_new.fy2[1,k,l]     = 0.0
                        bulkevols_new[1].B1[:,k,l] .= bb14(grid.x_min,grid.y_min)
                        bulkevols_new[1].B2[:,k,l] .= bb24(grid.x_min,grid.y_min)
                        bulkevols_new[1].G[:,k,l]  .= 0.0
                        bulkevols_new[1].phi[:,k,l].= pphi2(grid.x_min,grid.y_min)/phi0^3+xi0^2/phi0^2
                    else
                        for j in 1:Nu
                            bulkevols_new[i].B1[j,k,l]  = BB1(u[j],grid.x_min,grid.y_min)
                            bulkevols_new[i].B2[j,k,l]  = BB2(u[j],grid.x_min,grid.y_min)
                            bulkevols_new[i].G[j,k,l]   = 0.0
                            bulkevols_new[i].phi[j,k,l] = pphi(u[j],grid.x_min,grid.y_min)
                        end
                    end
                end
            end
        end
    end

    fill!(gauge_new.xi,xi0)

    return boundary_new, gauge_new, bulkevols_new
end

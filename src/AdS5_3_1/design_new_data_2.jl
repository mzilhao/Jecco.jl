using FFTW
using HDF5
import Base.Threads.@threads

function change_energy(io::InOut, e_new::T, potential::Potential;
                                        fix::Symbol=:no) where {T<:Real}
    read_dir     = io.recover_dir

    ts                = OpenPMDTimeSeries(read_dir, prefix="boundary_")
    boundary, chart2D = construct_boundary(ts, ts.iterations[end])

    ts       = OpenPMDTimeSeries(read_dir, prefix="gauge_")
    gauge, _ = construct_gauge(ts, ts.iterations[end])

    ts                = OpenPMDTimeSeries(read_dir, prefix="bulk_")
    bulkevols, charts = construct_bulkevols(ts, ts.iterations[end])

    phi0 = try
        ts.params["phi0"]
    catch e
        if isa(e, KeyError)
            0.0   # if "phi0" is not found in the params Dict, set phi0 = 0
            @warn "phi0 not found, setting it to 0.0"
        else
            throw(e)
        end
    end

    oophiM2 = try
        ts.params["oophiM2"]
    catch e
        if isa(e, KeyError)
            0.0   # if "oophiM2" is not found in the params Dict, set oophiM2 = 0
            @warn "oophhiM2 not found, setting it to 0.0"
        else
            throw(e)
        end
    end

    a4        = boundary.a4
    phi2      = BoundaryTimeSeries(read_dir, :phi2)[:,:,:]
    _, Nx, Ny = size(a4)

    e       = AdS5_3_1.compute_energy.(a4, phi2, phi0, oophiM2)
    plan    = plan_rfft(e)
    e0      = real(1/(Nx*Ny)*(plan*e)[1])
    a40     = real(1/(Nx*Ny)*(plan*a4)[1])

    if fix == :low
        k      = (-4/3*(e_new-e0)*phi0^4-maximum(a4)+a40)/(a40-maximum(a4))
        a4    .= k*(a4.-maximum(a4)).+maximum(a4)
    elseif fix == :high
        k      = (-4/3*(e_new-e0)*phi0^4-minimum(a4)+a40)/(a40-minimum(a4))
        a4    .= k*(a4.-minimum(a4)).+minimum(a4)
    elseif fix == :no
        a4    .= a4.+4/3*(e0-e_new)
    else
        @warn "Choose a correct option for fix: :high, :low or :no."
        return
    end

    evolvars = AdS5_3_1.EvolVars(boundary, gauge, bulkevols)
    create_outputs(io.out_dir, evolvars, chart2D, Tuple(charts), io, potential, phi0)

end

function to1plus1(grid::SpecCartGrid3D, boundary::Boundary, gauge::Gauge,
                                    bulkevols::BulkPartition)


    atlas         = Atlas(grid)
    systems       = SystemPartition(grid)
    Nsys          = length(systems)
    boundary_new  = Boundary(grid)
    gauge_new     = Gauge(grid)
    bulkevols_new = BulkEvolvedPartition(grid)
    empty         = Cartesian{1}("u", 0.0, 0.0, 1)
    chart2D       = Chart(empty, systems[1].xcoord, systems[1].ycoord)

    a4_new  = boundary_new.a4
    fx2_new = boundary_new.fx2
    fy2_new = boundary_new.fy2
    xi_new  = gauge_new.xi
    a4      = boundary.a4
    fx2     = boundary.fx2
    fy2     = boundary.fy2
    xi      = gauge.xi

    _, Nx, Ny         = size(a4)
    _, Nx_new, Ny_new = size(a4_new)

    if Nx != Nx_new
        @warn "New grid has to be identical in the x direction"
        return
    end

    y_index = Int(floor(Ny/2))

    #For when you want to freely specify the y coordinate of the slice.
    # y_index = findfirst(ycoord .>= y_specified)

    a4_new[1,:,:] .= view(a4,1,:,y_index)
    fx2_new[1,:,:].= view(fx2,1,:,y_index)
    xi_new[1,:,:] .= view(xi,1,:,y_index)

    fill!(fy2_new, 0.0)

    for n in 1:Nsys

        B1      = bulkevols[n].B1
        B2      = bulkevols[n].B2
        G       = bulkevols[n].G
        phi     = bulkevols[n].phi
        B1_new  = bulkevols_new[n].B1
        B2_new  = bulkevols_new[n].B2
        G_new   = bulkevols_new[n].G
        phi_new = bulkevols_new[n].phi

        B1_new  .= view(B1,:,:,y_index)
        B2_new  .= view(B2,:,:,y_index)
        G_new   .= view(G,:,:,y_index)
        phi_new .= view(phi,:,:,y_index)
        #fill!(G_new, 0.0)
    end

    boundary_new, gauge_new, bulkevols_new, chart2D, atlas

end

function to1plus1(grid::SpecCartGrid3D, io::InOut, potential::Potential)


    read_dir    = io.recover_dir

    ts           = OpenPMDTimeSeries(read_dir, prefix="boundary_")
    boundary, _  = construct_boundary(ts, ts.iterations[end])

    ts           = OpenPMDTimeSeries(read_dir, prefix="gauge_")
    gauge, _     = construct_gauge(ts, ts.iterations[end])

    ts           = OpenPMDTimeSeries(read_dir, prefix="bulk_")
    bulkevols, _ = construct_bulkevols(ts, ts.iterations[end])

    boundary, gauge, bulkevols, chart2D, atlas = to1plus1(grid, boundary, gauge, bulkevols)

    phi0 = try
        ts.params["phi0"]
    catch e
        if isa(e, KeyError)
            0.0   # if "phi0" is not found in the params Dict, set phi0 = 0
            @warn "phi0 not found, setting it to 0.0"
        else
            throw(e)
        end
    end

    oophiM2 = try
        ts.params["oophiM2"]
    catch e
        if isa(e, KeyError)
            0.0   # if "oophiM2" is not found in the params Dict, set oophiM2 = 0
            @warn "oophhiM2 not found, setting it to 0.0"
        else
            throw(e)
        end
    end

    evolvars = AdS5_3_1.EvolVars(boundary, gauge, bulkevols)
    create_outputs(io.out_dir, evolvars, chart2D, atlas.x, io, potential, phi0)
end

#Function to create a circularly symmetric configuration from a y-independent data.
#We are going to use the 2D linear interpolator, if you want to optimize create a 1D interpolator, as you
#only need the x direction
function to2plus1(boundary::Boundary, gauge::Gauge, bulkevols::BulkPartition, chart2D::Chart, charts::Array{Chart, 1})
    a4        = boundary.a4
    xi        = gauge.xi
    _, Nx, Ny = size(a4)
    Nsys      = length(bulkevols)
    xcoord    = chart2D.coords[2]
    ycoord    = CartesianCoord{3,typeof(xcoord[1])}("y", xcoord.min, xcoord.max, xcoord.nodes)
    x         = xcoord[:]
    y         = ycoord[:]
    interp    = Linear_Interpolator(chart2D)
    a4_inter  = interp(a4[1,:,:])
    xi_inter  = interp(xi[1,:,:])
    a4_new    = zeros(1,Nx,Nx)
    fx2_new   = zeros(1,Nx,Nx)
    fy2_new   = zeros(1,Nx,Nx)
    xi_new    = zeros(1,Nx,Nx)

    @fastmath @inbounds @threads for j in 1:Nx
        for i in 1:Nx
            r = sqrt(x[i]^2+y[j]^2)
            if r > x[end] r = x[1] end
            a4_new[1,i,j] = a4_inter(r, 0.)
            xi_new[1,i,j] = xi_inter(r, 0.)
        end
    end

    boundary_new   = Boundary{typeof(a4_new[1,1,1])}(a4_new, fx2_new, fy2_new)
    gauge_new      = Gauge{typeof(xi_new[1,1,1])}(xi_new)
    bulk_new       = Array{BulkEvolved,1}(undef, Nsys)
    charts_new     = Array{Chart, 1}(undef, Nsys)

    for n in 1:Nsys
        B1            = bulkevols[n].B1
        B2            = bulkevols[n].B2
        phi           = bulkevols[n].phi
        Nu, _, _      = size(B1)
        ucoord        = charts[n].coords[1]
        charts_new[n] = Chart(ucoord, xcoord, ycoord)
        B1_new        = zeros(Nu, Nx, Nx)
        B2_new        = zeros(Nu, Nx, Nx)
        G_new         = zeros(Nu, Nx, Nx)
        phi_new       = zeros(Nu, Nx, Nx)

        for i in 1:Nu
            B1_inter  = interp(B1[i,:,:])
            B2_inter  = interp(B2[i,:,:])
            phi_inter = interp(phi[i,:,:])
            @fastmath @inbounds @threads for k in 1:Nx
                for j in 1:Nx
                    r = sqrt(x[j]^2+y[k]^2)
                    if r > x[end] r = x[1] end
                    B1_new[i,j,k]   = B1_inter(r, 0.)
                    B2_new[i,j,k]   = B2_inter(r, 0.)
                    phi_new[i,j,k]  = phi_inter(r, 0.)
                end
            end
        end
        bulk_new[n] = BulkEvolved{typeof(B1_new[1,1,1])}(B1_new, B2_new, G_new, phi_new)
    end
    bulkevols_new   = AdS5_3_1.BulkPartition((bulk_new...))
    empty           = Cartesian{1}("u", 0.0, 0.0, 1)
    chart2D         = Chart(empty, xcoord, ycoord)
    evolvars        = AdS5_3_1.EvolVars(boundary_new, gauge_new, bulkevols_new)

    evolvars, chart2D, charts_new
end

function to2plus1(io::InOut, potential::Potential)
    read_dir          = io.recover_dir
    ts                = OpenPMDTimeSeries(read_dir, prefix="boundary_")
    boundary, chart2D = construct_boundary(ts, ts.iterations[end])
    ts                = OpenPMDTimeSeries(read_dir, prefix="gauge_")
    gauge, _          = construct_gauge(ts, ts.iterations[end])
    ts                = OpenPMDTimeSeries(read_dir, prefix="bulk_")
    bulkevols, charts = construct_bulkevols(ts, ts.iterations[end])

    phi0 = try
        ts.params["phi0"]
    catch e
        if isa(e, KeyError)
            0.0   # if "phi0" is not found in the params Dict, set phi0 = 0
            @warn "phi0 not found, setting it to 0.0"
        else
            throw(e)
        end
    end

    oophiM2 = try
        ts.params["oophiM2"]
    catch e
        if isa(e, KeyError)
            0.0   # if "oophiM2" is not found in the params Dict, set oophiM2 = 0
            @warn "oophhiM2 not found, setting it to 0.0"
        else
            throw(e)
        end
    end

    evolvars, chart2D, charts = to2plus1(boundary, gauge, bulkevols, chart2D, charts)
    create_outputs(io.out_dir, evolvars, chart2D, Tuple(charts), io, potential, phi0)
end

function cut_1D_hole(boundary::Boundary, R::T, chart2D::Chart, chart2D_new::Chart) where {T<:Real}
    a4        = boundary.a4
    fx2       = boundary.fx2
    fy2       = boundary.fy2
    interp    = Linear_Interpolator(chart2D)
    a4_inter  = interp(a4[1,:,:])
    fx2_inter = interp(fx2[1,:,:])
    fy2_inter = interp(fy2[1,:,:])
    x_old     = chart2D.coords[2][:]
    y_old     = chart2D.coords[3][:]
    x_new     = chart2D_new.coords[2][:]
    y_new     = chart2D_new.coords[3][:]
    Nx        = chart2D_new.coords[2].nodes
    Ny        = chart2D_new.coords[3].nodes
    a4_new    = zeros(1, Nx, Ny)
    fx2_new   = zeros(1, Nx, Ny)
    fy2_new   = zeros(1, Nx, Ny)

    @fastmath @inbounds @threads for j in 1:Ny
        for i in 1:Nx
            x              = x_new[i]+sign(x_new[i])*R
            y              = y_new[j]
            a4_new[1,i,j]  = a4_inter(x, y)
            fx2_new[1,i,j] = fx2_inter(x, y)
            fy2_new[1,i,j] = fy2_inter(x, y)
        end
    end
    Boundary{typeof(a4_new[1,1,1])}(a4_new, fx2_new, fy2_new)
end

function cut_1D_hole(gauge::Gauge, R::T, chart2D::Chart, chart2D_new::Chart) where {T<:Real}
    xi        = gauge.xi
    interp    = Linear_Interpolator(chart2D)
    xi_inter  = interp(xi[1,:,:])
    x_old     = chart2D.coords[2][:]
    y_old     = chart2D.coords[3][:]
    x_new     = chart2D_new.coords[2][:]
    y_new     = chart2D_new.coords[3][:]
    Nx        = chart2D_new.coords[2].nodes
    Ny        = chart2D_new.coords[3].nodes
    xi_new    = zeros(1, Nx, Ny)

    @fastmath @inbounds @threads for j in 1:Ny
        for i in 1:Nx
            x              = x_new[i]+sign(x_new[i])*R
            y              = y_new[j]
            xi_new[1,i,j]  = xi_inter(x, y)
        end
    end
    Gauge{typeof(xi_new[1,1,1])}(xi_new)
end

function cut_1D_hole(bulkevols::BulkPartition, R::T, chart2D::Chart, chart2D_new::Chart) where {T<:Real}
    Nsys      = length(bulkevols)
    interp    = Linear_Interpolator(chart2D)
    x_old     = chart2D.coords[2][:]
    y_old     = chart2D.coords[3][:]
    x_new     = chart2D_new.coords[2][:]
    y_new     = chart2D_new.coords[3][:]
    Nx        = chart2D_new.coords[2].nodes
    Ny        = chart2D_new.coords[3].nodes
    bulk      = Array{BulkEvolved, 1}(undef, Nsys)

    for n in 1:Nsys
        B1        = bulkevols[n].B1
        B2        = bulkevols[n].B2
        G         = bulkevols[n].G
        phi       = bulkevols[n].phi
        Nu, _, _  = size(B1)
        B1_new    = zeros(Nu, Nx, Ny)
        B2_new    = zeros(Nu, Nx, Ny)
        G_new     = zeros(Nu, Nx, Ny)
        phi_new   = zeros(Nu, Nx, Ny)
        for i in 1:Nu
            B1_inter  = interp(B1[i,:,:])
            B2_inter  = interp(B2[i,:,:])
            G_inter   = interp(G[i,:,:])
            phi_inter = interp(phi[i,:,:])
            @fastmath @inbounds @threads for k in 1:Ny
                for j in 1:Nx
                    x  = x_new[j]+sign(x_new[j])*R
                    y  = y_new[k]
                    B1_new[i,j,k]  = B1_inter(x, y)
                    B2_new[i,j,k]  = B2_inter(x, y)
                    G_new[i,j,k]   = G_inter(x, y)
                    phi_new[i,j,k] = phi_inter(x, y)
                end
            end
        end
        bulk[n] = BulkEvolved(B1_new, B2_new, G_new, phi_new)
    end

    AdS5_3_1.BulkPartition((bulk...))
end

function cut_1D_hole(io::InOut, potential::Potential, R::T, Nx::Int) where {T<:Real}
    read_dir     = io.recover_dir

    ts                = OpenPMDTimeSeries(read_dir, prefix="boundary_")
    it_boundary       = ts.iterations[end]
    boundary, chart2D = construct_boundary(ts, it_boundary)
    t_boundary        = ts.current_t

    ts       = OpenPMDTimeSeries(read_dir, prefix="gauge_")
    it_gauge = ts.iterations[end]
    gauge, _ = construct_gauge(ts, it_gauge)
    t_gauge  = ts.current_t

    ts                = OpenPMDTimeSeries(read_dir, prefix="bulk_")
    it_bulk           = ts.iterations[end]
    bulkevols, charts = construct_bulkevols(ts, it_bulk)
    t_bulk            = ts.current_t

    it    = max(it_bulk, it_gauge, it_bulk)
    t     = max(t_bulk, t_gauge, t_bulk)
    tinfo = Jecco.TimeInfo(it, t, 0.0, 0.0)

    empty       = Cartesian{1}("u", 0.0, 0.0, 1)
    xcoord      = Cartesian{2}("x", chart2D.coords[2].min+R, chart2D.coords[2].max-R, Nx)
    ycoord      = chart2D.coords[3]
    chart2D_new = Chart(empty, xcoord, ycoord)
    Nsys        = length(bulkevols)
    charts_new  = Array{Chart, 1}(undef, Nsys)
    for n in 1:Nsys
        charts_new[n] = Chart(charts[n].coords[1], xcoord, ycoord)
    end

    @time boundary_new  = cut_1D_hole(boundary, R, chart2D, chart2D_new)
    @time gauge_new     = cut_1D_hole(gauge, R, chart2D, chart2D_new)
    @time bulkevols_new = cut_1D_hole(bulkevols, R, chart2D, chart2D_new)

    phi0 = try
        ts.params["phi0"]
    catch e
        if isa(e, KeyError)
            0.0   # if "phi0" is not found in the params Dict, set phi0 = 0
            @warn "phi0 not found, setting it to 0.0"
        else
            throw(e)
        end
    end

    oophiM2 = try
        ts.params["oophiM2"]
    catch e
        if isa(e, KeyError)
            0.0   # if "oophiM2" is not found in the params Dict, set oophiM2 = 0
            @warn "oophhiM2 not found, setting it to 0.0"
        else
            throw(e)
        end
    end
    evolvars = AdS5_3_1.EvolVars(boundary_new, gauge_new, bulkevols_new)
    create_outputs(io.out_dir, evolvars, chart2D_new, Tuple(charts_new), io, potential, phi0,
                                    t=t, it=it)
end

function cut_circular_hole(boundary::Boundary, R::T, chart2D::Chart, chart2D_new::Chart) where {T<:Real}
    a4        = boundary.a4
    fx2       = boundary.fx2
    fy2       = boundary.fy2
    interp    = Linear_Interpolator(chart2D)
    a4_inter  = interp(a4[1,:,:])
    fx2_inter = interp(fx2[1,:,:])
    fy2_inter = interp(fy2[1,:,:])
    x_old     = chart2D.coords[2][:]
    y_old     = chart2D.coords[3][:]
    x_new     = chart2D_new.coords[2][:]
    y_new     = chart2D_new.coords[3][:]
    Nx        = chart2D_new.coords[2].nodes
    Ny        = chart2D_new.coords[3].nodes
    #xcoord    = Cartesian{2}("x", x_old[1]+R, x_old[end]-R, Nx)
    #ycoord    = Cartesian{3}("y", y_old[1]+R, y_old[end]-R, Ny)
    #x_new     = xcoord[:]
    #y_new     = ycoord[:]
    a4_new    = zeros(1, Nx, Ny)
    fx2_new   = zeros(1, Nx, Ny)
    fy2_new   = zeros(1, Nx, Ny)

    @fastmath @inbounds @threads for j in 1:Ny
        for i in 1:Nx
            x  = x_new[i]
            y  = y_new[j]
            r  = sqrt(x^2+y^2)
            if r == 0.0
                xx, yy = (R, 0.)
            else
                xx = (R+r)/r*x
                yy = (R+r)/r*y
            end
            a4_new[1,i,j]  = a4_inter(xx, yy)
            fx2_new[1,i,j] = fx2_inter(xx, yy)
            fy2_new[1,i,j] = fy2_inter(xx, yy)
        end
    end
    Boundary{typeof(a4_new[1,1,1])}(a4_new, fx2_new, fy2_new)
end

function cut_circular_hole(gauge::Gauge, R::T, chart2D::Chart, chart2D_new::Chart) where {T<:Real}
    xi        = gauge.xi
    interp    = Linear_Interpolator(chart2D)
    xi_inter  = interp(xi[1,:,:])
    x_old     = chart2D.coords[2][:]
    y_old     = chart2D.coords[3][:]
    x_new     = chart2D_new.coords[2][:]
    y_new     = chart2D_new.coords[3][:]
    Nx        = chart2D_new.coords[2].nodes
    Ny        = chart2D_new.coords[3].nodes
    #xcoord    = Cartesian{2}("x", x_old[1]+R, x_old[end]-R, Nx)
    #ycoord    = Cartesian{3}("y", y_old[1]+R, y_old[end]-R, Ny)
    #x_new     = xcoord[:]
    #y_new     = ycoord[:]
    xi_new    = zeros(1, Nx, Ny)

    @fastmath @inbounds @threads for j in 1:Ny
        for i in 1:Nx
            x  = x_new[i]
            y  = y_new[j]
            r  = sqrt(x^2+y^2)
            if r == 0.0
                xx, yy = (R, 0.)
            else
                xx = (R+r)/r*x
                yy = (R+r)/r*y
            end
            xi_new[1,i,j]  = xi_inter(xx, yy)
        end
    end
    Gauge{typeof(xi_new[1,1,1])}(xi_new)
end

function cut_circular_hole(bulkevols::BulkPartition, R::T, chart2D::Chart, chart2D_new::Chart) where {T<:Real}
    Nsys      = length(bulkevols)
    interp    = Linear_Interpolator(chart2D)
    x_old     = chart2D.coords[2][:]
    y_old     = chart2D.coords[3][:]
    x_new     = chart2D_new.coords[2][:]
    y_new     = chart2D_new.coords[3][:]
    Nx        = chart2D_new.coords[2].nodes
    Ny        = chart2D_new.coords[3].nodes
    #xcoord    = Cartesian{2}("x", x_old[1]+R, x_old[end]-R, Nx)
    #ycoord    = Cartesian{3}("y", y_old[1]+R, y_old[end]-R, Ny)
    #x_new     = xcoord[:]
    #y_new     = ycoord[:]
    bulk      = Array{BulkEvolved, 1}(undef, Nsys)

    for n in 1:Nsys
        B1        = bulkevols[n].B1
        B2        = bulkevols[n].B2
        G         = bulkevols[n].G
        phi       = bulkevols[n].phi
        Nu, _, _  = size(B1)
        B1_new    = zeros(Nu, Nx, Ny)
        B2_new    = zeros(Nu, Nx, Ny)
        G_new     = zeros(Nu, Nx, Ny)
        phi_new   = zeros(Nu, Nx, Ny)
        for i in 1:Nu
            B1_inter  = interp(B1[i,:,:])
            B2_inter  = interp(B2[i,:,:])
            G_inter   = interp(G[i,:,:])
            phi_inter = interp(phi[i,:,:])
            @fastmath @inbounds @threads for k in 1:Ny
                for j in 1:Nx
                    x  = x_new[j]
                    y  = y_new[k]
                    r  = sqrt(x^2+y^2)
                    if r == 0.0
                        xx, yy = (R, 0.)
                    else
                        xx = (R+r)/r*x
                        yy = (R+r)/r*y
                    end
                    B1_new[i,j,k]  = B1_inter(xx, yy)
                    B2_new[i,j,k]  = B2_inter(xx, yy)
                    G_new[i,j,k]   = G_inter(xx, yy)
                    phi_new[i,j,k] = phi_inter(xx, yy)
                end
            end
        end
        bulk[n] = BulkEvolved(B1_new, B2_new, G_new, phi_new)
    end

    AdS5_3_1.BulkPartition((bulk...))
end

function cut_circular_hole(io::InOut, potential::Potential, R::T, Nx::Int, Ny::Int) where {T<:Real}
    read_dir     = io.recover_dir

    ts                = OpenPMDTimeSeries(read_dir, prefix="boundary_")
    it_boundary       = ts.iterations[end]
    boundary, chart2D = construct_boundary(ts, it_boundary)
    t_boundary        = ts.current_t

    ts       = OpenPMDTimeSeries(read_dir, prefix="gauge_")
    it_gauge = ts.iterations[end]
    gauge, _ = construct_gauge(ts, it_gauge)
    t_gauge  = ts.current_t

    ts                = OpenPMDTimeSeries(read_dir, prefix="bulk_")
    it_bulk           = ts.iterations[end]
    bulkevols, charts = construct_bulkevols(ts, it_bulk)
    t_bulk            = ts.current_t

    it    = max(it_bulk, it_gauge, it_bulk)
    t     = max(t_bulk, t_gauge, t_bulk)
    tinfo = Jecco.TimeInfo(it, t, 0.0, 0.0)

    empty       = Cartesian{1}("u", 0.0, 0.0, 1)
    xcoord      = Cartesian{2}("x", chart2D.coords[2].min+R, chart2D.coords[2].max-R, Nx)
    ycoord      = Cartesian{3}("y", chart2D.coords[3].min+R, chart2D.coords[3].max-R, Ny)
    chart2D_new = Chart(empty, xcoord, ycoord)
    Nsys        = length(bulkevols)
    charts_new  = Array{Chart, 1}(undef, Nsys)
    for n in 1:Nsys
        charts_new[n] = Chart(charts[n].coords[1], xcoord, ycoord)
    end

    @time boundary_new  = cut_circular_hole(boundary, R, chart2D, chart2D_new)
    @time gauge_new     = cut_circular_hole(gauge, R, chart2D, chart2D_new)
    @time bulkevols_new = cut_circular_hole(bulkevols, R, chart2D, chart2D_new)

    phi0 = try
        ts.params["phi0"]
    catch e
        if isa(e, KeyError)
            0.0   # if "phi0" is not found in the params Dict, set phi0 = 0
            @warn "phi0 not found, setting it to 0.0"
        else
            throw(e)
        end
    end

    oophiM2 = try
        ts.params["oophiM2"]
    catch e
        if isa(e, KeyError)
            0.0   # if "oophiM2" is not found in the params Dict, set oophiM2 = 0
            @warn "oophhiM2 not found, setting it to 0.0"
        else
            throw(e)
        end
    end
    evolvars = AdS5_3_1.EvolVars(boundary_new, gauge_new, bulkevols_new)
    create_outputs(io.out_dir, evolvars, chart2D_new, Tuple(charts_new), io, potential, phi0,
                                    t=t, it=it)
end

#For the moment we use the same grid as in the phase separation file and we will use create_new_data to change anything at the end.
#We use the B's and G of the PS state, and we modify the scalar fields and boundary data. f2 to 0 for the moment. Center the low energy
#phase in the middle of the box, so that the metaestable state will lie outside
function bubble_expansion(grid::SpecCartGrid3D, io::InOut, potential::Potential, A_dir::String, B_dir::String, PS_dir::String;
                                    same_spacing::Symbol=:no, b_cold::Bool=false)
    atlas   = Atlas(grid)
    systems = SystemPartition(grid)
    Nsys    = length(systems)

    a4_A  = BoundaryTimeSeries(A_dir,:a4)[end,1,1]
    a4_B  = BoundaryTimeSeries(B_dir,:a4)[end,1,1]
    a4_PS = BoundaryTimeSeries(PS_dir,:a4)
    xi_A  = XiTimeSeries(A_dir)[end,1,1]
    xi_B  = XiTimeSeries(B_dir)[end,1,1]
    xi_PS = XiTimeSeries(PS_dir)
    _, chart2D = get_field(a4_PS.ts, it=0, field="a4")
    _, x, y = chart2D[:]
    Nx    = chart2D.coords[2].nodes
    Ny    = chart2D.coords[3].nodes
    hot   = CartesianIndex(1,1)
    cold  = CartesianIndex(Int(floor(Nx/2)),Int(floor(Ny/2)))
    if b_cold
        a4_B  = a4_PS[end,cold]
        xi_B  = xi_PS[end,cold]
    end
    k_a4  = (a4_A-a4_B)/(a4_PS[end,hot]-a4_PS[end,cold])
    k_xi  = (xi_A-xi_B)/(xi_PS[end,hot]-xi_PS[end,cold])

    a4  = zeros(1,Nx,Ny)
    xi  = zeros(1,Nx,Ny)
    fx2 = zeros(1,Nx,Ny)
    fy2 = zeros(1,Nx,Ny)

    a4[1,:,:] = k_a4*(a4_PS[end,:,:].-a4_PS[end,cold]).+a4_B
    xi[1,:,:] = k_xi*(xi_PS[end,:,:].-xi_PS[end,cold]).+xi_B

    boundary  = Boundary{typeof(a4[1,1,1])}(a4,fx2,fy2)
    gauge     = Gauge{typeof(xi[1,1,1])}(xi)
    bulk      = Array{BulkEvolved,1}(undef, Nsys)

    for n in 1:Nsys
        phi_A  = BulkTimeSeries(A_dir, :phi, n)[end,:,1,1]
        phi_B  = BulkTimeSeries(B_dir, :phi, n)[end,:,1,1]
        phi_PS = BulkTimeSeries(PS_dir, :phi, n)
        if b_cold
            phi_B .= phi_PS[end,:,cold]
        end
        B1     = BulkTimeSeries(PS_dir, :B1, n)[end,:,:,:]
        B2     = BulkTimeSeries(PS_dir, :B2, n)[end,:,:,:]
        G      = BulkTimeSeries(PS_dir, :G, n)[end,:,:,:]
        phi    = similar(phi_PS[end,:,:,:])
        k      = (phi_A-phi_B)./(phi_PS[end,:,hot]-phi_PS[end,:,cold])
        Nu     = length(phi[:,1,1])
        for i in 1:Nu
            phi[i,:,:] = k[i]*(phi_PS[end,i,:,:].-phi_PS[end,i,cold]).+phi_B[i]
        end
        bulk[n] = BulkEvolved{typeof(B1[1,1,1])}(B1, B2, G, phi)
    end

    bulkevols = AdS5_3_1.BulkPartition((bulk...))

    if same_spacing == :yes
        boundary_new, bulkevols_new, gauge_new = same_grid_spacing(grid, boundary, bulkevols, gauge, chart2D)
    else
        boundary_new, bulkevols_new, gauge_new = different_grid_spacing(grid, boundary, bulkevols, gauge, chart2D)
    end

    phi0 = try
        a4_PS.ts.params["phi0"]
    catch e
        if isa(e, KeyError)
            0.0   # if "phi0" is not found in the params Dict, set phi0 = 0
        else
            throw(e)
        end
    end
    empty   = Cartesian{1}("u", 0.0, 0.0, 1)
    chart2D = Chart(empty, systems[1].xcoord, systems[1].ycoord)
    evolvars = AdS5_3_1.EvolVars(boundary_new, gauge_new, bulkevols_new)
    create_outputs(io.out_dir, evolvars, chart2D, atlas.x, io, potential, phi0)
end

function initial_numerical_phi(grid::SpecCartGrid3D, io::InOut, potential::Potential)
    atlas   = Atlas(grid)
    systems = SystemPartition(grid)
    Nsys    = length(systems)

    boundary       = Boundary(grid)
    gauge          = Gauge(grid)
    bulkevols      = BulkEvolvedPartition(grid)
    tinfo          = Jecco.TimeInfo(0, 0.0, 0.0, 0.0)
    empty          = Cartesian{1}("u", 0.0, 0.0, 1)
    chart2D        = Chart(empty, systems[1].xcoord, systems[1].ycoord)
    _, Nx, Ny      = size(chart2D)

    fid    = h5open(string(io.recover_dir,"phi_mathematica.h5"), "r")
    f      = read(fid)

    phiug1 = f["phi c1"]
    phiug2 = f["phi c2"]
    phi0   = f["phi0"]
    a4     = f["a4"]
    xi     = f["xi"]

    fill!(boundary.a4, a4)
    fill!(boundary.fx2, 0.0)
    fill!(boundary.fy2, 0.0)
    fill!(gauge.xi, xi)
    fill!(bulkevols[1].B1, 0.0)
    fill!(bulkevols[1].B2, 0.0)
    fill!(bulkevols[1].G, 0.0)
    fill!(bulkevols[2].B1, 0.0)
    fill!(bulkevols[2].B2, 0.0)
    fill!(bulkevols[2].G, 0.0)

    phi1 = bulkevols[1].phi
    phi2 = bulkevols[2].phi

    @threads for j in 1:Ny
        for i in 1:Nx
            phi1[:,i,j] = phiug1
            phi2[:,i,j] = phiug2
        end
    end

    evolvars = AdS5_3_1.EvolVars(boundary, gauge, bulkevols)
    create_outputs(io.out_dir, evolvars, chart2D, atlas.x, io, potential, phi0)

    close(fid)

end

#We take 2 equilibrium configurations and we join them along some side and boost them at will.
#Both configuration should have the same u coordinate settings, at least until we implement changes in create_new_data()
#Suggestion: put the one with data in the edge closer to the low energy phase as domain 1.
function design_collision(grid::SpecCartGrid3D, io::InOut, new_parameters_coll::NewParameters)
    dirname1  = new_parameters_coll.dirname1
    dirname2  = new_parameters_coll.dirname2

    atlas     = Atlas(grid)
    systems   = SystemPartition(grid)
    Nsys      = length(systems)

    boundary       = Boundary(grid)
    gauge          = Gauge(grid)
    bulkevols      = BulkEvolvedPartition(grid)
    tinfo          = Jecco.TimeInfo(0, 0.0, 0.0, 0.0)
    empty          = Cartesian{1}("u", 0.0, 0.0, 1)
    chart2D        = Chart(empty, systems[1].xcoord, systems[1].ycoord)
    _,x_new,y_new  = chart2D[:]
    Nx_new         = length(x_new)
    Ny_new         = length(y_new)
    dx_new         = x_new[2]-x_new[1]
    dy_new         = y_new[2]-y_new[1]

    a41   = BoundaryTimeSeries(dirname1, :a4)
    a42   = BoundaryTimeSeries(dirname2, :a4)
    phi21 = BoundaryTimeSeries(dirname1, :phi2)
    phi22 = BoundaryTimeSeries(dirname2, :phi2)
    xi1   = XiTimeSeries(dirname1)
    xi2   = XiTimeSeries(dirname2)

    _, x1, y1 = get_coords(a41, 1, :, :)
    _, x2, y2 = get_coords(a42, 1, :, :)
    #Old runs do not have ts.params field, I set it by hand to what we usually use.
    phi0 = 1.0
    potential = Phi8Potential(oophiM2=-1.0, oophiQ=0.1,)
    try
        phi0  = a4.ts.params["phi0"]
        potential = Phi8Potential(oophiM2=a4.ts.params["oophiM2"], oophiQ=a4.ts.params["oophiQ"],)
    catch
        @warn "No ts.params field, setting phi0=1.0, oophiM2=-1,0 and oophiQ=0.1"
    end

    dx1 = x1[2]-x1[1]
    dy1 = y1[2]-y1[1]
    dx2 = x2[2]-x2[1]
    dy2 = y2[2]-y2[1]

    x1mid = (x1[end]+x1[2])/(x1[end]+x1[2]-2*x1[1])
    y1mid = (y1[end]+y1[2])/(y1[end]+y1[2]-2*y1[1])
    x2mid = (x2[end]+x2[2])/(x2[end]+x2[2]-2*x2[1])
    y2mid = (y2[end]+y2[2])/(y2[end]+y2[2]-2*y2[1])
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
    if x1max_new+dx1>x_new[end]+dx_new || x1min_new<x_new[1] || y1max_new+dy1>y_new[end]+dy_new || y1min_new<y_new[1]
        println("Box 1 out of bounds.")
        return nothing
    elseif x2max_new+dx2>x_new[end]+dx_new || x2min_new<x_new[1] || y2max_new+dy2>y_new[end]+dy_new || y2min_new<y_new[1]
        println("Box 2 out of bounds.")
        return nothing
    end

    interp1 = xy_interpolator(x1,y1)
    interp2 = xy_interpolator(x2,y2)

#Interpolation in the old grids
    a41_inter   = interp1(a41[:,:,:])
    a42_inter   = interp2(a42[:,:,:])
    phi21_inter = interp1(phi21[:,:,:])
    phi22_inter = interp2(phi22[:,:,:])
    xi1_inter   = interp1(xi1[:,:,:])
    xi2_inter   = interp2(xi2[:,:,:])
#We compute the energy density of the boxes
    e1 = zeros(1,length(x1),length(y1))
    e2 = zeros(1,length(x2),length(y2))

    e1[1,:,:] = compute_energy.(a41[end,:,:], phi21[end,:,:], phi0, potential.oophiM2)
    e2[1,:,:] = compute_energy.(a42[end,:,:], phi22[end,:,:], phi0, potential.oophiM2)

    e1_inter = interp1(e1)
    e2_inter = interp2(e2)
    e1_min   = minimum(e1)
    e2_min   = minimum(e2)
    #Read the desired new momenta
    fx21 = new_parameters_coll.fx21
    fy21 = new_parameters_coll.fy21
    fx22 = new_parameters_coll.fx22
    fy22 = new_parameters_coll.fy22
    #We give an initial boost that is below 0.1, and we increment step by step until the specified one.
    f2max = maximum(abs.([fx21,fy21,fx22,fy22]))
    if f2max <= 0.1
        kappa = 1.0
    else
        kappa = 0.1/f2max
    end
    fx2 = zeros(1,Nx_new,Ny_new)
    fy2 = zeros(1,Nx_new,Ny_new)
    boundary.fx2 .= 0.0
    boundary.fy2 .= 0.0

    u = Array{Array{Real,1},1}(undef,Nsys)
#Embedding in the new box
    for i in 1:Nsys
        B11  = BulkTimeSeries(dirname1, :B1, i)
        B12  = BulkTimeSeries(dirname2, :B1, i)
        B21  = BulkTimeSeries(dirname1, :B2, i)
        B22  = BulkTimeSeries(dirname2, :B2, i)
        G1   = BulkTimeSeries(dirname1, :G, i)
        G2   = BulkTimeSeries(dirname2, :G, i)
        phi1 = BulkTimeSeries(dirname1, :phi, i)
        phi2 = BulkTimeSeries(dirname2, :phi, i)

        _,u[i],_,_ = get_coords(B11,1,:,1,1)

        B11_inter  = interp1(B11[end,:,:,:])
        B12_inter  = interp2(B12[end,:,:,:])
        B21_inter  = interp1(B21[end,:,:,:])
        B22_inter  = interp2(B22[end,:,:,:])
        G1_inter   = interp1(G1[end,:,:,:])
        G2_inter   = interp2(G2[end,:,:,:])
        phi1_inter = interp1(phi1[end,:,:,:])
        phi2_inter = interp2(phi2[end,:,:,:])
        for l in 1:Ny_new
            for k in 1:Nx_new
                if x1min_new<=x_new[k]<=x1max_new && y1min_new<=y_new[l]<=y1max_new #we displace box 1
                    bulkevols[i].B1[:,k,l]  = B11_inter(x_new[k]-x1mid_new+x1mid,y_new[l]-y1mid_new+y1mid)
                    bulkevols[i].B2[:,k,l]  = B21_inter(x_new[k]-x1mid_new+x1mid,y_new[l]-y1mid_new+y1mid)
                    bulkevols[i].G[:,k,l]   = G1_inter(x_new[k]-x1mid_new+x1mid,y_new[l]-y1mid_new+y1mid)
                    bulkevols[i].phi[:,k,l] = phi1_inter(x_new[k]-x1mid_new+x1mid,y_new[l]-y1mid_new+y1mid)
                    if i == 1
                        e_aux               = e1_inter(x_new[k]-x1mid_new+x1mid,y_new[l]-y1mid_new+y1mid)[1]-e1_min
                        boundary.a4[1,k,l]  = a41_inter(x_new[k]-x1mid_new+x1mid,y_new[l]-y1mid_new+y1mid)[1]
                        fx2[1,k,l]          = fx21*e_aux
                        fy2[1,k,l]          = fy21*e_aux
                        gauge.xi[1,k,l]     = xi1_inter(x_new[k]-x1mid_new+x1mid,y_new[l]-y1mid_new+y1mid)[1]
                    end
                elseif x2min_new<=x_new[k]<=x2max_new && y2min_new<=y_new[l]<=y2max_new #we displace box 2
                    bulkevols[i].B1[:,k,l]  = B12_inter(x_new[k]-x2mid_new+x2mid,y_new[l]-y2mid_new+y2mid)
                    bulkevols[i].B2[:,k,l]  = B22_inter(x_new[k]-x2mid_new+x2mid,y_new[l]-y2mid_new+y2mid)
                    bulkevols[i].G[:,k,l]   = G2_inter(x_new[k]-x2mid_new+x2mid,y_new[l]-y2mid_new+y2mid)
                    bulkevols[i].phi[:,k,l] = phi2_inter(x_new[k]-x2mid_new+x2mid,y_new[l]-y2mid_new+y2mid)
                    if i == 1
                        e_aux               = e2_inter(x_new[k]-x2mid_new+x2mid,y_new[l]-y2mid_new+y2mid)[1]-e2_min
                        boundary.a4[1,k,l]  = a42_inter(x_new[k]-x2mid_new+x2mid,y_new[l]-y2mid_new+y2mid)[1]
                        fx2[1,k,l]          = fx22*e_aux
                        fy2[1,k,l]          = fy22*e_aux
                        gauge.xi[1,k,l]     = xi2_inter(x_new[k]-x2mid_new+x2mid,y_new[l]-y2mid_new+y2mid)[1]
                    end
                else #we extrapolate the low energy phase of box 1
                    bulkevols[i].B1[:,k,l]  = B11[end,:,1,1]
                    bulkevols[i].B2[:,k,l]  = B21[end,:,1,1]
                    bulkevols[i].G[:,k,l]   = G1[end,:,1,1]
                    bulkevols[i].phi[:,k,l] = phi1[end,:,1,1]
                    if i == 1
                        e_aux               = e1_min
                        boundary.a4[1,k,l]  = a41[1,1,1]
                        fx2[1,k,l]          = 0.0
                        fy2[1,k,l]          = 0.0
                        gauge.xi[1,k,l]     = xi1[1,1,1]
                    end
                end
            end
        end
    end

    #These lines for the test and comment the change gauge
    boundary.fx2 .= fx2
    boundary.fy2 .= fy2
    ############################


#TODO: Change gauge using an independent function and sigma as argument to create the new initial data
#with uAH = 1.0. Think on extending this to small boosts up to a final big one in an automatic way.
#=
if f2max != 0
    u_AH = new_parameters_coll.u_AH
    evoleq = AffineNull(phi0=phi0, potential=potential, gaugecondition = ConstantAH(u_AH = 1.0),)
    sigma  = 0.5*ones(1,Nx_new,Ny_new)
    change_gauge!(sigma, grid, boundary, gauge, bulkevols, evoleq, systems, u_AH)
    epsilon = 0
    evoleq  = AffineNull(phi0=phi0, potential=potential, gaugecondition = ConstantAH(u_AH = u_AH),)
    println("epsilon = $epsilon")
    while epsilon < 1.0
        n     = 0
        #kappa = kappa/0.8
        sigma = 1/(2.0*u[end][end])*ones(1,Nx_new,Ny_new)
        while maximum(1 ./sigma) >= u[end][end]+0.005
            if n != 0
                kappa = kappa*0.8
            end
            n += 1
            aux           = epsilon+kappa
            boundary.fx2 .= fx2*aux^0.5
            boundary.fy2 .= fy2*aux^0.5
            try
                change_gauge!(sigma, grid, boundary, gauge, bulkevols, evoleq, systems, u_AH)
            catch

            end
        end
        epsilon += kappa
        println("epsilon = $epsilon")
        println("max fx2 = $(maximum(boundary.fx2))")
    end
end
=#
#=
u_AH      = 0.9
sigma     = zeros(1,Nx_new,Ny_new)
change_gauge!(sigma, grid, boundary, gauge, bulkevols, evoleq, systems, u_AH)
epsilon   = 1.0
#gauge.xi .= (4/3*epsilon)^0.25 - 1/u_AH
#a4       = boundary.a4
a4_final = zeros(1,Nx_new,Ny_new)
copyto!(a4_final,boundary.a4)
boundary.fx2 .= fx2
boundary.fy2 .= fy2
evoleq  = AffineNull(phi0=phi0, potential=potential, gaugecondition = ConstantAH(u_AH = u_AH),)

while epsilon >= 0
    println("epsilon = $epsilon")
    boundary.a4 .= -4/3*epsilon.+a4_final
    change_gauge!(sigma, grid, boundary, gauge, bulkevols, evoleq, systems, u_AH)
    epsilon -= 0.02
end

change_gauge!(sigma, grid, boundary, gauge, bulkevols, evoleq, systems, 1.0)
=#
#We finally write the output files
evolvars = AdS5_3_1.EvolVars(boundary, gauge, bulkevols)
create_outputs(io.out_dir, evolvars, chart2D, atlas.x, io, potential, phi0)
#return sigma
end

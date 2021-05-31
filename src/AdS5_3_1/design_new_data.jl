using FFTW
using HDF5
import Base.Threads.@spawn
import Base.Threads.@threads
#using Interpolations
#Functions that manipulate the output to generate an new initial state from which to run Jecco.
#It creates a checkpoint file from where to run jecco and also normal files to plot it in mathematica and see
#if it is what you wanted

abstract type NewParameters end

Base.@kwdef struct new_parameters <: NewParameters
    e_new   :: Real = -1.0
    boostx  :: Bool = false
    boosty  :: Bool = false
    fx20    :: Real = 0.0
    fy20    :: Real = 0.0
    a4_ampx :: Real = 0.0
    a4_ampy :: Real = 0.0
    a4_kx   :: Int  = 1
    a4_ky   :: Int  = 1
    u_AH    :: Real = 1.0
end

Base.@kwdef struct new_parameters_coll <: NewParameters
    dirname1      :: String
    dirname2      :: String
    fx21          :: Real = 0.0
    fy21          :: Real = 0.0
    fx22          :: Real = 0.0
    fy22          :: Real = 0.0
    x1_center     :: Real
    y1_center     :: Real
    x2_center     :: Real
    y2_center     :: Real
    u_AH          :: Real = 1.0
end

abstract type XYInterpolator end
 #Interpolator of real functions.
struct xy_interpolator{T<:Real, TP<:Integer} <: XYInterpolator
    xmin      :: T
    xmax      :: T
    ymin      :: T
    ymax      :: T
    Nx        :: TP
    Ny        :: TP
    tolerance :: T
end

function xy_interpolator(x::Array{T,1}, y::Array{T,1}; tolerance::T = 0.0) where T<:Real
    xmin = x[1]
    xmax = x[end]+x[2]-x[1]
    ymin = y[1]
    ymax = y[end]+y[2]-y[1]
    Nx   = length(x)
    Ny   = length(y)

    xy_interpolator{typeof(xmin),typeof(Nx)}(xmin,xmax,ymin,ymax,Nx,Ny,tolerance)
end

#Careful, FFTW decomposes taking as origin the (xmin,ymin), not the middle of the box as I am used to think.
function (interpolator::xy_interpolator)(f::Array{T,3}) where {T<:Real}
    dx   = (interpolator.xmax-interpolator.xmin)/interpolator.Nx
    dy   = (interpolator.ymax-interpolator.ymin)/interpolator.Ny
    xmin = interpolator.xmin
    ymin = interpolator.ymin
    Nu   = length(f[:,1,1])
    plan = plan_fft(f[1,:,:])
    k_x  = 2π*fftfreq(interpolator.Nx,1/dx)
    k_y  = 2π*fftfreq(interpolator.Ny,1/dy)
    fk   = im*zeros(Nu,length(k_x),length(k_y))
    index = Array{Array{CartesianIndex{2},1},1}(undef, Nu)
    tol = interpolator.tolerance
    @fastmath @inbounds @threads for i in 1:Nu
        fk[i,:,:]  = 1/(interpolator.Nx*interpolator.Ny)*(plan*f[i,:,:])
        index[i]   = findall(abs.(fk[i,:,:]) .> tol)
    end
    function (x::TP, y::TP) where {TP<:Real}
        #@assert interpolator.xmin<=x<=interpolator.xmax
        #@assert interpolator.ymin<=y<=interpolator.ymax
        sum     = im*zeros(Nu)
        @fastmath @inbounds @threads for i in 1:Nu
            for I in index[i]
                sum[i] += fk[i,I]*exp(im*(k_x[I[1]]*(x-xmin)+k_y[I[2]]*(y-ymin)))
            end
        end
        real(sum)
    end
    #=
    function (x::Array{TP,2}, y::Array{TP,2}) where {TP<:Real}
        #@assert interpolator.xmin<=x<=interpolator.xmax
        #@assert interpolator.ymin<=y<=interpolator.ymax
        Nx, Ny  = size(x)
        sum     = im*zeros(Nu, Nx, Ny)
        @fastmath @inbounds @threads for i in 1:Nu
            for I in index[i]
                sum[i,:,:] += fk[i,I]*exp.(im*(k_x[I[1]]*(x.-xmin)+k_y[I[2]]*(y.-ymin)))
            end
        end
        real.(sum)
    end
    =#
end

function change_B!(B_old::T, u::T, dxi::T, gidx::Int, grid::Int) where T<:Real
    if grid == 1
        if gidx == 1
            return (1-dxi*u)^-4*B_old
        else
            return u^-4*B_old
        end
    else
        if gidx == 1
            return u^4/(1-dxi*u)^4*B_old
        else
            return B_old
        end
    end
end

function change_phi!(phi_old::T, u::T, xi::T, dxi::T, gidx::Int, grid::Int, phi0::T) where T<:Real
    if phi0 < 1.e-12 #There is no scalar field.
        return 0.0
    end
    if grid == 1
        if gidx == 1
            xi_old = xi+dxi
            return dxi*(dxi-u*dxi^2-2*xi_old+u*xi_old*dxi)/(phi0*(1-dxi*u))^2+(1-dxi*u)^-3*phi_old
        else
            return -1/(phi0*u)^2+xi/(phi0^2*u)+phi0^-3*u^-3*phi_old
        end
    else
        if gidx == 1
            xi_old = xi+dxi
            u_old  = u/(1-dxi*u)
            return phi0*u_old-phi0*xi_old*u_old^2+phi0^3*u_old^3*phi_old
        else
            return phi_old
        end
    end
end

function construct_boundary(ts::OpenPMDTimeSeries, it::Int)
    a4, chart2D = get_field(ts, it=it, field="a4")
    fx2, _      = get_field(ts, it=it, field="fx2")
    fy2, _      = get_field(ts, it=it, field="fy2")
    return Boundary(a4, fx2, fy2), chart2D
end

function construct_gauge(ts::OpenPMDTimeSeries, it::Int)
    xi, chart2D = get_field(ts, it=it, field="xi")
    return Gauge(xi), chart2D
end

function construct_bulkevols(ts::OpenPMDTimeSeries, it::Int)
    Nsys = 0
    c    = 0
    while Nsys == 0
        try
            c += 1
            get_field(ts, it=it, field="B1 c=$c")
        catch
            Nsys = c-1
        end
    end

    charts = Array{Chart,1}(undef, Nsys)
    bulk   = Array{BulkEvolved,1}(undef, Nsys)

    for n in 1:Nsys
        B1, chart  = get_field(ts, it=it, field="B1 c=$n")
        B2, _  = get_field(ts, it=it, field="B2 c=$n")
        G, _   = get_field(ts, it=it, field="G c=$n")
        phi, _ = get_field(ts, it=it, field="phi c=$n")

        charts[n] = chart
        bulk[n]  = BulkEvolved(B1, B2, G, phi)
    end

    return AdS5_3_1.BulkPartition((bulk...)), charts
end

function create_outputs(out_dir, evolvars, chart2D, charts, io, potential, phi0; it=0, t=0)
    tinfo      = Jecco.TimeInfo(it, t, 0.0, 0.0)
    try run(`rm -r $out_dir`) catch end
    run(`mkdir $out_dir`)
    checkpoint = AdS5_3_1.checkpoint_writer(evolvars, chart2D, charts, tinfo, io)
    out        = AdS5_3_1.output_writer(evolvars, chart2D, charts, tinfo, io, potential, phi0)

    checkpoint(evolvars)
    out(evolvars)
    nothing
end

#Creates a checkpoint file with the data of the last iteration in a given folder.
function create_checkpoint(io::InOut, potential::Potential)
    ts          = OpenPMDTimeSeries(io.recover_dir, prefix="boundary_")
    it_boundary = ts.iterations[end]
    boundary, _ = construct_boundary(ts, it_boundary)
    t_boundary  = ts.current_t

    ts             = OpenPMDTimeSeries(io.recover_dir, prefix="gauge_")
    it_gauge       = ts.iterations[end]
    gauge, chart2D = construct_gauge(ts, it_gauge)
    t_gauge        = ts.current_t

    ts                = OpenPMDTimeSeries(io.recover_dir, prefix="bulk_")
    it_bulk           = ts.iterations[end]
    bulkevols, charts = construct_bulkevols(ts, it_bulk)
    t_bulk            = ts.current_t

    phi0 = try
        ts.params["phi0"]
    catch e
        if isa(e, KeyError)
            0.0   # if "phi0" is not found in the params Dict, set phi0 = 0
        else
            throw(e)
        end
    end
    it = max(it_bulk, it_gauge, it_bulk)
    t  = max(t_bulk, t_gauge, t_bulk)

    evolvars = AdS5_3_1.EvolVars(boundary, gauge, bulkevols)
    tinfo    = Jecco.TimeInfo(it, t, 0.0, 0.0)

    try run(`rm -r $(io.out_dir)`) catch end
    run(`mkdir $(io.out_dir)`)
    checkpoint = AdS5_3_1.checkpoint_writer(evolvars, chart2D, charts, tinfo, io)
    checkpoint(evolvars)
end

function change_gauge!(sigma::Array{T,3}, grid::SpecCartGrid3D, boundary::Boundary, gauge::Gauge, bulkevols::BulkPartition,
                               evoleq::AffineNull, sys::SystemPartition, u_AH::T) where T<:Real
    phi0           = evoleq.phi0
    #u_AH           = evoleq.gaugecondition.u_AH
    bulkconstrains = BulkConstrainedPartition(grid)
    bulkderivs     = BulkDerivPartition(grid)
    horizoncache   = HorizonCache(sys[end], evoleq.gaugecondition.fd_order)
    ahf            = AHF()
    nested         = Nested(sys,bulkconstrains)
    nested(bulkevols,boundary,gauge,evoleq)

    xi             = getxi(gauge)
    #We find the AH surface
    fill!(sigma, 1/evoleq.gaugecondition.u_AH)
    find_AH!(sigma,bulkconstrains[end],bulkevols[end],bulkderivs[end],gauge,horizoncache,sys[end],ahf)

    if maximum(1 ./sigma) >= sys[end].ucoord[:][end]+0.005
        return
    end

    dxi = 1/u_AH.-sigma#The diff in the gauge parameter
    #We now change the bulk variables accordingly. As u will possibly map among grids, a possible solution is to build
    #a container for the bulk variables as an array of length Nsys of 3D arrays. Then, to update the bulk we map the u
    #coord to the old one and choose the array belonging to the correct grid. We interpolate and evaluate.
    #As diff grids have different definitions we choose to work in g1 like redefs for interpolations and when getting the new functions
    #we will reconstruct the correct redefinition.
    #change_bulk!(bulkevols,xi,dxi,sys,phi0)
    Nx      = sys[1].xcoord.nodes
    Ny      = sys[1].ycoord.nodes
    Nsys    = length(sys)
    B1      = Array{Array{Real,3},1}(undef,Nsys)
    B2      = Array{Array{Real,3},1}(undef,Nsys)
    G       = Array{Array{Real,3},1}(undef,Nsys)
    phi     = Array{Array{Real,3},1}(undef,Nsys)
    u       = Array{Array{Real,1},1}(undef,Nsys)
    uinterp = Array{ChebInterpolator,1}(undef,Nsys)
    u_max   = zeros(Nsys)

    #We are going to change directly the bulk quantites, so we construct the interpolators and interpolated
    #Functions already.
    warning = 0
    for i in 1:Nsys
        B1[i]      = bulkevols[i].B1
        B2[i]      = bulkevols[i].B2
        G[i]       = bulkevols[i].G
        phi[i]     = bulkevols[i].phi
        u[i]       = sys[i].ucoord[:]
        #uinterp[i] = ChebInterpolator(u[i][1],u[i][end],length(u[i]))
        uinterp[i]    = ChebInterpolator(u[i])
        u_max[i] = u[i][end]
    end
    for k in 1:Ny
        for j in 1:Nx
            #println("We are in thread number $(Threads.threadid())")
            dxxii     = dxi[1,j,k]
            xi[1,j,k]-= dxxii
            xxii      = xi[1,j,k]
            for l in 1:Nsys
                Nu = length(u[l])
                for i in 1:Nu
                    u_new  = u[l][i]
                    u_old  = u_new/(1-dxxii*u_new) #The new point maps to this old point
                    gidx   = findfirst(u_max.-u_old.>=0) #We find to what grid it belongs
                    if gidx == nothing
                        gidx  = Nsys
                        u_old = u_max[end]
                        warning = 1
                    end
                    bulkevols[l].B1[i,j,k]  = change_B!(uinterp[gidx](B1[gidx][:,j,k])(u_old),u_new,dxxii,gidx,l)
                    bulkevols[l].B2[i,j,k]  = change_B!(uinterp[gidx](B2[gidx][:,j,k])(u_old),u_new,dxxii,gidx,l)
                    bulkevols[l].G[i,j,k]   = change_B!(uinterp[gidx](G[gidx][:,j,k])(u_old),u_new,dxxii,gidx,l)
                    bulkevols[l].phi[i,j,k] = change_phi!(uinterp[gidx](phi[gidx][:,j,k])(u_old),u_new,xxii,dxxii,gidx,l,phi0)
                end
            end
        end
    end
    if warning == 1 @warn "We are extrapolating outside the original domain" end
    #fill!(sigma, 1/u_AH)
    #nested = Nested(sys,bulkconstrains)
    #nested(bulkevols,boundary,gauge,evoleq)
    #find_AH!(sigma,bulkconstrains[end],bulkevols[end],bulkderivs[end],gauge,horizoncache,sys[end],ahf)
    #return sigma, bulkevols, bulkconstrains, gauge
end

function same_grid_spacing!(f_new::Array{T,3}, f::Array{T,3}, x_indices::Tuple{Int,Int},
                                            y_indices::Tuple{Int,Int}) where {T<:Real}

    ifirst, ilast     = x_indices
    jfirst, jlast     = y_indices
    _, Nx, Ny         = size(f)
    _, Nx_new, Ny_new = size(f_new)

    if ilast-ifirst+1 != Nx || jlast-jfirst+1 != Ny
        @warn "Something is wrong with the x and y indices"
        return
    end

    f_new[:,1:ifirst,1:jfirst]        .= f[:,1,1]
    @fastmath @inbounds @threads for i in 1:ifirst
        f_new[:,i,jfirst:jlast]       .= @view f[:,1,:]
    end
    f_new[:,1:ifirst,jlast:end]       .= f[:,1,end]

    @fastmath @inbounds @threads for i in ifirst:ilast
        n = i-ifirst+1
        f_new[:,i,1:jfirst]           .= f[:,n,1]
        f_new[:,i,jlast:end]          .= f[:,n,end]
    end
    f_new[:,ifirst:ilast,jfirst:jlast] = @view f[:,:,:]

    f_new[:,ilast:end,1:jfirst]       .= f[:,end,1]
    @fastmath @inbounds @threads for i in ilast:Nx_new
        f_new[:,i,jfirst:jlast]       .= @view f[:,end,:]
    end
    f_new[:,ilast:end,jlast:end]      .= f[:,end,end]

    nothing
end

function same_grid_spacing(grid::SpecCartGrid3D, boundary::Boundary , bulkevols::BulkPartition ,gauge::Gauge, chart2D::Chart)
    systems = SystemPartition(grid)
    Nsys    = length(systems)
    x_new   = systems[1].xcoord[:]
    y_new   = systems[1].ycoord[:]
    Nx_new  = length(x_new)
    Ny_new  = length(y_new)
    _, x, y = chart2D[:]

    if x_new[1] == x[1] && x_new[end] == x[end] && Nx_new == length(x)
        if y_new[1] == y[1] && y_new[end] == y[end] && Ny_new == length(y)
            return boundary, bulkevols, gauge
        end
    end
    boundary_new  = Boundary(grid)
    gauge_new     = Gauge(grid)
    bulkevols_new = BulkEvolvedPartition(grid)

    ifirst = findfirst(x_new .>= x[1])
    jfirst = findfirst(y_new .>= y[1])
    ilast  = findfirst(x_new .>= x[end])
    jlast  = findfirst(y_new .>= y[end])

    x_indices = (ifirst, ilast)
    y_indices = (jfirst, jlast)

    a4      = boundary.a4
    fx2     = boundary.fx2
    fy2     = boundary.fy2
    xi      = gauge.xi
    a4_new  = boundary_new.a4
    fx2_new = boundary_new.fx2
    fy2_new = boundary_new.fy2
    xi_new  = gauge_new.xi

    @time same_grid_spacing!(a4_new, a4, x_indices, y_indices)
    @time same_grid_spacing!(fx2_new, fx2, x_indices, y_indices)
    @time same_grid_spacing!(fy2_new, fy2, x_indices, y_indices)
    @time same_grid_spacing!(xi_new, xi, x_indices, y_indices)

    for i in 1:Nsys
        B1      = bulkevols[i].B1
        B2      = bulkevols[i].B2
        G       = bulkevols[i].G
        phi     = bulkevols[i].phi
        B1_new  = bulkevols_new[i].B1
        B2_new  = bulkevols_new[i].B2
        G_new   = bulkevols_new[i].G
        phi_new = bulkevols_new[i].phi

        @time same_grid_spacing!(B1_new, B1, x_indices, y_indices)
        @time same_grid_spacing!(B2_new, B2, x_indices, y_indices)
        @time same_grid_spacing!(G_new, G, x_indices, y_indices)
        @time same_grid_spacing!(phi_new, phi, x_indices, y_indices)

    end

    boundary_new, bulkevols_new, gauge_new

end

function different_grid_spacing(grid::SpecCartGrid3D, boundary::Boundary , bulkevols::BulkPartition ,gauge::Gauge,
                              chart2D::Chart)

    systems = SystemPartition(grid)
    Nsys    = length(systems)
    x_new   = systems[1].xcoord[:]
    y_new   = systems[1].ycoord[:]
    Nx_new  = length(x_new)
    Ny_new  = length(y_new)
    _, x, y = chart2D[:]

    if x_new[1] == x[1] && x_new[end] == x[end] && Nx_new == length(x)
        if y_new[1] == y[1] && y_new[end] == y[end] && Ny_new == length(y)
            return boundary, bulkevols, gauge
        end
    end

    interp = xy_interpolator(x,y)
    a4     = interp(boundary.a4)
    fx2    = interp(boundary.fx2)
    fy2    = interp(boundary.fy2)
    xi     = interp(gauge.xi)

    boundary_new  = Boundary(grid)
    gauge_new     = Gauge(grid)
    bulkevols_new = BulkEvolvedPartition(grid)

    a4_new  = boundary_new.a4
    fx2_new = boundary_new.fx2
    fy2_new = boundary_new.fy2
    xi_new  = gauge_new.xi

    ifirst = findfirst(x_new .>= x[1])
    jfirst = findfirst(y_new .>= y[1])
    ilast  = findfirst(x_new .>= x[end])
    jlast  = findfirst(y_new .>= y[end])

    for i in 1:Nsys
        B1  = interp(bulkevols[i].B1)
        B2  = interp(bulkevols[i].B2)
        G   = interp(bulkevols[i].G)
        phi = interp(bulkevols[i].phi)

        B1_new  = bulkevols_new[i].B1
        B2_new  = bulkevols_new[i].B2
        G_new   = bulkevols_new[i].G
        phi_new = bulkevols_new[i].phi

        @time B1_new[:,1:ifirst-1,1:jfirst-1]  .= B1(x[1], y[1])
        @time B2_new[:,1:ifirst-1,1:jfirst-1]  .= B2(x[1], y[1])
        @time G_new[:,1:ifirst-1,1:jfirst-1]   .= G(x[1], y[1])
        @time phi_new[:,1:ifirst-1,1:jfirst-1] .= phi(x[1], y[1])

        @time B1_new[:,1:ifirst-1,jlast+1:end]  .= B1(x[1], y[1])
        @time B2_new[:,1:ifirst-1,jlast+1:end]  .= B2(x[1], y[1])
        @time G_new[:,1:ifirst-1,jlast+1:end]   .= G(x[1], y[1])
        @time phi_new[:,1:ifirst-1,jlast+1:end] .= phi(x[1], y[1])

        @time B1_new[:,ilast+1:end,1:jfirst-1]  .= B1(x[1], y[1])
        @time B2_new[:,ilast+1:end,1:jfirst-1]  .= B2(x[1], y[1])
        @time G_new[:,ilast+1:end,1:jfirst-1]   .= G(x[1], y[1])
        @time phi_new[:,ilast+1:end,1:jfirst-1] .= phi(x[1], y[1])

        @time B1_new[:,ilast+1:end,jlast+1:end]  .= B1(x[1], y[1])
        @time B2_new[:,ilast+1:end,jlast+1:end]  .= B2(x[1], y[1])
        @time G_new[:,ilast+1:end,jlast+1:end]   .= G(x[1], y[1])
        @time phi_new[:,ilast+1:end,jlast+1:end] .= phi(x[1], y[1])

        if i == 1
            @time a4_new[:,1:ifirst-1,1:jfirst-1]  .= a4(x[1], y[1])
            @time fx2_new[:,1:ifirst-1,1:jfirst-1] .= fx2(x[1], y[1])
            @time fy2_new[:,1:ifirst-1,1:jfirst-1] .= fy2(x[1], y[1])
            @time xi_new[:,1:ifirst-1,1:jfirst-1]  .= xi(x[1], y[1])

            @time a4_new[:,1:ifirst-1,jlast+1:end]  .= a4(x[1], y[1])
            @time fx2_new[:,1:ifirst-1,jlast+1:end] .= fx2(x[1], y[1])
            @time fy2_new[:,1:ifirst-1,jlast+1:end] .= fy2(x[1], y[1])
            @time xi_new[:,1:ifirst-1,jlast+1:end]  .= xi(x[1], y[1])

            @time a4_new[:,ilast+1:end,1:jfirst-1]  .= a4(x[1], y[1])
            @time fx2_new[:,ilast+1:end,1:jfirst-1] .= fx2(x[1], y[1])
            @time fy2_new[:,ilast+1:end,1:jfirst-1] .= fy2(x[1], y[1])
            @time xi_new[:,ilast+1:end,1:jfirst-1]  .= xi(x[1], y[1])

            @time a4_new[:,ilast+1:end,jlast+1:end]  .= a4(x[1], y[1])
            @time fx2_new[:,ilast+1:end,jlast+1:end] .= fx2(x[1], y[1])
            @time fy2_new[:,ilast+1:end,jlast+1:end] .= fy2(x[1], y[1])
            @time xi_new[:,ilast+1:end,jlast+1:end]  .= xi(x[1], y[1])
        end

        @time @fastmath @inbounds @threads for j in ifirst:ilast
            B1_new[:,j,1:jfirst-1]  .= B1(x_new[j], y[1])
            B2_new[:,j,1:jfirst-1]  .= B2(x_new[j], y[1])
            G_new[:,j,1:jfirst-1]   .= G(x_new[j], y[1])
            phi_new[:,j,1:jfirst-1] .= phi(x_new[j], y[1])

            B1_new[:,j,jlast+1:end]  .= B1(x_new[j], y[1])
            B2_new[:,j,jlast+1:end]  .= B2(x_new[j], y[1])
            G_new[:,j,jlast+1:end]   .= G(x_new[j], y[1])
            phi_new[:,j,jlast+1:end] .= phi(x_new[j], y[1])
            if i == 1
                a4_new[:,j,1:jfirst-1]  .= a4(x_new[j], y[1])
                fx2_new[:,j,1:jfirst-1] .= fx2(x_new[j], y[1])
                fy2_new[:,j,1:jfirst-1] .= fy2(x_new[j], y[1])
                xi_new[:,j,1:jfirst-1]  .= xi(x_new[j], y[1])

                a4_new[:,j,jlast+1:end]  .= a4(x_new[j], y[1])
                fx2_new[:,j,jlast+1:end] .= fx2(x_new[j], y[1])
                fy2_new[:,j,jlast+1:end] .= fy2(x_new[j], y[1])
                xi_new[:,j,jlast+1:end]  .= xi(x_new[j], y[1])
            end
        end

        @time @fastmath @inbounds @threads for k in jfirst:jlast
            B1_new[:,1:ifirst-1,k]  .= B1(x[1], y_new[k])
            B2_new[:,1:ifirst-1,k]  .= B2(x[1], y_new[k])
            G_new[:,1:ifirst-1,k]   .= G(x[1], y_new[k])
            phi_new[:,1:ifirst-1,k] .= phi(x[1], y_new[k])

            B1_new[:,ilast+1:end,k]  .= B1(x[1], y_new[k])
            B2_new[:,ilast+1:end,k]  .= B2(x[1], y_new[k])
            G_new[:,ilast+1:end,k]   .= G(x[1], y_new[k])
            phi_new[:,ilast+1:end,k] .= phi(x[1], y_new[k])
            if i == 1
                a4_new[:,1:ifirst-1,k]  .= a4(x[1], y_new[k])
                fx2_new[:,1:ifirst-1,k] .= fx2(x[1], y_new[k])
                fy2_new[:,1:ifirst-1,k] .= fy2(x[1], y_new[k])
                xi_new[:,1:ifirst-1,k]  .= xi(x[1], y_new[k])

                a4_new[:,ilast+1:end,k]  .= a4(x[1], y_new[k])
                fx2_new[:,ilast+1:end,k] .= fx2(x[1], y_new[k])
                fy2_new[:,ilast+1:end,k] .= fy2(x[1], y_new[k])
                xi_new[:,ilast+1:end,k]  .= xi(x[1], y_new[k])
            end
        end

        @time @fastmath @inbounds
         @threads for k in jfirst:jlast
            for j in ifirst:ilast
                B1_new[:,j,k]  = B1(x_new[j], y_new[k])
                B2_new[:,j,k]  = B2(x_new[j], y_new[k])
                G_new[:,j,k]   = G(x_new[j], y_new[k])
                phi_new[:,j,k] = phi(x_new[j], y_new[k])
                if i == 1
                    a4_new[:,j,k]  = a4(x_new[j], y_new[k])
                    fx2_new[:,j,k] = fx2(x_new[j], y_new[k])
                    fy2_new[:,j,k] = fy2(x_new[j], y_new[k])
                    xi_new[:,j,k]  = xi(x_new[j], y_new[k])
                end
            end
        end
     end
    return boundary_new, bulkevols_new, gauge_new
end

function new_box(grid::SpecCartGrid3D, io::InOut, potential::Potential;
                            same_spacing::Symbol = :no)
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

    it       = max(it_bulk, it_gauge, it_bulk)
    t        = max(t_bulk, t_gauge, t_bulk)
    tinfo    = Jecco.TimeInfo(it, t, 0.0, 0.0)

    if same_spacing == :yes
        boundary, bulkevols, gauge = same_grid_spacing(grid, boundary, bulkevols, gauge, chart2D)
    else
        boundary, bulkevols, gauge = different_grid_spacing(grid, boundary, bulkevols, gauge, chart2D)
    end

    phi0 = try
        ts.params["phi0"]
    catch e
        if isa(e, KeyError)
            0.0   # if "phi0" is not found in the params Dict, set phi0 = 0
        else
            throw(e)
        end
    end
    oophiM2 = try
        ts.params["oophiM2"]
    catch e
        if isa(e, KeyError)
            0.0   # if "oophiM2" is not found in the params Dict, set oophiM2 = 0
        else
            throw(e)
        end
    end


    atlas   = Atlas(grid)
    systems = SystemPartition(grid)
    empty   = Cartesian{1}("u", 0.0, 0.0, 1)
    chart2D = Chart(empty, systems[1].xcoord, systems[1].ycoord)

    evolvars = AdS5_3_1.EvolVars(boundary, gauge, bulkevols)
    create_outputs(io.out_dir, evolvars, chart2D, Tuple(charts), io, potential, phi0
                                ,t=t, it=it)

end

function join_boxes(boundary_1::Boundary, boundary_2::Boundary)
    _, Nx_1, Ny = size(boundary_1.a4)
    _, Nx_2, _  = size(boundary_2.a4)
    Nx       = Nx_1 + Nx_2

    a4  = zeros(1, Nx, Ny)
    fy2 = zeros(1, Nx, Ny)
    fx2 = zeros(1, Nx, Ny)

    a4[1,1:Nx_1,:]      = boundary_1.a4
    a4[1,Nx_1+1:end,:]  = boundary_2.a4
    fx2[1,1:Nx_1,:]     = boundary_1.fx2
    fx2[1,Nx_1+1:end,:] = boundary_2.fx2
    fy2[1,1:Nx_1,:]     = boundary_1.fy2
    fy2[1,Nx_1+1:end,:] = boundary_2.fy2

    Boundary(a4, fx2, fy2)
end

function join_boxes(gauge_1::Gauge, gauge_2::Gauge)
    _, Nx_1, Ny = size(gauge_1.xi)
    _, Nx_2, _  = size(gauge_2.xi)
    Nx       = Nx_1 + Nx_2

    xi = zeros(1, Nx, Ny)

    xi[1,1:Nx_1,:]     = gauge_1.xi
    xi[1,Nx_1+1:end,:] = gauge_2.xi

    Gauge(xi)
end

function join_boxes(bulkevols_1::BulkPartition, bulkevols_2::BulkPartition)

    if length(bulkevols_1) != length(bulkevols_2)
        @warn "Both bulk grids must be identical"
        return
    end

    Nsys          = length(bulkevols_1)
    _, Nx_1, Ny   = size(bulkevols_1[1].B1)
    _, Nx_2, _    = size(bulkevols_2[1].B1)
    Nx            = Nx_1 + Nx_2
    bulk          = Array{BulkEvolved,1}(undef, Nsys)

    for n in 1:Nsys
        Nu, _, _  = size(bulkevols_1[n].B1)
        B1        = zeros(Nu, Nx, Ny)
        B2        = zeros(Nu, Nx, Ny)
        G         = zeros(Nu, Nx, Ny)
        phi       = zeros(Nu, Nx, Ny)

        B1[:,1:Nx_1,:]      = bulkevols_1[n].B1
        B1[:,Nx_1+1:end,:]  = bulkevols_2[n].B1
        B2[:,1:Nx_1,:]      = bulkevols_1[n].B2
        B2[:,Nx_1+1:end,:]  = bulkevols_2[n].B2
        G[:,1:Nx_1,:]       = bulkevols_1[n].G
        G[:,Nx_1+1:end,:]   = bulkevols_2[n].G
        phi[:,1:Nx_1,:]     = bulkevols_1[n].phi
        phi[:,Nx_1+1:end,:] = bulkevols_2[n].phi

        bulk[n] = BulkEvolved(B1, B2, G, phi)
    end

    AdS5_3_1.BulkPartition((bulk...))
end

function join_boxes(chart2D_1::Chart, chart2D_2::Chart)
    xmin_1   = chart2D_1.coords[2].min
    xmax_1   = chart2D_1.coords[2].max
    dx_1     = Jecco.delta(chart2D_1.coords[2])
    xmin_2   = chart2D_2.coords[2].min
    xmax_2   = chart2D_2.coords[2].max
    Nx_1     = chart2D_1.coords[2].nodes
    Nx_2     = chart2D_2.coords[2].nodes
    xmax     = xmax_2-xmin_2
    xmin     = xmin_1-(xmax_1+dx_1)
    Nx       = Nx_1+Nx_2
    xcoord   = Cartesian{2}("x", xmin, xmax, Nx)
    ycoord   = chart2D_1.coords[3]
    empty    = Cartesian{1}("u", 0.0, 0.0, 1)

    Chart(empty, xcoord, ycoord)
end

function join_boxes(charts_1::Array{Chart,1}, chart2D::Chart)
    Nsys   = length(charts_1)
    charts = Array{Chart,1}(undef, Nsys)
    for n in 1:Nsys
        ucoord   = charts_1[n].coords[1]
        xcoord   = chart2D.coords[2]
        ycoord   = chart2D.coords[3]
        charts[n] = Chart(ucoord, xcoord, ycoord)
    end

    Tuple(charts)
end

function join_boxes(io::InOut, potential::Potential, dir1::String, dir2::String)
    ts                    = OpenPMDTimeSeries(dir1, prefix="boundary_")
    boundary_1, chart2D_1 = construct_boundary(ts, ts.iterations[end])
    ts                    = OpenPMDTimeSeries(dir2, prefix="boundary_")
    boundary_2, chart2D_2 = construct_boundary(ts, ts.iterations[end])

    ts         = OpenPMDTimeSeries(dir1, prefix="gauge_")
    gauge_1, _ = construct_gauge(ts, ts.iterations[end])
    ts         = OpenPMDTimeSeries(dir2, prefix="gauge_")
    gauge_2, _ = construct_gauge(ts, ts.iterations[end])

    ts                    = OpenPMDTimeSeries(dir1, prefix="bulk_")
    bulkevols_1, charts_1 = construct_bulkevols(ts, ts.iterations[end])
    ts                    = OpenPMDTimeSeries(dir2, prefix="bulk_")
    bulkevols_2, charts_2 = construct_bulkevols(ts, ts.iterations[end])
#=
    if charts_1[1].coords[1] != charts_2[1].coords[1] || charts_1[2].coords[1] != chart_2[2].coords[1]
        @warn "Bulk grids must coincide for both data sets"
        return
    end
=#
    phi0 = try
        ts.params["phi0"]
    catch e
        if isa(e, KeyError)
            0.0   # if "phi0" is not found in the params Dict, set phi0 = 0
        else
            throw(e)
        end
    end
    chart2D   = join_boxes(chart2D_1, chart2D_2)
    charts    = join_boxes(charts_1, chart2D)
    boundary  = join_boxes(boundary_1, boundary_2)
    gauge     = join_boxes(gauge_1, gauge_2)
    bulkevols = join_boxes(bulkevols_1, bulkevols_2)
    evolvars  = AdS5_3_1.EvolVars(boundary, gauge, bulkevols)

    create_outputs(io.out_dir, evolvars, chart2D, charts, io, potential, phi0)
end

function shift!(boundary::Boundary, bulkevols::BulkPartition, gauge::Gauge, chart2D::Chart, new_center::Tuple{T,T}) where {T<:Real}

    _, x, y    = chart2D[:]
    Nx      = length(x)
    Ny      = length(y)
    Nsys    = length(bulkevols)
    x0      = (x[end]+x[2])/Nx
    y0      = (y[end]+y[2])/Ny
    Lx      = x[end]-2*x[1]+x[2]
    Ly      = y[end]-2*y[1]+y[2]
    x0p, y0p= new_center
    dx      = x0p-x0
    dy      = y0p-y0
    xmin    = 2*x0-x0p-Lx/2
    xmax    = 2*x0-x0p+Lx/2
    ymin    = 2*y0-y0p-Ly/2
    ymax    = 2*y0-y0p+Ly/2
    interp  = xy_interpolator(x,y)
    a4      = boundary.a4
    fx2     = boundary.fx2
    fy2     = boundary.fy2
    xi      = gauge.xi

    a4_inter  = interp(a4)
    fx2_inter = interp(fx2)
    fy2_inter = interp(fy2)
    xi_inter  = interp(xi)

    qx(x) = -sign(x+x0p-2*x0)*Int(floor(2*abs(x+x0p-2*x0)/Lx))*Lx
    qy(y) = -sign(y+y0p-2*y0)*Int(floor(2*abs(y+y0p-2*y0)/Ly))*Ly

    for n in 1:Nsys
        B1  = bulkevols[n].B1
        B2  = bulkevols[n].B2
        G   = bulkevols[n].G
        phi = bulkevols[n].phi

        B1_inter  = interp(B1)
        B2_inter  = interp(B2)
        G_inter   = interp(G)
        phi_inter = interp(phi)
        Threads.@threads for j in 1:Ny
            for i in 1:Nx
                B1[:,i,j]  = B1_inter(x[i]+dx, y[j]+dy)
                B2[:,i,j]  = B2_inter(x[i]+dx, y[j]+dy)
                G[:,i,j]   = G_inter(x[i]+dx, y[j]+dy)
                phi[:,i,j] = phi_inter(x[i]+dx, y[j]+dy)
                if n==1
                    a4[:,i,j]  = a4_inter(x[i]+dx, y[j]+dy)
                    fx2[:,i,j] = fx2_inter(x[i]+dx, y[j]+dy)
                    fy2[:,i,j] = fy2_inter(x[i]+dx, y[j]+dy)
                    xi[:,i,j]  = xi_inter(x[i]+dx, y[j]+dy)
                end
            end
        end
    end
end

function shift(io::InOut, potential::Potential; new_center::Tuple{T,T}=(0,0)) where {T<:Real}
    read_dir     = io.recover_dir

    ts                = OpenPMDTimeSeries(read_dir, prefix="boundary_")
    boundary, chart2D = construct_boundary(ts, ts.iterations[end])

    ts       = OpenPMDTimeSeries(read_dir, prefix="gauge_")
    gauge, _ = construct_gauge(ts, ts.iterations[end])

    ts                = OpenPMDTimeSeries(read_dir, prefix="bulk_")
    bulkevols, charts = construct_bulkevols(ts, ts.iterations[end])

    if new_center == (0,0)
        new_center = (chart2D[:][2][1], chart2D[:][3][1])
    end

    shift!(boundary, bulkevols, gauge, chart2D, new_center)

    phi0 = try
        ts.params["phi0"]
    catch e
        if isa(e, KeyError)
            0.0   # if "phi0" is not found in the params Dict, set phi0 = 0
        else
            throw(e)
        end
    end

    oophiM2 = try
        ts.params["oophiM2"]
    catch e
        if isa(e, KeyError)
            0.0   # if "oophiM2" is not found in the params Dict, set oophiM2 = 0
        else
            throw(e)
        end
    end
    evolvars = AdS5_3_1.EvolVars(boundary, gauge, bulkevols)
    create_outputs(io.out_dir, evolvars, chart2D, Tuple(charts), io, potential, phi0)
end

function change_energy(io::InOut, e_new::T, potential::Potential) where {T<:Real}
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
        else
            throw(e)
        end
    end

    oophiM2 = try
        ts.params["oophiM2"]
    catch e
        if isa(e, KeyError)
            0.0   # if "oophiM2" is not found in the params Dict, set oophiM2 = 0
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
    #k      = (-4/3*(e_new-e0)*phi0^4-maximum(a4)+a40)/(a40-maximum(a4))
    #a4 = k*(a4.-maximum(a4)).+maximum(a4)
    a4  .= a4.+4/3*(e0-e_new)

    evolvars = AdS5_3_1.EvolVars(boundary, gauge, bulkevols)
    create_outputs(io.out_dir, evolvars, chart2D, Tuple(charts), io, potential, phi0)

end

#For the moment we use the same grid as in the phase separation file and we will use create_new_data to change anything at the end.
#We use the B's and G of the PS state, and we modify the scalar fields and boundary data. f2 to 0 for the moment. Center the low energy
#phase in the middle of the box, so that the metaestable state will lie outside
function bubble_expansion(grid::SpecCartGrid3D, io::InOut, potential::Potential, A_dir::String, B_dir::String, PS_dir::String)
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
    a4_B  = a4_PS[end,cold]
    xi_B  = xi_PS[end,cold]
    k_a4  = (a4_A-a4_B)/(a4_PS[end,hot]-a4_PS[end,cold])
    k_xi  = (xi_A-xi_B)/(xi_PS[end,hot]-xi_PS[end,cold])

    a4  = zeros(1,Nx,Ny)
    xi  = zeros(1,Nx,Ny)
    fx2 = zeros(1,Nx,Ny)
    fy2 = zeros(1,Nx,Ny)

    a4[1,:,:] = k_a4.*(a4_PS[end,:,:].-a4_PS[end,cold]).+a4_B
    xi[1,:,:] = k_xi.*(xi_PS[end,:,:].-xi_PS[end,cold]).+xi_B

    boundary  = Boundary{typeof(a4[1,1,1])}(a4,fx2,fy2)
    gauge     = Gauge{typeof(xi[1,1,1])}(xi)
    bulk      = Array{BulkEvolved,1}(undef, Nsys)

    for n in 1:Nsys
        phi_A  = BulkTimeSeries(A_dir, :phi, n)[end,:,1,1]
        phi_B  = BulkTimeSeries(B_dir, :phi, n)[end,:,1,1]
        phi_PS = BulkTimeSeries(PS_dir, :phi, n)
        phi_B .= phi_PS[end,:,cold]
        B1     = BulkTimeSeries(PS_dir, :B2, n)[end,:,:,:]
        B2     = BulkTimeSeries(PS_dir, :B2, n)[end,:,:,:]
        G      = BulkTimeSeries(PS_dir, :G, n)[end,:,:,:]
        phi    = similar(phi_PS[end,:,:,:])
        k      = (phi_A-phi_B)./(phi_PS[end,:,hot]-phi_PS[end,:,cold])
        Nu     = length(phi[:,1,1])
        for i in 1:Nu
            phi[i,:,:] = k[i].*(phi_PS[end,i,:,:].-phi_PS[end,i,cold]).+phi_B[i]
        end
        bulk[n] = BulkEvolved{typeof(B1[1,1,1])}(B1, B2, G, phi)
    end

    bulkevols = AdS5_3_1.BulkPartition((bulk...))
    boundary_new, bulkevols_new, gauge_new = new_box(grid, boundary, bulkevols, gauge, chart2D)

    phi0 = try
        a4_PS.ts.params["phi0"]
    catch e
        if isa(e, KeyError)
            0.0   # if "phi0" is not found in the params Dict, set phi0 = 0
        else
            throw(e)
        end
    end

    empty = Cartesian{1}("u", 0.0, 0.0, 1)
    chart2D = Chart(empty, systems[1].xcoord, systems[1].ycoord)
    evolvars = AdS5_3_1.EvolVars(boundary_new, gauge_new, bulkevols_new)
    create_outputs(io.out_dir, evolvars, chart2D, atlas.charts, io, potential, phi0)
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
    create_outputs(io.out_dir, evolvars, chart2D, atlas.charts, io, potential, phi0)

    close(fid)

end

#For the moment we change a4 and fx2, fy2 by multiplying the energy by some factor. We can change this so that we decide what
#profile is used for these functions. We mantain the cold phase to its original value. Enbed in a bigger box by extending the edge value
function create_new_data(grid::SpecCartGrid3D, io::InOut, new_parameters::NewParameters, potential::Potential)
    dirname   = io.recover_dir
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

    a4    = BoundaryTimeSeries(dirname,:a4)
    fx2   = BoundaryTimeSeries(dirname,:fx2)
    fy2   = BoundaryTimeSeries(dirname,:fy2)
    phi2  = BoundaryTimeSeries(dirname,:phi2)
    xi    = XiTimeSeries(dirname)
    _,x,y = get_coords(a4,1,:,:)
    xmin  = x[1]
    xmax  = x[end]
    ymin  = y[1]
    ymax  = y[end]
    Nx    = length(x)
    Ny    = length(y)
    #Old runs do not have ts.params field, I set it by hand to what we usually use.
    phi0 = try
        a4.ts.params["phi0"]
    catch e
        if isa(e, KeyError)
            0.0   # if "phi0" is not found in the params Dict, set phi0 = 0
        else
            throw(e)
        end
    end
    oophiM2 = try
        a4.ts.params["oophiM2"]
    catch e
        if isa(e, KeyError)
            0.0   # if "oophiM2" is not found in the params Dict, set oophiM2 = 0
        else
            throw(e)
        end
    end
    #potential = Phi8Potential(oophiM2=-1.0, oophiQ=0.1,)
    #try
    #    phi0  = a4.ts.params["phi0"]
    #    potential = Phi8Potential(oophiM2=a4.ts.params["oophiM2"], oophiQ=a4.ts.params["oophiQ"],)
    #catch
    #    @warn "No ts.params field, setting phi0=1.0, oophiM2=-1,0 and oophiQ=0.1"
    #end

    interp = xy_interpolator(x,y)

    a4_inter   = interp(a4[:,:,:])
    fx2_inter  = interp(fx2[:,:,:])
    fy2_inter  = interp(fy2[:,:,:])
    phi2_inter = interp(phi2[:,:,:])
    xi_inter   = interp(xi[:,:,:])

    a4_new   = zeros(Nx_new,Ny_new)
    fx2_new  = zeros(Nx_new,Ny_new)
    fy2_new  = zeros(Nx_new,Ny_new)
    phi2_new = zeros(Nx_new,Ny_new)

    u = Array{Array{Real,1},1}(undef,Nsys)
    println("We start the insertion in the new box")

    for i in 1:Nsys
        B1  = BulkTimeSeries(dirname, :B1, i)
        B2  = BulkTimeSeries(dirname, :B2, i)
        G   = BulkTimeSeries(dirname, :G, i)
        phi = BulkTimeSeries(dirname, :phi, i)
        u[i]   = systems[i].ucoord[:]

        B1_inter  = interp(B1[end,:,:,:])
        B2_inter  = interp(B2[end,:,:,:])
        G_inter   = interp(G[end,:,:,:])
        phi_inter = interp(phi[end,:,:,:])
        for l in 1:Ny_new
            for k in 1:Nx_new
                if xmin<=x_new[k]<=xmax && ymin<=y_new[l]<=ymax
                    bulkevols[i].B1[:,k,l]  = B1_inter(x_new[k],y_new[l])
                    bulkevols[i].B2[:,k,l]  = B2_inter(x_new[k],y_new[l])
                    bulkevols[i].G[:,k,l]   = G_inter(x_new[k],y_new[l])
                    bulkevols[i].phi[:,k,l] = phi_inter(x_new[k],y_new[l])
                    if i == 1
                        a4_new[k,l]     = a4_inter(x_new[k],y_new[l])[1]
                        fx2_new[k,l]    = fx2_inter(x_new[k],y_new[l])[1]
                        fy2_new[k,l]    = fy2_inter(x_new[k],y_new[l])[1]
                        phi2_new[k,l]   = phi2_inter(x_new[k],y_new[l])[1]
                        gauge.xi[1,k,l] = xi_inter(x_new[k],y_new[l])[1]
                    end
                elseif xmin<=x_new[k]<=xmax
                    bulkevols[i].B1[:,k,l]  = B1_inter(x_new[k],ymin)
                    bulkevols[i].B2[:,k,l]  = B2_inter(x_new[k],ymin)
                    bulkevols[i].G[:,k,l]   = G_inter(x_new[k],ymin)
                    bulkevols[i].phi[:,k,l] = phi_inter(x_new[k],ymin)
                    if i == 1
                        a4_new[k,l]     = a4_inter(x_new[k],ymin)[1]
                        fx2_new[k,l]    = fx2_inter(x_new[k],ymin)[1]
                        fy2_new[k,l]    = fy2_inter(x_new[k],ymin)[1]
                        phi2_new[k,l]   = phi2_inter(x_new[k],ymin)[1]
                        gauge.xi[1,k,l] = xi_inter(x_new[k],ymin)[1]
                    end
                elseif ymin<=y_new[l]<=ymax
                    bulkevols[i].B1[:,k,l]  = B1_inter(xmin,y_new[l])
                    bulkevols[i].B2[:,k,l]  = B2_inter(xmin,y_new[l])
                    bulkevols[i].G[:,k,l]   = G_inter(xmin,y_new[l])
                    bulkevols[i].phi[:,k,l] = phi_inter(xmin,y_new[l])
                    if i == 1
                        a4_new[k,l]     = a4_inter(xmin,y_new[l])[1]
                        fx2_new[k,l]    = fx2_inter(xmin,y_new[l])[1]
                        fy2_new[k,l]    = fy2_inter(xmin,y_new[l])[1]
                        phi2_new[k,l]   = phi2_inter(xmin,y_new[l])[1]
                        gauge.xi[1,k,l] = xi_inter(xmin,y_new[l])[1]
                    end
                else
                    bulkevols[i].B1[:,k,l]  = B1[end,:,1,1]
                    bulkevols[i].B2[:,k,l]  = B2[end,:,1,1]
                    bulkevols[i].G[:,k,l]   = G[end,:,1,1]
                    bulkevols[i].phi[:,k,l] = phi[end,:,1,1]
                    if i == 1
                        a4_new[k,l]   = a4[end,1,1]
                        fx2_new[k,l]  = fx2[end,1,1]
                        fy2_new[k,l]  = fy2[end,1,1]
                        phi2_new[k,l] = phi2[end,1,1]
                        gauge.xi[1,k,l] = xi[end,1,1]
                    end
                end
            end
        end
    end
    println("Already set in the new box")
    e       = AdS5_3_1.compute_energy.(a4_new, phi2_new, phi0, oophiM2)
    e_new   = new_parameters.e_new
    fx20    = new_parameters.fx20
    fy20    = new_parameters.fy20
#change of energy
    if e_new != -1.0
        plan   = plan_rfft(e)
        e0     = real(1/(Nx_new*Ny_new)*(plan*e)[1])
        a40    = real(1/(Nx_new*Ny_new)*(plan*a4_new)[1])
        #k      = (-4/3*(e_new-e0)*phi0^4-maximum(a4_new)+a40)/(a40-maximum(a4_new))
        #a4_new = k*(a4_new.-maximum(a4_new)).+maximum(a4_new)
        a4_new  = a4_new.+4/3*(e0-e_new)
    end
#boosting, with static cold phase (edge).
    if new_parameters.boostx
        boundary.fx2[1,:,:] .= fx2_new
        fx2 = fx20*(e)#.-minimum(e))
        kappa = 0.1/fx20
    else
        boundary.fx2[1,:,:] .= fx2_new
    end
    if new_parameters.boosty
        boundary.fy2[1,:,:] .= fy20.*(e.-minimum(e))
    else
        boundary.fy2[1,:,:] .= fy2_new
    end
#We add the specified a4 cos perturbation
    ampx    = new_parameters.a4_ampx
    ampy    = new_parameters.a4_ampy
    kx      = new_parameters.a4_kx
    ky      = new_parameters.a4_ky
    xmin    = grid.x_min
    xmax    = grid.x_max
    ymin    = grid.y_min
    ymax    = grid.y_max
    xmid    = (xmax+xmin)/2
    ymid    = (ymax+ymin)/2

    plan   = plan_rfft(a4_new)
    a40    = real(1/(Nx_new*Ny_new)*(plan*a4_new)[1])
    for j in 1:Ny_new
        for i in 1:Nx_new
            x = x_new[i]
            y = y_new[j]

            boundary.a4[1,i,j] = a4_new[i,j]-a40*(ampx*cos(2*π*kx*(x-xmid)/(xmax-xmin))
                                                  +ampy*cos(2*π*ky*(y-ymid)/(ymax-ymin)))
        end
    end

#=
    #TODO: Change gauge
    u_AH = new_parameters.u_AH
    evoleq = AffineNull(phi0=phi0, potential=potential, gaugecondition = ConstantAH(u_AH = 1.0),)

    sigma  = 0.5*ones(1,Nx_new,Ny_new)
    change_gauge!(sigma, grid, boundary, gauge, bulkevols, evoleq, systems, u_AH)
    epsilon = 0
    println("epsilon = $epsilon")
    evoleq  = AffineNull(phi0=phi0, potential=potential, gaugecondition = ConstantAH(u_AH = u_AH),)

    while epsilon < 1.0
        n     = 0
        #kappa = kappa/0.8
        sigma = 1/(2.0*u[end][end])*ones(1,Nx_new,Ny_new)
        while maximum(1 ./sigma) >= u[end][end]+0.005
            if n != 0
                kappa = kappa*0.8
            end
            n += 1
            aux                 = epsilon+kappa
            boundary.fx2[1,:,:] = fx2*aux^0.4
#            boundary.fy2 .= fy2*aux^0.5
            try
                change_gauge!(sigma, grid, boundary, gauge, bulkevols, evoleq, systems, u_AH)
            catch

            end
        end
        epsilon += kappa
        println("epsilon = $epsilon")
        println("max fx2 = $(maximum(boundary.fx2))")
    end
=#
    evolvars = AdS5_3_1.EvolVars(boundary, gauge, bulkevols)
    create_outputs(io.out_dir, evolvars, chart2D, atlas.charts, io, potential, phi0)
    #return bulkevols, boundary, gauge
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
create_outputs(io.out_dir, evolvars, chart2D, atlas.charts, io, potential, phi0)
#return sigma
end

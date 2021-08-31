using FFTW
using HDF5

#It takes folders of continuation of the same runs in bigger and bigger boxes and piuts them inside the newest, and biggest box
#It puts the final result in the same folder, only boundary data
#TODO: Check possible discrepancy when enlarging the box. It might be that you are shifting
#TODO: the whole thing avery time you enlarge. It would not affect the physics though
#But check it and make sure it does not.
function same_box_all_runs(outdir::T, directories::Array{T,1}) where {T<:String, TP<:Real}
    try run(`mkdir $outdir`) catch end
    #run(`cp $(directories[end])/boundary_* $outdir`)
    tinfo           = Jecco.TimeInfo()
    ts              = OpenPMDTimeSeries(directories[end], prefix="boundary_")
    _, chart2D      = get_field(ts, it=ts.iterations[end], field="a4")
    _, x_new, y_new = chart2D[:]
    Nx_new, Ny_new  = (length(x_new), length(y_new))

    for (index, directory) in enumerate(directories[1:end-1])
        ts           = OpenPMDTimeSeries(directory, prefix="boundary_")
        it_max       = OpenPMDTimeSeries(directories[index+1], prefix="boundary_").iterations[1]
        iterations   = ts.iterations
        _, chart_old = get_field(ts, it=iterations[1], field="a4")
        xmin         = chart_old.coords[2].min
        xmax         = chart_old.coords[2].max
        ymin         = chart_old.coords[3].min
        ymax         = chart_old.coords[3].max
        ifirst       = findfirst(x_new .>= xmin)
        ilast        = findfirst(x_new .>= xmax)
        jfirst       = findfirst(y_new .>= ymin)
        jlast        = findfirst(y_new .>= ymax)
        x_indices    = (ifirst, ilast)
        y_indices    = (jfirst, jlast)

        for (idx, it) in enumerate(iterations)
            if it >= it_max break end
            a4       = get_field(ts, it=it, field="a4")[1]
            b14      = get_field(ts, it=it, field="b14")[1]
            b24      = get_field(ts, it=it, field="b24")[1]
            fx2      = get_field(ts, it=it, field="fx2")[1]
            fy2      = get_field(ts, it=it, field="fy2")[1]
            g4       = get_field(ts, it=it, field="g4")[1]
            phi2     = get_field(ts, it=it, field="phi2")[1]

            tinfo.it = it
            tinfo.t  = ts.current_t

            #=
            f_new              = zeros(1, Nx_new, Ny_new)
            boundary_fields    = Array{Jecco.Field, 1}(undef, 7)
            AdS5_3_1.same_grid_spacing!(f_new, a4, x_indices, y_indices)
            boundary_fields[1] = Jecco.Field("a4", f_new, chart2D)
            AdS5_3_1.same_grid_spacing!(f_new, b14, x_indices, y_indices)
            boundary_fields[2] = Jecco.Field("b14", f_new, chart2D)
            AdS5_3_1.same_grid_spacing!(f_new, b24, x_indices, y_indices)
            boundary_fields[3] = Jecco.Field("b24", f_new, chart2D)
            AdS5_3_1.same_grid_spacing!(f_new, fx2, x_indices, y_indices)
            boundary_fields[4] = Jecco.Field("fx2", f_new, chart2D)
            AdS5_3_1.same_grid_spacing!(f_new, fy2, x_indices, y_indices)
            boundary_fields[5] = Jecco.Field("fy2", f_new, chart2D)
            AdS5_3_1.same_grid_spacing!(f_new, g4, x_indices, y_indices)
            boundary_fields[6] = Jecco.Field("g4", f_new, chart2D)
            AdS5_3_1.same_grid_spacing!(f_new, phi2, x_indices, y_indices)
            boundary_fields[7] = Jecco.Field("phi2", f_new, chart2D)
            bdry_out           = Jecco.Output(outdir, "boundary_", tinfo)


            bdry_out((boundary_fields...,), params=ts.params)

            =#

            #TODO: You have to optimize this, a single f_new should be enough.
            #it doesn't work above...
            a4_new   = zeros(1, Nx_new, Ny_new)
            b14_new  = zeros(1, Nx_new, Ny_new)
            b24_new  = zeros(1, Nx_new, Ny_new)
            fx2_new  = zeros(1, Nx_new, Ny_new)
            fy2_new  = zeros(1, Nx_new, Ny_new)
            g4_new   = zeros(1, Nx_new, Ny_new)
            phi2_new = zeros(1, Nx_new, Ny_new)

            AdS5_3_1.same_grid_spacing!(a4_new, a4, x_indices, y_indices)
            AdS5_3_1.same_grid_spacing!(b14_new, b14, x_indices, y_indices)
            AdS5_3_1.same_grid_spacing!(b24_new, b24, x_indices, y_indices)
            AdS5_3_1.same_grid_spacing!(fx2_new, fx2, x_indices, y_indices)
            AdS5_3_1.same_grid_spacing!(fy2_new, fy2, x_indices, y_indices)
            AdS5_3_1.same_grid_spacing!(g4_new, g4, x_indices, y_indices)
            AdS5_3_1.same_grid_spacing!(phi2_new, phi2, x_indices, y_indices)

            bdry_out = Jecco.Output(outdir, "boundary_", tinfo)

            boundary_fields = (
                Jecco.Field("a4", a4_new, chart2D),
                Jecco.Field("b14", b14_new, chart2D),
                Jecco.Field("b24", b24_new, chart2D),
                Jecco.Field("fx2", fx2_new, chart2D),
                Jecco.Field("fy2", fy2_new, chart2D),
                Jecco.Field("g4", g4_new, chart2D),
                Jecco.Field("phi2", phi2_new, chart2D)
            )
            bdry_out(boundary_fields, params=ts.params)
        end
    end
end

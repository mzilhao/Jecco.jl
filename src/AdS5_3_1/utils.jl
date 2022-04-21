
"""
Need to have method `get_data`
"""
abstract type TimeSeries{N,T} end

struct BoundaryTimeSeries{T} <: TimeSeries{2,T}
    ts    :: T
    field :: Symbol

    function BoundaryTimeSeries(foldername::String, field::Symbol)
        ts = OpenPMDTimeSeries(foldername, "boundary_")
        new{typeof(ts)}(ts, field)
    end
end

struct XiTimeSeries{T} <: TimeSeries{2,T}
    ts :: T

    function XiTimeSeries(foldername::String)
        ts = OpenPMDTimeSeries(foldername, "gauge_")
        new{typeof(ts)}(ts)
    end
end

struct BulkTimeSeries{T} <: TimeSeries{3,T}
    ts        :: T
    field     :: Symbol
    component :: Int

    function BulkTimeSeries(foldername::String, field::Symbol, c::Int)
        ts = OpenPMDTimeSeries(foldername, "bulk_")
        new{typeof(ts)}(ts, field, c)
    end
end

struct ConstrainedTimeSeries{T} <: TimeSeries{3,T}
    ts        :: T
    field     :: Symbol
    component :: Int

    function ConstrainedTimeSeries(foldername::String, field::Symbol, c::Int)
        ts = OpenPMDTimeSeries(foldername, "constrained_")
        new{typeof(ts)}(ts, field, c)
    end
end

struct VEVTimeSeries{T} <: TimeSeries{2,T}
    ts  :: T
    vev :: Symbol

    function VEVTimeSeries(foldername::String, vev::Symbol)
        ts = OpenPMDTimeSeries(foldername, "boundary_")
        new{typeof(ts)}(ts, vev)
    end
end

struct AHTimeSeries{T} <: TimeSeries{2,T}
    ts            :: T
    ts_diag       :: T
    property      :: Symbol

    function AHTimeSeries(foldername::String, property::Symbol)
        ts       = OpenPMDTimeSeries(foldername, "constrained_")
        ts_diag  = OpenPMDTimeSeries(foldername, "diagnostics_")
        new{OpenPMDTimeSeries}(ts, ts_diag, property)
    end
end

struct LocalVEVsTimeSeries{T} <: TimeSeries{2,T}
    ts    :: T
    field :: Symbol
    function LocalVEVsTimeSeries(foldername::String, field::Symbol)
        ts = OpenPMDTimeSeries(foldername, "boundary_")
        new{typeof(ts)}(ts, field)
    end
end

struct IdealHydroTimeSeries{T} <: TimeSeries{2,T}
    ts  :: T
    vev :: Symbol
    eos :: String

    function IdealHydroTimeSeries(foldername::String, path_to_eos::String, vev::Symbol)
        ts = OpenPMDTimeSeries(foldername, "boundary_")
        new{typeof(ts)}(ts, vev, path_to_eos)
    end
end


function get_data(ff::BoundaryTimeSeries, it::Int)
    f, chart = get_field(ff.ts, it=it, field=String(ff.field))
    _, x, y  = chart[:]
    f[1,:,:], [x, y]
end

function get_data(xi::XiTimeSeries, it::Int)
    f, chart = get_field(xi.ts, it=it, field="xi")
    _, x, y  = chart[:]
    f[1,:,:], [x, y]
end

function get_data(ff::BulkTimeSeries, it::Int)
    field = "$(ff.field) c=$(ff.component)"
    f, chart = get_field(ff.ts, it=it, field=field)
    u, x, y = chart[:]
    f, [u, x, y]
end

function get_data(ff::ConstrainedTimeSeries, it::Int)
    field = "$(ff.field) c=$(ff.component)"
    f, chart = get_field(ff.ts, it=it, field=field)
    u, x, y = chart[:]
    f, [u, x, y]
end

function get_data(ff::VEVTimeSeries, it::Int)
    if ff.vev == :energy
        f, chart = get_energy(ff.ts, it=it)
    elseif ff.vev == :Jx
        f, chart = get_Jx(ff.ts, it=it)
    elseif ff.vev == :Jy
        f, chart = get_Jy(ff.ts, it=it)
    elseif ff.vev == :px
        f, chart = get_px(ff.ts, it=it)
    elseif ff.vev == :py
        f, chart = get_py(ff.ts, it=it)
    elseif ff.vev == :pz
        f, chart = get_pz(ff.ts, it=it)
    elseif ff.vev == :Jx
        f, chart = get_Jx(ff.ts, it=it)
    elseif ff.vev == :pxy
        f, chart = get_pxy(ff.ts, it=it)
    elseif ff.vev == :Ophi
        f, chart = get_Ophi(ff.ts, it=it)
    else
        error("Unknown VEV")
    end
    _, x, y  = chart[:]
    f[1,:,:], [x, y]
end

function get_data(ff::AHTimeSeries, it::Int)
    if ff.property == :temperature
        f, chart = get_temperature(ff.ts, ff.ts_diag, it=it)
    elseif ff.property == :entropy
        f, chart = get_entropy(ff.ts, ff.ts_diag, it=it)
    else
        error("Unknown property")
    end
    _, x, y = chart[:]
    f[1,:,:], [x,y]
end

function Base.size(ff::TimeSeries)
    Nt  = length(ff.ts.iterations)
    it0 = ff.ts.iterations[1]
    f0, = get_data(ff, it0)
    size_ = size(f0)
    Nt, size_...
end

function get_data(ff::LocalVEVsTimeSeries, it::Int)
    e, chart  = get_energy(ff.ts, it=it)
    Jx, _     = get_Jx(ff.ts, it=it)
    Jy, _     = get_Jy(ff.ts, it=it)
    px, _     = get_px(ff.ts, it=it)
    pxy, _    = get_pxy(ff.ts, it=it)
    py, _     = get_py(ff.ts, it=it)

    ut, ux, uy, el, p1, p2 = compute_local_VEVs(e, Jx, Jy, px, pxy, py)

    if ff.field == :ut
        f = ut
    elseif ff.field == :ux
        f = ux
    elseif ff.field == :uy
        f = uy
    elseif ff.field == :vx
        f = ux./ut
    elseif ff.field == :vy
        f = uy./ut
    elseif ff.field == :energy
        f = el
    elseif ff.field == :p1
        f = p1
    elseif ff.field == :p2
        f = p2
    else
        error("Unknown local field")
    end

    x, y = chart[:]

    f[:,:], [x, y]
end

function get_data(ff::IdealHydroTimeSeries, it::Int)
    ee        = h5read(ff.eos, "Energy")
    pp        = h5read(ff.eos, "Pressure")
    e, chart  = get_energy(ff.ts, it=it)
    Jx, _     = get_Jx(ff.ts, it=it)
    Jy, _     = get_Jy(ff.ts, it=it)
    px, _     = get_px(ff.ts, it=it)
    pxy, _    = get_pxy(ff.ts, it=it)
    py, _     = get_py(ff.ts, it=it  )
    _, Nx, Ny = size(e)
    p         = zeros(Nx, Ny)

    ut, ux, uy, el, _, _ = compute_local_VEVs(e, Jx, Jy, px, pxy, py)

    @inbounds @threads for j in 1:Ny
        for i in 1:Nx
            e1 = findlast(ee .<= el[i,j])
            e2 = findfirst(ee .>= el[i,j])
            if e1 != e2
                p[i,j] = (pp[e2]-pp[e1])/(ee[e2]-ee[e1])*(el[i,j]-ee[e1])+pp[e1]
            else
                p[i,j] = pp[e1]
            end
        end
    end

    if ff.vev == :energy
        f = e_Ideal.(ut, el, p)
    elseif ff.vev == :Jx
        f = Jx_Ideal.(ut, ux, el, p)
    elseif ff.vev == :Jy
        f = Jy_Ideal.(ut, uy, el, p)
    elseif ff.vev == :px
        f = px_Ideal.(ux, el, p)
    elseif ff.vev == :pxy
        f = pxy_Ideal.(ux, uy, el, p)
    elseif ff.vev == :py
        f = py_Ideal.(uy, el, p)
    elseif ff.vev == :pz
        f = p
    else
        error("Unknown VEV")
    end

    x, y = chart[:]

    f[:,:], [x, y]
end

@inline Base.firstindex(ff::TimeSeries, d::Int) = 1
@inline Base.lastindex(ff::TimeSeries, d::Int) = size(ff)[d]

function Base.getindex(ff::TimeSeries, a::Int, idx::Vararg)
    it = ff.ts.iterations[a]
    f, = get_data(ff, it)
    f[idx...]
end

function Base.getindex(ff::TimeSeries, aa, idx::Vararg)
    it  = ff.ts.iterations[aa[1]]
    f0, = get_data(ff, it)

    Na    = length(aa)
    size_ = size(f0[idx...])
    f     = zeros(Na, size_...)

    slicer = [Colon() for _ in 1:ndims(f0)]
    for (i,a) in enumerate(aa)
        it  = ff.ts.iterations[a]
        f0, = get_data(ff, it)
        f[i,slicer...] .= f0[idx...]
    end
    f
end

function Base.getindex(ff::TimeSeries, ::Colon, idx::Vararg)
    Na  = length(ff.ts.iterations)
    getindex(ff, 1:Na, idx...)
end


function Jecco.get_coords(ff::TimeSeries{N}, a::Int, idx::Vararg) where{N}
    it = ff.ts.iterations[a]
    f, coords = get_data(ff, it)
    t = ff.ts.current_t
    coord = [coords[a][idx[a]] for a in 1:N]
    t, coord...
end

function Jecco.get_coords(ff::TimeSeries{N}, aa, idx::Vararg) where{N}
    it = ff.ts.iterations[aa[1]]
    f0, coords = get_data(ff, it)
    t0 = ff.ts.current_t
    coord = [coords[a][idx[a]] for a in 1:N]

    Na = length(aa)
    t  = Vector{typeof(t0)}(undef, Na)
    for (i,a) in enumerate(aa)
        it   = ff.ts.iterations[a]
        get_data(ff, it)
        t0   = ff.ts.current_t
        t[i] = t0
    end
    t, coord...
end

function Jecco.get_coords(ff::TimeSeries{N}, ::Colon, idx::Vararg) where{N}
    Na = length(ff.ts.iterations)
    get_coords(ff, 1:Na, idx...)
end


function convert_to_mathematica(dirname::String; outfile::String="data_mathematica.h5")
    ts = OpenPMDTimeSeries(dirname, prefix="boundary_")

    iterations = ts.iterations

    Nt = length(iterations)
    t  = zeros(Nt)

    it = 0
    en0, chart = get_energy(ts, it=it)

    _, x, y = chart[:]
    Nx = length(x)
    Ny = length(y)

    en   = zeros(Nt,Nx,Ny)
    Jx   = zeros(Nt,Nx,Ny)
    Jy   = zeros(Nt,Nx,Ny)
    px   = zeros(Nt,Nx,Ny)
    py   = zeros(Nt,Nx,Ny)
    pz   = zeros(Nt,Nx,Ny)
    pxy  = zeros(Nt,Nx,Ny)
    Ophi = zeros(Nt,Nx,Ny)

    for (idx,it) in enumerate(iterations)
        en[idx,:,:]   .= get_energy(ts, it=it)[1][1,:,:]
        Jx[idx,:,:]   .= get_Jx(ts, it=it)[1][1,:,:]
        Jy[idx,:,:]   .= get_Jy(ts, it=it)[1][1,:,:]
        px[idx,:,:]   .= get_px(ts, it=it)[1][1,:,:]
        py[idx,:,:]   .= get_py(ts, it=it)[1][1,:,:]
        pz[idx,:,:]   .= get_pz(ts, it=it)[1][1,:,:]
        pxy[idx,:,:]  .= get_pxy(ts, it=it)[1][1,:,:]
        Ophi[idx,:,:] .= get_Ophi(ts, it=it)[1][1,:,:]

        t[idx] = ts.current_t
    end

    # store in an array suitable for Mathematica
    T_m = zeros(11, Nt*Nx*Ny)

    n = 1
    for i in 1:Nt
        for j in 1:Nx
            for k in 1:Ny
                T_m[:,n] = [t[i] x[j] y[k] en[i,j,k] Jx[i,j,k] Jy[i,j,k] px[i,j,k] pxy[i,j,k] py[i,j,k] pz[i,j,k] Ophi[i,j,k]]
                n += 1
            end
        end
    end

    output   = abspath(dirname, outfile)
    fid      = h5open(output, "w")
    group_st = g_create(fid, "data")
    Jecco.write_dataset(group_st, "VEVs", T_m)
    close(fid)
end

function convert_to_mathematica_y0(dirname::String; outfile::String="data_mathematica_y0.h5")
    ts = OpenPMDTimeSeries(dirname, prefix="boundary_")

    iterations = ts.iterations

    Nt = length(iterations)
    t  = zeros(Nt)

    it = 0
    en0, chart = get_energy(ts, it=it)

    _, x, y = chart[:]
    Nx      = length(x)
    Ny      = length(y)
    idy     = Int(floor(Ny/2))

    en   = zeros(Nt,Nx)
    Jx   = zeros(Nt,Nx)
    Jy   = zeros(Nt,Nx)
    px   = zeros(Nt,Nx)
    py   = zeros(Nt,Nx)
    pz   = zeros(Nt,Nx)
    pxy  = zeros(Nt,Nx)
    Ophi = zeros(Nt,Nx)

    for (idx,it) in enumerate(iterations)
        en[idx,:]   .= get_energy(ts, it=it)[1][1,:,idy]
        Jx[idx,:]   .= get_Jx(ts, it=it)[1][1,:,idy]
        Jy[idx,:]   .= get_Jy(ts, it=it)[1][1,:,idy]
        px[idx,:]   .= get_px(ts, it=it)[1][1,:,idy]
        py[idx,:]   .= get_py(ts, it=it)[1][1,:,idy]
        pz[idx,:]   .= get_pz(ts, it=it)[1][1,:,idy]
        pxy[idx,:]  .= get_pxy(ts, it=it)[1][1,:,idy]
        Ophi[idx,:] .= get_Ophi(ts, it=it)[1][1,:,idy]

        t[idx] = ts.current_t
    end

    # store in an array suitable for Mathematica
    T_m = zeros(10, Nt*Nx)

    n = 1
    for i in 1:Nt
        for j in 1:Nx
            T_m[:,n] = [t[i] x[j] en[i,j] Jx[i,j] Jy[i,j] px[i,j] pxy[i,j] py[i,j] pz[i,j] Ophi[i,j]]
            n += 1
        end
    end

    output   = abspath(dirname, outfile)
    fid      = h5open(output, "w")
    group_st = g_create(fid, "data")
    Jecco.write_dataset(group_st, "VEVs", T_m)
    close(fid)
end

function convert_to_mathematica_diagonal(dirname::String; outfile::String="data_mathematica_diagonal.h5")
    ts = OpenPMDTimeSeries(dirname, prefix="boundary_")

    iterations = ts.iterations

    Nt = length(iterations)
    t  = zeros(Nt)

    it = 0
    en0, chart = get_energy(ts, it=it)

    _, x, y = chart[:]
    Nx      = length(x)
    Ny      = length(y)

    # store in an array suitable for Mathematica
    T_m = zeros(10, Nt*Nx)
    n   = 1
    @inbounds for (idx,it) in enumerate(iterations)
        en    = interpolate((x,y,), get_energy(ts, it=it)[1][1,:,:],Gridded(Linear()))
        Jx    = interpolate((x,y,), get_Jx(ts, it=it)[1][1,:,:],Gridded(Linear()))
        Jy    = interpolate((x,y,), get_Jy(ts, it=it)[1][1,:,:],Gridded(Linear()))
        px    = interpolate((x,y,), get_px(ts, it=it)[1][1,:,:],Gridded(Linear()))
        py    = interpolate((x,y,), get_py(ts, it=it)[1][1,:,:],Gridded(Linear()))
        pz    = interpolate((x,y,), get_pz(ts, it=it)[1][1,:,:],Gridded(Linear()))
        pxy   = interpolate((x,y,), get_pxy(ts, it=it)[1][1,:,:],Gridded(Linear()))
        Ophi  = interpolate((x,y,), get_Ophi(ts, it=it)[1][1,:,:],Gridded(Linear()))
        t     = ts.current_t
        for j in 1:Nx
            r        = x[j]/sqrt(2)
            T_m[:,n] = [t x[j] en(r,r) Jx(r,r) Jy(r,r) px(r,r) pxy(r,r) py(r,r) pz(r,r) Ophi(r,r)]
            n       += 1
        end
    end

    output   = abspath(dirname, outfile)
    fid      = h5open(output, "w")
    group_st = g_create(fid, "data")
    Jecco.write_dataset(group_st, "VEVs", T_m)
    close(fid)
end

function convert_to_mathematica_local(dirname::String; outfile::String="local_data_mathematica.h5")#,
                                #phi0, oophiM2)

    ts = OpenPMDTimeSeries(dirname, prefix="boundary_")

    iterations = ts.iterations

    Nt = length(iterations)
    t  = zeros(Nt)

    it = 0
    en0, chart = get_energy(ts, it=it)

    _, x, y = chart[:]
    Nx = length(x)
    Ny = length(y)

    ut        = zeros(Nt,Nx,Ny)
    ux        = zeros(Nt,Nx,Ny)
    uy        = zeros(Nt,Nx,Ny)
    en_local  = zeros(Nt,Nx,Ny)
    p1_local  = zeros(Nt,Nx,Ny)
    p2_local  = zeros(Nt,Nx,Ny)

    for (idx,it) in enumerate(iterations)
        en   = get_energy(ts, it=it)[1][1,:,:]
        Jx   = get_Jx(ts, it=it)[1][1,:,:]
        Jy   = get_Jy(ts, it=it)[1][1,:,:]
        px   = get_px(ts, it=it)[1][1,:,:]
        py   = get_py(ts, it=it)[1][1,:,:]
        pxy  = get_pxy(ts, it=it)[1][1,:,:]

        t[idx] = ts.current_t

        for j in 1:Ny
            for i in 1:Nx
                ut[idx,i,j],ux[idx,i,j],uy[idx,i,j],en_local[idx,i,j],p1_local[idx,i,j],p2_local[idx,i,j] = AdS5_3_1.compute_local_VEVs(en[i,j],Jx[i,j],Jy[i,j],px[i,j],py[i,j],pxy[i,j])
            end
        end
    end

    # store in an array suitable for Mathematica
    T_m = zeros(9, Nt*Nx*Ny)

    n = 1
    for i in 1:Nt
        for j in 1:Nx
            for k in 1:Ny
                T_m[:,n] = [t[i] x[j] y[k] ut[i,j,k] ux[i,j,k] uy[i,j,k] en_local[i,j,k] p1_local[i,j,k] p2_local[i,j,k]]
                n += 1
            end
        end
    end

    output   = abspath(dirname, outfile)
    fid      = h5open(output, "w")
    group_st = g_create(fid, "data")
    Jecco.write_dataset(group_st, "Local VEVs", T_m)
    close(fid)


end

function convert_to_mathematica_local_y0(dirname::String; outfile::String="local_data_mathematica_y0.h5")
    ts         = OpenPMDTimeSeries(dirname, prefix="boundary_")
    iterations = ts.iterations
    Nt         = length(iterations)
    t          = zeros(Nt)
    it         = 0
    en0, chart = get_energy(ts, it=it)

    _, x, y = chart[:]
    Nx      = length(x)
    Ny      = length(y)
    idy     = Int(floor(Ny/2))

    ut       = zeros(Nt,Nx)
    ux       = zeros(Nt,Nx)
    uy       = zeros(Nt,Nx)
    en_local = zeros(Nt,Nx)
    p1_local = zeros(Nt,Nx)
    p2_local = zeros(Nt,Nx)

    @inbounds for (idx,it) in enumerate(iterations)
        en   = get_energy(ts, it=it)[1][1,:,idy]
        Jx   = get_Jx(ts, it=it)[1][1,:,idy]
        Jy   = get_Jy(ts, it=it)[1][1,:,idy]
        px   = get_px(ts, it=it)[1][1,:,idy]
        py   = get_py(ts, it=it)[1][1,:,idy]
        pxy  = get_pxy(ts, it=it)[1][1,:,idy]

        t[idx] = ts.current_t

        for i in 1:Nx
            ut[idx,i],ux[idx,i],uy[idx,i],en_local[idx,i],p1_local[idx,i],p2_local[idx,i] = AdS5_3_1.compute_local_VEVs(en[i],Jx[i],Jy[i],px[i],py[i],pxy[i])
        end
    end

    # store in an array suitable for Mathematica
    T_m = zeros(8, Nt*Nx*Ny)

    n = 1
    for i in 1:Nt
        for j in 1:Nx
            T_m[:,n] = [t[i] x[j] ut[i,j] ux[i,j] uy[i,j] en_local[i,j] p1_local[i,j] p2_local[i,j]]
            n += 1
        end
    end

    output   = abspath(dirname, outfile)
    fid      = h5open(output, "w")
    group_st = g_create(fid, "data")
    Jecco.write_dataset(group_st, "Local VEVs", T_m)
    close(fid)
end

function convert_to_mathematica_local_diagonal(dirname::String; outfile::String="local_data_mathematica_diagonal.h5")
    ts         = OpenPMDTimeSeries(dirname, prefix="boundary_")
    iterations = ts.iterations
    Nt         = length(iterations)
    t          = zeros(Nt)
    it         = 0
    en0, chart = get_energy(ts, it=it)
    _, x, y    = chart[:]
    Nx         = length(x)
    Ny         = length(y)
    ut         = zeros(Nt,Nx)
    ux         = zeros(Nt,Nx)
    uy         = zeros(Nt,Nx)
    en_local   = zeros(Nt,Nx)
    p1_local   = zeros(Nt,Nx)
    p2_local   = zeros(Nt,Nx)

    @inbounds for (idx,it) in enumerate(iterations)
        en     = interpolate((x,y,), get_energy(ts, it=it)[1][1,:,:],Gridded(Linear()))
        Jx     = interpolate((x,y,), get_Jx(ts, it=it)[1][1,:,:],Gridded(Linear()))
        Jy     = interpolate((x,y,), get_Jy(ts, it=it)[1][1,:,:],Gridded(Linear()))
        px     = interpolate((x,y,), get_px(ts, it=it)[1][1,:,:],Gridded(Linear()))
        py     = interpolate((x,y,), get_py(ts, it=it)[1][1,:,:],Gridded(Linear()))
        pz     = interpolate((x,y,), get_pz(ts, it=it)[1][1,:,:],Gridded(Linear()))
        pxy    = interpolate((x,y,), get_pxy(ts, it=it)[1][1,:,:],Gridded(Linear()))
        Ophi   = interpolate((x,y,), get_Ophi(ts, it=it)[1][1,:,:],Gridded(Linear()))
        t[idx] = ts.current_t

        @inbounds for i in 1:Nx
            ut[idx,i],ux[idx,i],uy[idx,i],en_local[idx,i],p1_local[idx,i],p2_local[idx,i] = AdS5_3_1.compute_local_VEVs(en(r,r),Jx(r,r),Jy(r,r),px(r,r),py(r,r),pxy(r,r))
        end
    end
    # store in an array suitable for Mathematica
    T_m = zeros(8, Nt*Nx*Ny)
    n   = 1
    for i in 1:Nt
        for j in 1:Nx
            r        = x[j]/sqrt(2)
            T_m[:,n] = [t[i] x[j] ut[i,j] ux[i,j] uy[i,j] en_local[i,j] p1_local[i,j] p2_local[i,j]]
            n       += 1
        end
    end

    output   = abspath(dirname, outfile)
    fid      = h5open(output, "w")
    group_st = g_create(fid, "data")
    Jecco.write_dataset(group_st, "Local VEVs", T_m)
    close(fid)
end

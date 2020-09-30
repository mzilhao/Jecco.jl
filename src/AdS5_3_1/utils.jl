using HDF5

"""
Need to have method `get_data`
"""
abstract type TimeSeries{N,T} end

struct BoundaryTimeSeries{N,T} <: TimeSeries{N,T}
    ts    :: T
    field :: Symbol

    function BoundaryTimeSeries(foldername::String, field::Symbol)
        ts = OpenPMDTimeSeries(foldername, "boundary_")
        new{3,typeof(ts)}(ts, field)
    end
end

struct XiTimeSeries{N,T} <: TimeSeries{N,T}
    ts :: T

    function XiTimeSeries(foldername::String)
        ts = OpenPMDTimeSeries(foldername, "gauge_")
        new{3,typeof(ts)}(ts)
    end
end

struct BulkTimeSeries{N,T} <: TimeSeries{N,T}
    ts        :: T
    field     :: Symbol
    component :: Int

    function BulkTimeSeries(foldername::String, field::Symbol, c::Int)
        ts = OpenPMDTimeSeries(foldername, "bulk_")
        new{3,typeof(ts)}(ts, field, c)
    end
end

struct VEVTimeSeries{N,T} <: TimeSeries{N,T}
    ts  :: T
    vev :: Symbol

    function VEVTimeSeries(foldername::String, vev::Symbol)
        ts = OpenPMDTimeSeries(foldername, "boundary_")
        new{3,typeof(ts)}(ts, vev)
    end
end


function get_data(ff::BoundaryTimeSeries; it=Int, verbose::Bool=false)
    get_field(ff.ts, it=it, verbose=verbose, field=String(ff.field))
end

function get_data(xi::XiTimeSeries; it=Int, verbose::Bool=false)
    get_field(xi.ts, it=it, verbose=verbose, field="xi")
end

function get_data(ff::BulkTimeSeries; it=Int, verbose::Bool=false)
    field = "$(ff.field) c=$(ff.component)"
    get_field(ff.ts, it=it, verbose=verbose, field=field)
end

function get_data(ff::VEVTimeSeries; it=Int, verbose::Bool=false)
    if ff.vev == :energy
        return get_energy(ff.ts, it=it, verbose=verbose)
    elseif ff.vev == :Jx
        return get_Jx(ff.ts, it=it, verbose=verbose)
    elseif ff.vev == :Jy
        return get_Jy(ff.ts, it=it, verbose=verbose)
    elseif ff.vev == :px
        return get_px(ff.ts, it=it, verbose=verbose)
    elseif ff.vev == :py
        return get_py(ff.ts, it=it, verbose=verbose)
    elseif ff.vev == :pz
        return get_pz(ff.ts, it=it, verbose=verbose)
    elseif ff.vev == :Jx
        return get_Jx(ff.ts, it=it, verbose=verbose)
    elseif ff.vev == :pxy
        return get_pxy(ff.ts, it=it, verbose=verbose)
    elseif ff.vev == :Ophi
        return get_Ophi(ff.ts, it=it, verbose=verbose)
    else
        error("Unknown VEV")
    end
end


function Base.getindex(ff::TimeSeries, a::Int, idx::Vararg)
    it = ff.ts.iterations[a]
    f, = get_data(ff, it=it)
    f[idx...]
end

function Base.getindex(ff::TimeSeries, aa::UnitRange, idx::Vararg)
    it  = ff.ts.iterations[aa[1]]
    f0, = get_data(ff, it=it)

    Na    = length(aa)
    size_ = size(f0[idx...])
    f     = zeros(Na, size_...)

    slicer = [Colon() for _ in 1:ndims(f0)]
    for a in aa
        it  = ff.ts.iterations[a]
        f0, = get_data(ff, it=it)
        f[a,slicer...] .= f0[idx...]
    end
    f
end

function Base.getindex(ff::TimeSeries, ::Colon, idx::Vararg)
    it  = ff.ts.iterations[1]
    f0, = get_data(ff, it=it)

    Na    = length(ff.ts.iterations)
    size_ = size(f0[idx...])
    f     = zeros(Na, size_...)

    slicer = [Colon() for _ in 1:ndims(f0)]
    for a in 1:Na
        it  = ff.ts.iterations[a]
        f0, = get_data(ff, it=it)
        f[a,slicer...] .= f0[idx...]
    end
    f
end


function Jecco.get_coords(ff::TimeSeries, a::Int, idx::Vararg)
    it = ff.ts.iterations[a]
    f, chart = get_data(ff, it=it)
    t = ff.ts.current_t
    t, chart[idx...]...
end

function Jecco.get_coords(ff::TimeSeries, aa::UnitRange, idx::Vararg)
    it = ff.ts.iterations[aa[1]]
    f0, chart = get_data(ff, it=it)
    t0 = ff.ts.current_t

    Na = length(aa)
    t  = Vector{typeof(t0)}(undef, Na)
    for a in aa
        it   = ff.ts.iterations[a]
        get_data(ff, it=it)
        t0   = ff.ts.current_t
        t[a] = t0
    end
    t, chart[idx...]...
end

function Jecco.get_coords(ff::TimeSeries, ::Colon, idx::Vararg)
    it = ff.ts.iterations[1]
    f0, chart = get_data(ff, it=it)
    t0 = ff.ts.current_t

    Na = length(ff.ts.iterations)
    t  = Vector{typeof(t0)}(undef, Na)
    for a in 1:Na
        it = ff.ts.iterations[a]
        get_data(ff, it=it)
        t0 = ff.ts.current_t
        t[a] = t0
    end
    t, chart[idx...]...
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

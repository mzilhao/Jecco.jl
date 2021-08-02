using HDF5

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
    property :: Symbol

    function AHTimeSeries(foldername::String, property::Symbol)
        ts       = OpenPMDTimeSeries(foldername, "constrained_")
        ts_diag  = OpenPMDTimeSeries(foldername, "diagnostics_")
        new{OpenPMDTimeSeries}(ts, ts_diag, property)
    end
end

struct GWTimeSeries{T} <: TimeSeries{2,T}
    ts    :: T
    field :: Symbol

    function GWTimeSeries(foldername::String, field::Symbol)
        ts = OpenPMDTimeSeries(foldername, "perturbation_")
        new{typeof(ts)}(ts, field)
    end
end

struct TTTimeSeries{T} <: TimeSeries{2,T}
    ts    :: T
    field :: Symbol

    function TTTimeSeries(foldername::String, field::Symbol)
        ts = OpenPMDTimeSeries(foldername, "TT_")
        new{typeof(ts)}(ts, field)
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

function get_data(ff::GWTimeSeries, it::Int)
    f, chart = get_field(ff.ts, it=it, field=String(ff.field))
    x, y  = chart[:]
    f[:,:], [x, y]
end

function get_data(ff::TTTimeSeries, it::Int)
    f, chart = get_field(ff.ts, it=it, field=String(ff.field))
    x, y  = chart[:]
    f[:,:], [x, y]
end

function Base.size(ff::TimeSeries)
    Nt  = length(ff.ts.iterations)
    it0 = ff.ts.iterations[1]
    f0, = get_data(ff, it0)
    size_ = size(f0)
    Nt, size_...
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

function e_to_mathematica(ts::OpenPMDTimeSeries, group::HDF5.Group, dit::Int)
    iterations = ts.iterations

    Nt = length(iterations)
    t  = zeros(Nt)

    it = 0
    _, chart = get_energy(ts, it=it)

    _, x, y = chart[:]
    Nx = length(x)
    Ny = length(y)

    en   = zeros(Nt,Nx,Ny)

    for (idx,it) in enumerate(iterations)
        en[idx,:,:] = get_energy(ts, it=it)[1][1,:,:]
        t[idx]      = ts.current_t
    end

    T_m = zeros(4, Int(floor(Nt/dit))*Nx*Ny)
    n = 1
    @fastmath @inbounds for i in 1:dit:Nt-Nt%dit
        for j in 1:Nx
            for k in 1:Ny
                T_m[:,n] = [t[i] x[j] y[k] en[i,j,k]]
                n += 1
            end
        end
    end
    Jecco.write_dataset(group, "energy", T_m)
    nothing
end

function Jx_to_mathematica(ts::OpenPMDTimeSeries, group::HDF5.Group, dit::Int)
    iterations = ts.iterations

    Nt = length(iterations)
    t  = zeros(Nt)

    it = 0
    _, chart = get_Jx(ts, it=it)

    _, x, y = chart[:]
    Nx = length(x)
    Ny = length(y)

    Jx   = zeros(Nt,Nx,Ny)

    for (idx,it) in enumerate(iterations)
        Jx[idx,:,:] = get_Jx(ts, it=it)[1][1,:,:]
        t[idx]      = ts.current_t
    end

    T_m = zeros(4, Int(floor(Nt/dit))*Nx*Ny)
    n = 1
    @fastmath @inbounds for i in 1:dit:Nt-Nt%dit
        for j in 1:Nx
            for k in 1:Ny
                T_m[:,n] = [t[i] x[j] y[k] Jx[i,j,k]]
                n += 1
            end
        end
    end
    Jecco.write_dataset(group, "Jx", T_m)
    nothing
end

function Jy_to_mathematica(ts::OpenPMDTimeSeries, group::HDF5.Group, dit::Int)
    iterations = ts.iterations

    Nt = length(iterations)
    t  = zeros(Nt)

    it = 0
    _, chart = get_Jy(ts, it=it)

    _, x, y = chart[:]
    Nx = length(x)
    Ny = length(y)

    Jy   = zeros(Nt,Nx,Ny)

    for (idx,it) in enumerate(iterations)
        Jy[idx,:,:] = get_Jy(ts, it=it)[1][1,:,:]
        t[idx]      = ts.current_t
    end

    T_m = zeros(4, Int(floor(Nt/dit))*Nx*Ny)
    n = 1
    @fastmath @inbounds for i in 1:dit:Nt-Nt%dit
        for j in 1:Nx
            for k in 1:Ny
                T_m[:,n] = [t[i] x[j] y[k] Jy[i,j,k]]
                n += 1
            end
        end
    end
    Jecco.write_dataset(group, "Jy", T_m)
    nothing
end

function px_to_mathematica(ts::OpenPMDTimeSeries, group::HDF5.Group, dit::Int)
    iterations = ts.iterations

    Nt = length(iterations)
    t  = zeros(Nt)

    it = 0
    _, chart = get_px(ts, it=it)

    _, x, y = chart[:]
    Nx = length(x)
    Ny = length(y)

    px   = zeros(Nt,Nx,Ny)

    for (idx,it) in enumerate(iterations)
        px[idx,:,:] = get_px(ts, it=it)[1][1,:,:]
        t[idx]      = ts.current_t
    end

    T_m = zeros(4, Int(floor(Nt/dit))*Nx*Ny)
    n = 1
    @fastmath @inbounds for i in 1:dit:Nt-Nt%dit
        for j in 1:Nx
            for k in 1:Ny
                T_m[:,n] = [t[i] x[j] y[k] px[i,j,k]]
                n += 1
            end
        end
    end
    Jecco.write_dataset(group, "px", T_m)
    nothing
end

function pxy_to_mathematica(ts::OpenPMDTimeSeries, group::HDF5.Group, dit::Int)
    iterations = ts.iterations

    Nt = length(iterations)
    t  = zeros(Nt)

    it = 0
    _, chart = get_pxy(ts, it=it)

    _, x, y = chart[:]
    Nx = length(x)
    Ny = length(y)

    pxy   = zeros(Nt,Nx,Ny)

    for (idx,it) in enumerate(iterations)
        pxy[idx,:,:] = get_pxy(ts, it=it)[1][1,:,:]
        t[idx]      = ts.current_t
    end

    T_m = zeros(4, Int(floor(Nt/dit))*Nx*Ny)
    n = 1
    @fastmath @inbounds for i in 1:dit:Nt-Nt%dit
        for j in 1:Nx
            for k in 1:Ny
                T_m[:,n] = [t[i] x[j] y[k] pxy[i,j,k]]
                n += 1
            end
        end
    end
    Jecco.write_dataset(group, "pxy", T_m)
    nothing
end

function py_to_mathematica(ts::OpenPMDTimeSeries, group::HDF5.Group, dit::Int)
    iterations = ts.iterations

    Nt = length(iterations)
    t  = zeros(Nt)

    it = 0
    _, chart = get_py(ts, it=it)

    _, x, y = chart[:]
    Nx = length(x)
    Ny = length(y)

    py   = zeros(Nt,Nx,Ny)

    for (idx,it) in enumerate(iterations)
        py[idx,:,:] = get_py(ts, it=it)[1][1,:,:]
        t[idx]      = ts.current_t
    end

    T_m = zeros(4, Int(floor(Nt/dit))*Nx*Ny)
    n = 1
    @fastmath @inbounds for i in 1:dit:Nt-Nt%dit
        for j in 1:Nx
            for k in 1:Ny
                T_m[:,n] = [t[i] x[j] y[k] py[i,j,k]]
                n += 1
            end
        end
    end
    Jecco.write_dataset(group, "py", T_m)
    nothing
end

function pz_to_mathematica(ts::OpenPMDTimeSeries, group::HDF5.Group, dit::Int)
    iterations = ts.iterations

    Nt = length(iterations)
    t  = zeros(Nt)

    it = 0
    _, chart = get_pz(ts, it=it)

    _, x, y = chart[:]
    Nx = length(x)
    Ny = length(y)

    pz   = zeros(Nt,Nx,Ny)

    for (idx,it) in enumerate(iterations)
        pz[idx,:,:] = get_pz(ts, it=it)[1][1,:,:]
        t[idx]      = ts.current_t
    end

    T_m = zeros(4, Int(floor(Nt/dit))*Nx*Ny)
    n = 1
    @fastmath @inbounds for i in 1:dit:Nt-Nt%dit
        for j in 1:Nx
            for k in 1:Ny
                T_m[:,n] = [t[i] x[j] y[k] pz[i,j,k]]
                n += 1
            end
        end
    end
    Jecco.write_dataset(group, "pz", T_m)
    nothing
end

function convert_to_mathematica_2(dirname::String; dit::Int = 1, outfile::String="data_mathematica.h5")
    output     = abspath(dirname, outfile)
    fid        = h5open(output, "w")
    group_st   = g_create(fid, "data")
    ts         = OpenPMDTimeSeries(dirname, prefix="boundary_")

    e_to_mathematica(ts, group_st, dit)
    Jx_to_mathematica(ts, group_st, dit)
    Jy_to_mathematica(ts, group_st, dit)
    px_to_mathematica(ts, group_st, dit)
    pxy_to_mathematica(ts, group_st, dit)
    py_to_mathematica(ts, group_st, dit)
    pz_to_mathematica(ts, group_st, dit)

    close(fid)
end

#If you have memory issues for big datasets run this routine for the different VEVs you want. it does not load that much.
function Energy_to_mathematica(dirname::String; dit::Int = 1, outfile::String="energy_mathematica.h5")
    ts = OpenPMDTimeSeries(dirname, prefix="boundary_")

    iterations = ts.iterations

    Nt = length(iterations)
    t  = zeros(Nt)

    it = 0
    _, chart = get_energy(ts, it=it)

    _, x, y = chart[:]
    Nx = length(x)
    Ny = length(y)

    en   = zeros(Nt,Nx,Ny)

    for (idx,it) in enumerate(iterations)
        en[idx,:,:] = get_energy(ts, it=it)[1][1,:,:]
        t[idx]      = ts.current_t
    end

    # store in an array suitable for Mathematica
    T_m = zeros(4, Int(floor(Nt/dit))*Nx*Ny)

    n = 1
    @fastmath @inbounds for i in 1:dit:Nt-Nt%dit
        for j in 1:Nx
            for k in 1:Ny
                T_m[:,n] = [t[i] x[j] y[k] en[i,j,k]]
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

function GW_to_mathematica(dirname::String; dit::Int = 1, outfile::String="GW_mathematica.h5")

    ts         = OpenPMDTimeSeries(dirname, prefix="perturbation_")
    iterations = ts.iterations
    Nt         = length(iterations)
    t          = zeros(Nt)
    it         = 0
    _, chart   = get_field(ts, it=it, field="hxx")
    x, y       = chart[:]
    Nx, Ny     = size(chart)
    hdot2      = zeros(Nt, Nx, Ny)

    for (idx,it) in enumerate(iterations)
        hdxx           = get_field(ts, it=it, field="hdxx")[1]
        hdxy           = get_field(ts, it=it, field="hdxy")[1]
        hdyy           = get_field(ts, it=it, field="hdyy")[1]
        hdzz           = get_field(ts, it=it, field="hdzz")[1]
        t[idx]         = ts.current_t
        hdot2[idx,:,:] = hdxx.^2+hdxy.^2+hdyy.^2+hdzz.^2
    end

    # store in an array suitable for Mathematica
    hdot2_m = zeros(4, Int(floor(Nt/dit))*Nx*Ny)

    n = 1
    @fastmath @inbounds for i in 1:dit:Nt-Nt%dit
        for j in 1:Nx
            for k in 1:Ny
                hdot2_m[:,n] = [t[i] x[j] y[k] hdot2[i,j,k]]
                n += 1
            end
        end
    end

    output   = abspath(dirname, outfile)
    fid      = h5open(output, "w")
    group_st = g_create(fid, "data")
    Jecco.write_dataset(group_st, "hdot2", hdot2_m)
    close(fid)
end

function TT_to_mathematica(dirname::String; dit::Int = 1, outfile::String="TTT_mathematica.h5")

    ts         = OpenPMDTimeSeries(dirname, prefix="TT_")
    iterations = ts.iterations
    Nt         = length(iterations)
    t          = zeros(Nt)
    it         = 0
    _, chart   = get_field(ts, it=it, field="Txx")
    x, y       = chart[:]
    Nx, Ny     = size(chart)
    T2         = zeros(Nt, Nx, Ny)

    for (idx,it) in enumerate(iterations)
        Txx         = get_field(ts, it=it, field="Txx")[1]
        Txy         = get_field(ts, it=it, field="Txy")[1]
        Tyy         = get_field(ts, it=it, field="Tyy")[1]
        Tzz         = get_field(ts, it=it, field="Tzz")[1]
        t[idx]      = ts.current_t
        T2[idx,:,:] = Txx.^2+Txy.^2+Tyy.^2+Tzz.^2
    end

    # store in an array suitable for Mathematica
    T2_m = zeros(4, Int(floor(Nt/dit))*Nx*Ny)

    n = 1
    @fastmath @inbounds for i in 1:dit:Nt-Nt%dit
        for j in 1:Nx
            for k in 1:Ny
                T2_m[:,n] = [t[i] x[j] y[k] T2[i,j,k]]
                n += 1
            end
        end
    end

    output   = abspath(dirname, outfile)
    fid      = h5open(output, "w")
    group_st = g_create(fid, "data")
    Jecco.write_dataset(group_st, "T2", T2_m)
    close(fid)
end

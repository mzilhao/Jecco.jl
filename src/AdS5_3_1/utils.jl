using HDF5
using Interpolations

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

function get_data(ff::GWTimeSeries, it::Int)
    if String(ff.field) == "hd2"
        hdxx, chart = get_field(ff.ts, it=it, field="hdxx")
        hdxy, _     = get_field(ff.ts, it=it, field="hdxy")
        hdyy, _     = get_field(ff.ts, it=it, field="hdyy")
        hdzz, _     = get_field(ff.ts, it=it, field="hdzz")

        f   = hdxx.^2+hdxy.^2+hdyy.^2+hdzz.^2
    else
        f, chart = get_field(ff.ts, it=it, field=String(ff.field))
    end
    x, y  = chart[:]

    f[:,:], [x, y]
end

function get_data(ff::TTTimeSeries, it::Int)
    if String(ff.field) == "T2"
        Txx, chart = get_field(ff.ts, it=it, field="Txx")
        Txy, _     = get_field(ff.ts, it=it, field="Txy")
        Tyy, _     = get_field(ff.ts, it=it, field="Tyy")
        Tzz, _     = get_field(ff.ts, it=it, field="Tzz")

        f   = Txx.^2+Txy.^2+Tyy.^2+Tzz.^2
    elseif String(ff.field)[2] == 'd'
        field = split(String(ff.field), 'd')[1]*split(String(ff.field), 'd')[2]
        if it == ff.ts.iterations[1]
            T0, chart = get_field(ff.ts, it=it, field=field)
            i         = ff.ts.current_i
            t0        = ff.ts.current_t
            it1       = ff.ts.iterations[i+1]
            it2       = ff.ts.iterations[i+2]
            T1, _     = get_field(ff.ts, it=it1, field=field)
            t1        = ff.ts.current_t
            T2, _     = get_field(ff.ts, it=it2, field=field)
            f         = (-1.5 .*T0+2 .*T1-0.5 .*T2)./(t1-t0)
        elseif it == ff.ts.iterations[end]
            T0, chart = get_field(ff.ts, it=it, field=field)
            i         = ff.ts.current_i
            t0        = ff.ts.current_t
            it1       = ff.ts.iterations[i-1]
            it2       = ff.ts.iterations[i-2]
            T1, _     = get_field(ff.ts, it=it1, field=field)
            t1        = ff.ts.current_t
            T2, _     = get_field(ff.ts, it=it2, field=field)
            f         = (1.5 .*T0-2 .*T1+0.5 .*T2)./(t0-t1)
        else
            i         = findfirst(ff.ts.iterations .== it)
            it1       = ff.ts.iterations[i-1]
            it2       = ff.ts.iterations[i+1]
            T1, chart = get_field(ff.ts, it=it1, field=field)
            t1        = ff.ts.current_t
            T2, _     = get_field(ff.ts, it=it2, field=field)
            t2        = ff.ts.current_t
            f         = (T2-T1)./(t2-t1)
        end
    else
        f, chart = get_field(ff.ts, it=it, field=String(ff.field))
    end
    x, y  = chart[:]

    f[:,:], [x, y]
end

function get_Td2(Tdxx::T, Tdxy::T, Tdyy::T, Tdzz::T, it::Int) where {T<:TTTimeSeries}
    Tddxx, _ = get_data(Tdxx, it)
    Tddxy, _ = get_data(Tdxy, it)
    Tddyy, _ = get_data(Tdyy, it)
    Tddzz, _ = get_data(Tdzz, it)

    get_Td2(Tddxx, Tddxy, Tddyy, Tddzz)
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

function GW_to_mathematica(dirname::String; dit::Int = 1, outfile::String="GW_mathematica.h5")

    ts         = OpenPMDTimeSeries(dirname, prefix="perturbation_")
    iterations = ts.iterations
    Nt         = length(iterations)
    t          = zeros(Nt)
    it         = 0
    _, chart   = get_field(ts, it=it, field="hxx")
    x, y       = chart[:]
    Nx, Ny     = size(chart)
    # hdot2      = zeros(Nt, Nx, Ny)
    hdxx      = zeros(Nt, Nx, Ny)
    hdxy      = zeros(Nt, Nx, Ny)
    hdyy      = zeros(Nt, Nx, Ny)
    hdzz      = zeros(Nt, Nx, Ny)
    for (idx,it) in enumerate(iterations)
        hdxx[idx,:,:]  = get_field(ts, it=it, field="hdxx")[1]
        hdxy[idx,:,:]  = get_field(ts, it=it, field="hdxy")[1]
        hdyy[idx,:,:]  = get_field(ts, it=it, field="hdyy")[1]
        hdzz[idx,:,:]  = get_field(ts, it=it, field="hdzz")[1]
        t[idx]         = ts.current_t
        # hdot2[idx,:,:] = hdxx.^2+hdxy.^2+hdyy.^2+hdzz.^2
    end

    # store in an array suitable for Mathematica
    hdij = zeros(7, Int(floor(Nt/dit))*Nx*Ny)

    n = 1
    @inbounds for i in 1:dit:Nt-Nt%dit
        for j in 1:Nx
            for k in 1:Ny
                hdij[:,n] = [t[i] x[j] y[k] hdxx[i,j,k] hdxy[i,j,k] hdyy[i,j,k] hdzz[i,j,k]]
                n += 1
            end
        end
    end

    output   = abspath(dirname, outfile)
    fid      = h5open(output, "w")
    group_st = g_create(fid, "data")
    Jecco.write_dataset(group_st, "hdij", hdij)
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

function Pk_to_mathematica(dirname::String; outfile::String="Pk_mathematica.h5")

    ts         = OpenPMDTimeSeries(dirname, prefix="boundary_")
    iterations = ts.iterations
    Nt         = length(iterations)
    t          = zeros(Nt)
    aux, chart = get_energy(ts, it=iterations[1])
    _, δx, δy     = Jecco.delta(chart)
    _, Nx, Ny  = size(aux)
    plan       = 1/(Nx*Ny)*plan_rfft(aux[1,:,:])
    kx         = 2*π.*rfftfreq(Nx, 1/δx)
    ky         = 2*π.*fftfreq(Ny, 1/δy)
    Nkx, Nky   = (length(kx),length(ky))
    px         = im.*zeros(Nt,Nkx,Nky)
    pxy        = im.*zeros(Nt,Nkx,Nky)
    py         = im.*zeros(Nt,Nkx,Nky)
    pz         = im.*zeros(Nt,Nkx,Nky)
    for (idx,it) in enumerate(iterations)
        px[idx,:,:]  .= (plan * get_px(ts, it=it)[1][1,:,:])
        pxy[idx,:,:] .= (plan * get_pxy(ts, it=it)[1][1,:,:])
        py[idx,:,:]  .= (plan * get_py(ts, it=it)[1][1,:,:])
        pz[idx,:,:]  .= (plan * get_pz(ts, it=it)[1][1,:,:])
        t[idx]        = ts.current_t
    end

    # store in an array suitable for Mathematica
    T_m = zeros(11, Nt*Nkx*Nky)

    n = 1
    @fastmath @inbounds for i in 1:Nt
        for j in 1:Nkx
            for k in 1:Nky
                T_m[:,n] = [t[i] kx[j] ky[k] real(px[i,j,k]) imag(px[i,j,k]) real(pxy[i,j,k]) imag(pxy[i,j,k]) real(py[i,j,k]) imag(py[i,j,k]) real(pz[i,j,k]) imag(pz[i,j,k])]
                n += 1
            end
        end
    end

    output   = abspath(dirname, outfile)
    fid      = h5open(output, "w")
    group_st = g_create(fid, "data")
    Jecco.write_dataset(group_st, "Pk", T_m)
    close(fid)
end

function IdealPk_to_mathematica(dirname::String, dir_eos::String; outfile::String="Pk_Ideal_mathematica.h5")
    pxIdeal    = IdealHydroTimeSeries(dirname, dir_eos, :px)
    pxyIdeal   = IdealHydroTimeSeries(dirname, dir_eos, :pxy)
    pyIdeal    = IdealHydroTimeSeries(dirname, dir_eos, :py)
    pzIdeal    = IdealHydroTimeSeries(dirname, dir_eos, :pz)
    ts         = pxIdeal.ts
    iterations = ts.iterations
    Nt         = length(iterations)
    t          = zeros(Nt)
    aux, chart = get_energy(ts, it=iterations[1])
    _, δx, δy  = Jecco.delta(chart)
    _, Nx, Ny  = size(aux)
    plan       = 1/(Nx*Ny)*plan_rfft(aux[1,:,:])
    kx         = 2*π.*rfftfreq(Nx, 1/δx)
    ky         = 2*π.*fftfreq(Ny, 1/δy)
    Nkx, Nky   = (length(kx),length(ky))
    px         = im.*zeros(Nt,Nkx,Nky)
    pxy        = im.*zeros(Nt,Nkx,Nky)
    py         = im.*zeros(Nt,Nkx,Nky)
    pz         = im.*zeros(Nt,Nkx,Nky)
    for (idx,it) in enumerate(iterations)
        px[idx,:,:]  = (plan * pxIdeal[idx,:,:])
        pxy[idx,:,:] = (plan * pxyIdeal[idx,:,:])
        py[idx,:,:]  = (plan * pyIdeal[idx,:,:])
        pz[idx,:,:]  = (plan * pzIdeal[idx,:,:])
        t[idx]        = ts.current_t
    end

    # store in an array suitable for Mathematica
    T_m = zeros(11, Nt*Nkx*Nky)

    n = 1
    @inbounds for i in 1:Nt
        for j in 1:Nkx
            for k in 1:Nky
                T_m[:,n] = [t[i] kx[j] ky[k] real(px[i,j,k]) imag(px[i,j,k]) real(pxy[i,j,k]) imag(pxy[i,j,k]) real(py[i,j,k]) imag(py[i,j,k]) real(pz[i,j,k]) imag(pz[i,j,k])]
                n += 1
            end
        end
    end

    output   = abspath(dirname, outfile)
    fid      = h5open(output, "w")
    group_st = g_create(fid, "data")
    Jecco.write_dataset(group_st, "PIdealk", T_m)
    close(fid)
end

function Tk_to_mathematica(dirname::String; outfile::String="Tk_mathematica.h5")

    ts         = OpenPMDTimeSeries(dirname, prefix="TT_")
    iterations = ts.iterations
    Nt         = length(iterations)
    t          = zeros(Nt)
    aux, chart = get_field(ts, it=iterations[1], field="Txx")
    δx, δy     = Jecco.delta(chart)
    Nx, Ny     = size(aux)
    plan       = plan_rfft(aux)
    kx         = 2*π.*rfftfreq(Nx, 1/δx)
    ky         = 2*π.*fftfreq(Ny, 1/δy)
    Nkx, Nky   = (length(kx),length(ky))
    Txx        = zeros(Nt,Nkx,Nky)
    Txy        = zeros(Nt,Nkx,Nky)
    Tyy        = zeros(Nt,Nkx,Nky)
    Tzz        = zeros(Nt,Nkx,Nky)

    it         = 0
    for (idx,it) in enumerate(iterations)
        Txx[idx,:,:]  = 1/(Nx*Ny).*abs.(plan * get_field(ts, it=it, field="Txx")[1])
        Txy[idx,:,:]  = 1/(Nx*Ny).*abs.(plan * get_field(ts, it=it, field="Txy")[1])
        Tyy[idx,:,:]  = 1/(Nx*Ny).*abs.(plan * get_field(ts, it=it, field="Tyy")[1])
        Tzz[idx,:,:]  = 1/(Nx*Ny).*abs.(plan * get_field(ts, it=it, field="Tzz")[1])
        t[idx]        = ts.current_t
    end

    # store in an array suitable for Mathematica
    T_m = zeros(7, Nt*Nkx*Nky)

    n = 1
    @fastmath @inbounds for i in 1:Nt
        for j in 1:Nkx
            for k in 1:Nky
                if j==1 && k==1
                    T_m[:,n] = [t[i] kx[j] ky[k] 0. 0. 0. 0.]
                else
                    T_m[:,n] = [t[i] kx[j] ky[k] Txx[i,j,k] Txy[i,j,k] Tyy[i,j,k] Tzz[i,j,k]]
                end
                n += 1
            end
        end
    end

    output   = abspath(dirname, outfile)
    fid      = h5open(output, "w")
    group_st = g_create(fid, "data")
    Jecco.write_dataset(group_st, "Tk", T_m)
    close(fid)
end

function hdk_to_mathematica(dirname::String; outfile::String="hdk_mathematica.h5")

    ts         = OpenPMDTimeSeries(dirname, prefix="perturbation_")
    iterations = ts.iterations
    Nt         = length(iterations)
    t          = zeros(Nt)
    aux, chart = get_field(ts, it=iterations[1], field="hdxx")
    δx, δy     = Jecco.delta(chart)
    Nx, Ny     = size(aux)
    plan       = plan_rfft(aux)
    kx         = 2*π.*rfftfreq(Nx, 1/δx)
    ky         = 2*π.*fftfreq(Ny, 1/δy)
    # println(kx)
    Nkx, Nky   = (length(kx),length(ky))
    hdxx       = zeros(Nt,Nkx,Nky)
    hdxy       = zeros(Nt,Nkx,Nky)
    hdyy       = zeros(Nt,Nkx,Nky)
    hdzz       = zeros(Nt,Nkx,Nky)

    it         = 0
    for (idx,it) in enumerate(iterations)
        hdxx[idx,:,:]  = 1/(Nx*Ny).*abs.(plan * get_field(ts, it=it, field="hdxx")[1])
        hdxy[idx,:,:]  = 1/(Nx*Ny).*abs.(plan * get_field(ts, it=it, field="hdxy")[1])
        hdyy[idx,:,:]  = 1/(Nx*Ny).*abs.(plan * get_field(ts, it=it, field="hdyy")[1])
        hdzz[idx,:,:]  = 1/(Nx*Ny).*abs.(plan * get_field(ts, it=it, field="hdzz")[1])
        t[idx]         = ts.current_t
    end

    # store in an array suitable for Mathematica
    h_m = zeros(8, Nt*Nkx*Nky)

    n = 1
    @fastmath @inbounds for i in 1:Nt
        for j in 1:Nkx
            for k in 1:Nky
                if j==1 && k==1
                    h_m[:,n] = [t[i] kx[j] ky[k] 0. 0. 0. 0. 0.]
                else
                    hd2      = hdxx[i,j,k]^2+hdxy[i,j,k]^2+hdyy[i,j,k]^2+hdzz[i,j,k]^2
                    h_m[:,n] = [t[i] kx[j] ky[k] hdxx[i,j,k] hdxy[i,j,k] hdyy[i,j,k] hdzz[i,j,k] hd2]
                end
                n += 1
            end
        end
    end
    output   = abspath(dirname, outfile)
    fid      = h5open(output, "w")
    group_st = g_create(fid, "data")
    Jecco.write_dataset(group_st, "hd", h_m)
    close(fid)
end

function drhodlogk_to_mathematica_1(dirname::String, δk::Real; outfile::String="drhoGWdlogk_mathematica.h5")

    ts         = OpenPMDTimeSeries(dirname, prefix="perturbation_")
    iterations = ts.iterations
    Nt         = length(iterations)
    t          = zeros(Nt)
    aux, chart = get_field(ts, it=iterations[1], field="hdxx")
    δx, δy     = Jecco.delta(chart)
    Nx, Ny     = size(aux)
    plan       = 1/(Nx*Ny)*plan_rfft(aux)
    # plan       = δx*δy*plan_rfft(aux)
    kx         = sort(2*π.*rfftfreq(Nx, 1/δx)[1:end-1])
    ky         = sort(2*π.*fftfreq(Ny, 1/δy)[1:end-1])
    δkx        = kx[2]-kx[1]
    δky        = ky[2]-ky[1]
    Nkx, Nky   = (length(kx),length(ky))
    k          = Array(0:δk:kx[end])
    Nk         = length(k)
    hdxx       = zeros(Nkx,Nky)
    hdxy       = zeros(Nkx,Nky)
    hdyy       = zeros(Nkx,Nky)
    hdzz       = zeros(Nkx,Nky)
    dρdlogk    = zeros(Nt,Nk)
    ρ          = zeros(Nt)

    for (idx,it) in enumerate(iterations)
        hdxx[:,1:Nkx-1] = abs.(plan * get_field(ts, it=it, field="hdxx")[1])[1:end-1,Nkx+2:end]
        hdxx[:,Nkx:end] = abs.(plan * get_field(ts, it=it, field="hdxx")[1])[1:end-1,1:Nkx]
        hdxy[:,1:Nkx-1] = abs.(plan * get_field(ts, it=it, field="hdxy")[1])[1:end-1,Nkx+2:end]
        hdxy[:,Nkx:end] = abs.(plan * get_field(ts, it=it, field="hdxy")[1])[1:end-1,1:Nkx]
        hdyy[:,1:Nkx-1] = abs.(plan * get_field(ts, it=it, field="hdyy")[1])[1:end-1,Nkx+2:end]
        hdyy[:,Nkx:end] = abs.(plan * get_field(ts, it=it, field="hdyy")[1])[1:end-1,1:Nkx]
        hdzz[:,1:Nkx-1] = abs.(plan * get_field(ts, it=it, field="hdzz")[1])[1:end-1,Nkx+2:end]
        hdzz[:,Nkx:end] = abs.(plan * get_field(ts, it=it, field="hdzz")[1])[1:end-1,1:Nkx]
        dρdkxdky        = 2/π*(hdxx.^2+hdxy.^2+hdyy.^2+hdzz.^2)
        t[idx]          = ts.current_t
        for j in 1:Nky
            for i in 1:Nkx
                if kx[i] == 0.0 && ky[j] == 0.0
                    nothing
                else
                    ρ[idx] += 2*dρdkxdky[i,j]*δkx*δky
                end
            end
        end
        for n in 1:Nk
            kk = k[n]
            for j in 1:Nky
                for i in 1:Nkx
                    if kk-δk/2 <= sqrt(kx[i]^2+ky[j]^2) < kk+δk/2
                        dρdlogk[idx,n] += 2*kk*dρdkxdky[i,j]*δkx*δky/δk
                    end
                end
            end
            # GW_m[:,n,idx] = [t[idx] k[n] dρdlogk[n]]
        end
    end
    # store in an array suitable for Mathematica
    # hd2_m = zeros(3, Nky, Nkx, Nt)
    GW_m  = zeros(3, Nk, Nt)
    ρ_m   = zeros(2, Nt)
    # n = 1
    # @inbounds for i in 1:Nt
    #     for j in 1:Nkx
    #         for k in 1:Nky
    #             hd2_m[:,k,j,i] = [t[i] kx[j] ky[k] dρdkxdky[]]
    #         end
    #     end
    # end
    @fastmath @inbounds for i in 1:Nt
        for j in 1:Nk
            GW_m[:,j,i] = [t[i] k[j] dρdlogk[i,j]]
        end
        ρ_m[:,i] = [t[i] ρ[i]]
    end

    output   = abspath(dirname, outfile)
    fid      = h5open(output, "w")
    group_st = g_create(fid, "data")
    Jecco.write_dataset(group_st, "drhodlogk", GW_m)
    Jecco.write_dataset(group_st, "rho", ρ_m)
    close(fid)
end

function ρGW(dirname::String; outfile::String="rhoGW_mathematica.h5")
    hd2        = GWTimeSeries(dirname, :hd2)
    t, _, _    = get_coords(hd2, :, :, :)
    Nt, Nx, Ny = size(hd2)
    plan       = 1/(Nx*Ny) * plan_rfft(hd2[1,:,:])
    ρ_m        = zeros(2, Nt)
    @inbounds for n in 1:Nt
        f        = (plan * hd2[n,:,:])[1,1]
        ρ_m[:,n] = [t[n] 8*π*f]
    end
    output   = abspath(dirname, outfile)
    fid      = h5open(output, "w")
    group_st = g_create(fid, "data")
    Jecco.write_dataset(group_st, "rho", ρ_m)
    close(fid)
end

# function drhodlogk_to_mathematica_2(dirname::String, δk::Real; outfile::String="drhoGWdlogk_mathematica_2.h5")
#
#     ts         = OpenPMDTimeSeries(dirname, prefix="perturbation_")
#     iterations = ts.iterations
#     Nt         = length(iterations)
#     t          = zeros(Nt)
#     aux, chart = get_field(ts, it=iterations[1], field="hdxx")
#     δx, δy     = Jecco.delta(chart)
#     Nx, Ny     = size(aux)
#     plan       = 1/(Nx*Ny)*plan_rfft(aux)
#     kx         = sort(2*π.*rfftfreq(Nx, 1/δx)[1:end-1])
#     ky         = sort(2*π.*fftfreq(Ny, 1/δy)[1:end-1])
#     Nkx, Nky   = (length(kx),length(ky))
#     k          = Array(0:δk:kx[end])
#     Nk         = length(k)
#     hdxx       = zeros(Nkx,Nky)
#     hdxy       = zeros(Nkx,Nky)
#     hdyy       = zeros(Nkx,Nky)
#     hdzz       = zeros(Nkx,Nky)
#     dρdlogk    = zeros(Nt,Nk)
#
#     for (idx,it) in enumerate(iterations)
#         hdxx[:,1:Nkx-1] = abs.(plan * get_field(ts, it=it, field="hdxx")[1])[1:end-1,Nkx+2:end]
#         hdxx[:,Nkx:end] = abs.(plan * get_field(ts, it=it, field="hdxx")[1])[1:end-1,1:Nkx]
#         hdxy[:,1:Nkx-1] = abs.(plan * get_field(ts, it=it, field="hdxy")[1])[1:end-1,Nkx+2:end]
#         hdxy[:,Nkx:end] = abs.(plan * get_field(ts, it=it, field="hdxy")[1])[1:end-1,1:Nkx]
#         hdyy[:,1:Nkx-1] = abs.(plan * get_field(ts, it=it, field="hdyy")[1])[1:end-1,Nkx+2:end]
#         hdyy[:,Nkx:end] = abs.(plan * get_field(ts, it=it, field="hdyy")[1])[1:end-1,1:Nkx]
#         hdzz[:,1:Nkx-1] = abs.(plan * get_field(ts, it=it, field="hdzz")[1])[1:end-1,Nkx+2:end]
#         hdzz[:,Nkx:end] = abs.(plan * get_field(ts, it=it, field="hdzz")[1])[1:end-1,1:Nkx]
#         dρdkxdky        = 2/(2*π)^2 .*(hdxx.^2+hdxy.^2+hdyy.^2+hdzz.^2)
#         t[idx]          = ts.current_t
#         for n in 1:Nk
#             kk = k[n]
#             for j in 1:Nky
#                 for i in 1:Nkx
#                     if log(kk-δk) <= sqrt(kx[i]^2+ky[j]^2) <= log(kk+δk)
#                         dρdlogk[idx,n] +=  dρdkxdky[i,j]
#                     end
#                 end
#             end
#         end
#     end
#     # store in an array suitable for Mathematica
#     GW_m = zeros(3, Nk, Nt)
#     # n = 1
#     @fastmath @inbounds for i in 1:Nt
#         for j in 1:Nk
#             GW_m[:,j,i] = [t[i] k[j] dρdlogk[i,j]]
#             # n += 1
#         end
#     end
#
#     output   = abspath(dirname, outfile)
#     fid      = h5open(output, "w")
#     group_st = g_create(fid, "data")
#     Jecco.write_dataset(group_st, "drhodlogk", GW_m)
#     close(fid)
# end



#From exponential basis to cos sin one for real functions.
function Exp_to_cos_sin(f::Array{T,2}) where {T<:Complex}
    nkx, nky = size(f)
    Nkx      = nkx-1
    Nky      = Int(floor(nky/2))
    a        = zeros(Nkx, Nky)
    b        = zeros(Nkx, Nky)
    c        = zeros(Nkx, Nky)
    d        = zeros(Nkx, Nky)
    # fk       = im.*zeros(Nkx, 2*Nky)

    # fk[:,1:Nky] = f[2:end-1,2:Nky+1]
    # for j in 1:Nky
    #     fk[:,Nky+j] = f[2:end-1,end-j+1]
    # end
    # a = real.(f[:,1:Nky]+f[:,Nky+1:2*Nky])
    # b = imag.(-f[:,1:Nky]+f[:,Nky+1:2*Nky])
    # c = imag.(-f[:,1:Nky]-f[:,Nky+1:2*Nky])
    # d = real.(-f[:,1:Nky]+f[:,Nky+1:2*Nky])
    # The (-1)^(i+j) factor compensate as the FFTW does the transform from the (-L/2,-L/2) while our modes where centered in the box.
    # It does not really matter... the abs value will remain unchanged
    @inbounds for j in 1:Nky
        for i in 1:Nkx
            if j == 1
                a[i,1] = (-1)^(i+j)*2*real(f[i,1])
                b[i,1] = 0.
                c[i,1] = -(-1)^(i+j)*2*imag(f[i,1])
                d[i,1] = 0.
            elseif i == 1
                a[i,j] = (-1)^(i+j)*real(f[i,j]+f[i,end-(j-2)])
                b[i,j] = (-1)^(i+j)*imag(-f[i,j]+f[i,end-(j-2)])
                c[i,j] = -(-1)^(i+j)*imag(f[i,j]+f[i,end-(j-2)])
                d[i,j] = -(-1)^(i+j)*real(-f[i,j]+f[i,end-(j-2)])
            else
                a[i,j] = (-1)^(i+j)*2*real(f[i,j]+f[i,end-(j-2)])
                b[i,j] = (-1)^(i+j)*2*imag(-f[i,j]+f[i,end-(j-2)])
                c[i,j] = -(-1)^(i+j)*2*imag(f[i,j]+f[i,end-(j-2)])
                d[i,j] = (-1)^(i+j)*2*real(-f[i,j]+f[i,end-(j-2)])
            end
        end
    end

    a[1,1] = real(f[1,1])

    a, b, c, d
end

#Fourier decomposition in coscos, cossin, sincos and sinsin basis
function Fourier_cos_sin(f::Array{T,2}) where {T<:Real}
    Nx, Ny     = size(f)
    plan       = plan_rfft(f)
    fk         = 1/(Nx*Ny).*(plan * f)

    Exp_to_cos_sin(fk)
end

function Fourier_cos_sin(dir::String, Quantity::Symbol, path_to_eos::String)
    if String(Quantity)[1] == 'T'
        f = TTTimeSeries(dir, Quantity)
    elseif String(Quantity)[1] == 'h'
        f = GWTimeSeries(dir, Quantity)
    elseif String(Quantity)[1] == 'v'
        f = LocalVEVsTimeSeries(dir, Quantity)
    elseif Quantity == :a4
        f = BoundaryTimeSeries(dir, Quantity)
    elseif String(Quantity)[1] == 'I'
        if Quantity == :Ideal_Txx
            f = IdealHydroTimeSeries(dir, path_to_eos, :px)
        elseif Quantity == :Ideal_Txy
            f = IdealHydroTimeSeries(dir, path_to_eos, :pxy)
        elseif Quantity == :Ideal_Tyy
            f = IdealHydroTimeSeries(dir, path_to_eos, :py)
        elseif Quantity == :Ideal_Tzz
            f = IdealHydroTimeSeries(dir, path_to_eos, :pz)
        end
    else
        f = VEVTimeSeries(dir, Quantity)
    end
    if String(Quantity)[1] == 'I' || String(Quantity)[1] == 'v'
        t, x, y = get_coords(BoundaryTimeSeries(dir, :a4),:, :, :)
    else
        t, x, y = get_coords(f, :, :, :)
    end
    Nt, Nx, Ny = size(f)
    dx         = x[2]-x[1]
    dy         = y[2]-y[1]
    kx         = 2*π.*(rfftfreq(Nx,1/dx)[1:end-1])
    ky         = 2*π.*(rfftfreq(Ny,1/dy)[1:end-1])
    Nkx        = length(kx)
    Nky        = length(ky)

    a = zeros(Nt, Nkx, Nky)
    b = zeros(Nt, Nkx, Nky)
    c = zeros(Nt, Nkx, Nky)
    d = zeros(Nt, Nkx, Nky)

    for n in 1:Nt
        a[n,:,:], b[n,:,:], c[n,:,:], d[n,:,:] = Fourier_cos_sin(f[n,:,:])
    end

    a, b, c, d, kx, ky
end

function Fourier_cos_sin_vivj_TT(dir::String)
    ut         = LocalVEVsTimeSeries(dir, :ut)
    ux         = LocalVEVsTimeSeries(dir, :ux)
    uy         = LocalVEVsTimeSeries(dir, :uy)
    t, x, y    = get_coords(BoundaryTimeSeries(dir,:a4),:,:,:)
    Nt, Nx, Ny = size(ut)
    dx         = x[2]-x[1]
    dy         = y[2]-y[1]
    kx         = 2*π.*(rfftfreq(Nx,1/dx)[1:end-1])
    ky         = 2*π.*(rfftfreq(Ny,1/dy)[1:end-1])
    Nkx        = length(kx)
    Nky        = length(ky)

    a = zeros(Nt, Nkx, Nky)
    b = zeros(Nt, Nkx, Nky)
    c = zeros(Nt, Nkx, Nky)
    d = zeros(Nt, Nkx, Nky)

    for n in 1:Nt
        vx2  = ux[n,:,:].^2 ./(ut[n,:,:])
        vy2  = uy[n,:,:].^2 ./(ut[n,:,:])
        vxvy = ux[n,:,:].*uy[n,:,:] ./(ut[n,:,:])

        avx2, bvx2, cvx2, dvx2     = Fourier_cos_sin(vx2)
        avy2, bvy2, cvy2, dvy2     = Fourier_cos_sin(vy2)
        avxvy, bvxvy, cvxvy, dvxvy = Fourier_cos_sin(vxvy)

        @inbounds for j in 1:Nky
            for i in 1:Nkx
                kkx      = kx[i]
                kky      = ky[j]
                kx2      = kkx^2
                ky2      = kky^2
                kxky     = kkx*kky
                a[n,i,j] = (ky2*avx2[i,j]+kx2*avy2[i,j]+2*kxky*avxvy[i,j])
                b[n,i,j] = (ky2*bvx2[i,j]+kx2*bvy2[i,j]+2*kxky*bvxvy[i,j])
                c[n,i,j] = (ky2*cvx2[i,j]+kx2*cvy2[i,j]+2*kxky*cvxvy[i,j])
                d[n,i,j] = (ky2*dvx2[i,j]+kx2*dvy2[i,j]+2*kxky*dvxvy[i,j])
            end
        end
    end

    a, b, c, d, kx, ky
end

function Fourier_cos_sin_dρd2k(dir::String)
    hdxx       = GWTimeSeries(dir, :hdxx)
    hdxy       = GWTimeSeries(dir, :hdxy)
    hdyy       = GWTimeSeries(dir, :hdyy)
    hdzz       = GWTimeSeries(dir, :hdzz)
    t, x, y    = get_coords(hdxx, :, :, :)
    Nt, Nx, Ny = size(hdxx)
    dx         = x[2]-x[1]
    dy         = y[2]-y[1]
    kx         = 2*π.*(rfftfreq(Nx,1/dx)[1:end-1])
    ky         = 2*π.*(rfftfreq(Ny,1/dy)[1:end-1])
    Nkx        = length(kx)
    Nky        = length(ky)
    plan       = 1/(Nx*Ny) * plan_rfft(hdxx[1,:,:])

    a = zeros(Nt, Nkx, Nky)
    b = zeros(Nt, Nkx, Nky)
    c = zeros(Nt, Nkx, Nky)
    d = zeros(Nt, Nkx, Nky)

    for n in 1:Nt
        hdkxx = abs.(plan * hdxx[n,:,:])
        hdkxy = abs.(plan * hdxy[n,:,:])
        hdkyy = abs.(plan * hdyy[n,:,:])
        hdkzz = abs.(plan * hdzz[n,:,:])
        f     = 2/π*(hdkxx.^2 + hdkxy.^2 + hdkyy.^2 + hdkzz.^2) .+ im*0

        a[n,:,:], b[n,:,:], c[n,:,:], d[n,:,:] = Exp_to_cos_sin(f)
    end

    a, b, c, d, kx, ky
end

function modes_to_mathematica(dirname::String, Quantity::Symbol, outfile::String; path_to_eos::String = "a")
    if Quantity == :drhod2k
        a, b, c, d, _, _ = Fourier_cos_sin_dρd2k(dirname)
    elseif Quantity == :vivj
        a, b, c, d, _, _ = Fourier_cos_sin_vivj_TT(dirname)
    else
        a, b, c, d, _, _ = Fourier_cos_sin(dirname, Quantity, path_to_eos)
    end
    Nt, Nkx, Nky     = size(a)
    if String(Quantity)[1] == 'h' || String(Quantity)[1] == 'd'
        t, _, _          = get_coords(GWTimeSeries(dirname, :hxx),:,1,1)
    else
        t, _, _          = get_coords(BoundaryTimeSeries(dirname, :a4),:,1,1)
    end

    a_m  = zeros(7, Nt, Nky, Nkx)

    n = 1
    @fastmath @inbounds for j in 1:Nkx
        for k in 1:Nky
            for i in 1:Nt
                a_m[:,i,k,j] = [j-1 k-1 t[i] abs(a[i,j,k]) abs(b[i,j,k]) abs(c[i,j,k]) abs(d[i,j,k])]
                n += 1
            end
        end
    end

    output   = abspath(dirname, outfile)
    fid      = h5open(output, "w")
    group_st = g_create(fid, "data")
    Jecco.write_dataset(group_st, "modes", a_m)
    close(fid)
end

function Fourier2D(dir::String, VEV::Symbol)
    f          = VEVTimeSeries(dir, VEV)
    t, x, y    = get_coords(f,:,:,:)
    Nt, Nx, Ny = size(f)
    dx         = x[2]-x[1]
    dy         = y[2]-y[1]
    kx         = 2*π.*(rfftfreq(Nx,1/dx))
    ky         = 2*π.*(fftfreq(Ny,1/dy))
    Nkx        = length(kx)
    Nky        = length(ky)

    plan = plan_rfft(f[1,:,:])
    fk   = im.*zeros(Nt, Nkx, Nky)

    for n in 1:Nt
        fk[n,:,:] = 1/(Nx*Ny).*(plan * f[n,:,:])
    end

    fk, kx, ky
end

function Fourier_ϕ(f::Array{T,2}) where {T<:Real}
    Nr, Nϕ = size(f)
    δϕ     = 2π/Nϕ
    a      = zeros(Nr, Int(floor(Nϕ/2))+1)
    b      = similar(a)
    @inbounds for n in 1:Nr
        plan   = plan_rfft(f[n,:])
        fk     = 1/Nϕ .* (plan * f[n,:])
        a[n,:] = 0.5 .*real.(fk)
        b[n,:] = -0.5 .*imag.(fk)
    end
    a, b
end

function Fourier_ϕ(dir::String, VEV::Symbol; time::Real = 0.0)
    f         = VEVTimeSeries(dir, VEV)
    t, x, y   = get_coords(f, :, :, :)
    δx, δy    = (x[2]-x[1], y[2]-y[1])
    Lx, Ly    = (x[end]+δx-x[1], y[end]+δy-y[1])
    n         = findfirst(t .>= time)
    _, Nx, Ny = size(f)
    rmax      = minimum((Lx/2-δx, Ly/2-δy))
    N         = minimum((Nx, Ny))
    r         = [(i-1)/N*rmax for i in 1:N]
    ϕ         = [2π*(i-1)/N for i in 1:N]
    f_inter   = interpolate((x, y, ), f[n,:,:], Gridded(Linear()))
    f_rϕ      = zeros(N, N)

    @inbounds @threads for j in 1:N
        for i in 1:N
            f_rϕ[i,j] = f_inter(r[i]*cos(ϕ[j]), r[i]*sin(ϕ[j]))
        end
    end

    r, Fourier_ϕ(f_rϕ)
end

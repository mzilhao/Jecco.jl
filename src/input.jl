
abstract type AbstractTimeSeries{N,T} end

struct FieldTimeSeries{N,T} <: AbstractTimeSeries{N,T}
    ts     :: T
    field  :: String
end

# slightly adapted from openPMD-viewer
mutable struct OpenPMDTimeSeries{S}
    iterations        :: Array{Int64}
    files             :: Array{String}
    current_i         :: Int64
    current_iteration :: Int64
    current_t         :: Float64
    params            :: S
end

function OpenPMDTimeSeries(foldername::String, prefix::String)
    iterations, files = try
        list_h5_files(foldername, prefix=prefix)
    catch e
        if isa(e, SystemError) && e.errnum == 2 # "No such file or directory"
            throw(ErrorException("No files found."))
        else
            throw(e)
        end
    end

    if length(iterations) == 0
        throw(ErrorException("No files found."))
    end

    # open first file to extract parameters
    fid    = h5open(files[1], "r")
    it     = iterations[1]
    grp, t = read_openpmd_file(fid, it)

    # read group attributes
    params = read_group_attributes(grp)

    OpenPMDTimeSeries{typeof(params)}(iterations, files, 1, it, t, params)
end

"""
    OpenPMDTimeSeries(foldername::String; prefix::String="")

Initialize an openPMD time series, ie, scan the directory and extract the
openPMD files starting with the given prefix (can be left empty)

# Example
```
julia> ts = OpenPMDTimeSeries("./data"; prefix="wave_")
```
"""
OpenPMDTimeSeries(foldername::String; prefix::String="") =
    OpenPMDTimeSeries(foldername, prefix)


FieldTimeSeries{N}(ts::OpenPMDTimeSeries, field::String) where{N} =
    FieldTimeSeries{N,typeof(ts)}(ts, field)

"""
    FieldTimeSeries(foldername::String; prefix::String, field::String)

Initialize a (openPMD) time series for the given `N`-dimensional `field`, ie,
scan the directory and extract the openPMD files starting with the given
`prefix`. The data corresponding to the given `field` will then be extracted
from the corresponding hdf5 file when requested.

# Example
```
julia> xi_ts = FieldTimeSeries("./", prefix="gauge_", field="xi")

julia> xi_ts[1:2,1,2:4,4]
2×3 Array{Float64,2}:
 0.0466351  0.0466351  0.0466351
 0.0738048  0.0738148  0.0738312
```
"""
function FieldTimeSeries(foldername::String; prefix::String, field::String)
    ts = OpenPMDTimeSeries(foldername, prefix)
    f, chart = get_field(ts, it=1, field=field)
    N = ndims(chart)
    FieldTimeSeries{N}(ts, field)
end

function Base.size(ff::FieldTimeSeries)
    Nt  = length(ff.ts.iterations)
    it0 = ff.ts.iterations[1]
    f0, = get_field(ff.ts, it=it0, field=ff.field)
    size_ = size(f0)
    Nt, size_...
end

@inline Base.firstindex(ff::FieldTimeSeries, d::Int) = 1
@inline Base.lastindex(ff::FieldTimeSeries, d::Int) = size(ff)[d]

function Base.getindex(ff::FieldTimeSeries, a::Int, idx::Vararg)
    it = ff.ts.iterations[a]
    f, = get_field(ff.ts, it=it, field=ff.field)
    f[idx...]
end

function Base.getindex(ff::FieldTimeSeries, aa, idx::Vararg)
    it  = ff.ts.iterations[aa[1]]
    f0, = get_field(ff.ts, it=it, field=ff.field)

    Na    = length(aa)
    size_ = size(f0[idx...])
    f     = zeros(Na, size_...)

    slicer = [Colon() for _ in 1:ndims(f0)]
    for (i,a) in enumerate(aa)
        it  = ff.ts.iterations[a]
        f0, = get_field(ff.ts, it=it, field=ff.field)
        f[i,slicer...] .= f0[idx...]
    end
    f
end

function Base.getindex(ff::FieldTimeSeries, ::Colon, idx::Vararg)
    Na  = length(ff.ts.iterations)
    getindex(ff, 1:Na, idx...)
end


"""
    get_coords(ff::FieldTimeSeries, a::Int, idx::Vararg)

# Example
```
julia> t, u, x, y = Jecco.get_coords(xi_ts,1:2,1,2:4,4)
([0.0, 0.10214880285375176], NaN, -4.6875:0.3125:-4.0625, -4.0625)
```
"""
function get_coords(ff::FieldTimeSeries, a::Int, idx::Vararg)
    it = ff.ts.iterations[a]
    f, chart = get_field(ff.ts, it=it, field=ff.field)
    t = ff.ts.current_t
    t, chart[idx...]...
end

"""
    get_coords(ff::FieldTimeSeries, aa, idx::Vararg)

# Example
```
julia> t, u, x, y = Jecco.get_coords(xi_ts,1,1,2,4)
(0.0, NaN, -4.6875, -4.0625)
```
"""
function get_coords(ff::FieldTimeSeries, aa, idx::Vararg)
    it = ff.ts.iterations[aa[1]]
    f0, chart = get_field(ff.ts, it=it, field=ff.field)
    t0 = ff.ts.current_t

    Na = length(aa)
    t  = Vector{typeof(t0)}(undef, Na)
    for (i,a) in enumerate(aa)
        it   = ff.ts.iterations[a]
        get_field(ff.ts, it=it, field=ff.field)
        t0   = ff.ts.current_t
        t[i] = t0
    end
    t, chart[idx...]...
end

function Jecco.get_coords(ff::FieldTimeSeries, ::Colon, idx::Vararg)
    Na = length(ff.ts.iterations)
    get_coords(ff, 1:Na, idx...)
end


function list_h5_files(foldername::String; prefix::String="")
    path     = abspath(foldername)
    allfiles = readdir(path)

    Ns = length(prefix)

    its_names = Tuple[]
    # append only the files whose names start with the given prefix
    for file in allfiles
        try
            if (file[1:Ns] == prefix && (file[end-2:end] == ".h5" ||
                                         file[end-4:end] == ".hdf5"))
                fullname = joinpath(path, file)
                # extract all iterations from the file
                fid      = h5open(fullname, "r")
                its      = keys(fid["/data"])
                close(fid)
                # for each iteration add to list of tuples with iteration and
                # name
                for it in its
                    push!(its_names, (parse(Int64, it), fullname))
                end
            end
        catch ex
            if isa(ex, BoundsError)
                # probably triggered by string comparison; do nothing
            else
                println("Error reading file $file")
                throw(ex)
            end
        end
    end

    # sort according to iteration
    sort!(its_names)
    # and extract the list of filenames and iterations
    filenames = [name for (it, name) in its_names]
    its       = [it for (it, name) in its_names]

    (its, filenames)
end

"""
    get_field(ts::OpenPMDTimeSeries; it::Int, field::String)

Given a time series, extract the requested field (and corresponding chart) from
an HDF5 file in the openPMD format. As side-effects, ```ts.current_i```,
```ts.current_iteration``` and ```ts.current_t``` are correspondingly modified.

# Example
```
julia> psi, chart=get_field(ts, it=20, field="psi");

julia> ts.current_t
20.0

julia> ts.current_iteration
20

julia> ts.current_i
21

```
"""
function get_field(ts::OpenPMDTimeSeries; it::Int, field::String, verbose::Bool=false)
    # index that corresponds to the closest iteration requested
    ts.current_i = argmin(abs.(it .- ts.iterations))
    # the closest iteration found (it need not be the requested one)
    ts.current_iteration = ts.iterations[ts.current_i]
    # and corresponding file
    filename = ts.files[ts.current_i]

    if verbose
        println("Reading file $filename")
    end
    # open file
    fid = h5open(filename, "r")

    # read in openPMD structure
    grp, ts.current_t = read_openpmd_file(fid, ts.current_iteration)

    # read actual data
    data, chart = read_dataset(grp, field)

    # close file
    close(fid)

    data, chart
end

function read_openpmd_file(fid::HDF5.File, it::Integer)
    basePath   = read_attribute(fid, "basePath")
    basePath   = replace(basePath, "%T" => it)
    meshesPath = read_attribute(fid, "meshesPath")

    # pointer to base group (ie, with the information for the requested time
    # level) within the given hdf5 file
    grp_base = fid[basePath]

    time = read_attribute(grp_base, "time")

    # pointer to mesh group (with the actual chart function data)
    grp_mesh = grp_base[meshesPath]

    grp_mesh, time
end

function read_group_attributes(grp::HDF5.Group)
    grp_attrs = attributes(grp)
    ks        = keys(grp_attrs)
    vals      = read_attribute.(Ref(grp), ks)
    Dict(ks[i] => vals[i] for i in 1:length(ks))
end

function read_dataset(grp::HDF5.Group, var::String)
    dset  = grp[var]

    func  = read(dset)
    nodes = size(func)
    T     = eltype(func)
    dim_  = length(nodes)

    if dim_ == 1
        names         = read_attribute(dset, "axisLabels")
        mins          = read_attribute(dset, "gridGlobalOffset")
        maxs          = read_attribute(dset, "gridMax")
        gridtypes     = read_attribute(dset, "gridType")
    else
        # remember to flip the order since HDF5 uses row-major order to store
        # arrays, as opposed to Julia's column-major order
        names         = read_attribute(dset, "axisLabels")[end:-1:1]
        mins          = read_attribute(dset, "gridGlobalOffset")[end:-1:1]
        maxs          = read_attribute(dset, "gridMax")[end:-1:1]
        gridtypes     = read_attribute(dset, "gridType")[end:-1:1]
    end

    chart = Chart(gridtypes, names, mins, maxs, nodes)

    func, chart
end

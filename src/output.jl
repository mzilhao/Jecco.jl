
using Dates
using Printf

mutable struct TimeInfo{T<:Real}
    it      :: Int
    t       :: T
    dt      :: T
    runtime :: T
end
TimeInfo(it::Int, t::T, dt::T, runtime::T) where {T} = TimeInfo{T}(it, t, dt, runtime)
TimeInfo() = TimeInfo(0, 0.0, 0.0, 0.0)

mutable struct Field{A,G}
    name  :: String
    data  :: A
    chart :: G
end
Field(name::String, data, chart::Chart) = Field{typeof(data),typeof(chart)}(name, data, chart)

function out_info(it::Integer, t::Real, time_per_hour::Real, f, label::String, info_every::Integer,
                  header_every::Integer)

    if it % header_every == 0
        println("-------------------------------------------------------------")
        println("Iteration      Time | Time per hour  |           $label")
        println("                    |                |   minimum      maximum")
        println("-------------------------------------------------------------")
    end

    if it % info_every == 0
        @printf "%9d %9.3f | %9.3f      | %9.4g    %9.4g\n" it t time_per_hour minimum(f) maximum(f)
    end

    nothing
end


struct Output{T}
    dir              :: String
    prefix           :: String
    software         :: String
    software_version :: String
    tinfo            :: TimeInfo{T}

    function Output{T}(dir::String, prefix::String,
                       software::String, software_version::String,
                       tinfo::TimeInfo{T}) where {T}
        # if no name specified error out
        if dir == ""
            error("Output folder cannot be empty string.")
        end

        # create folder if it doesn't exist already
        if !isdir(dir)
            mkdir(dir)
        end

        if prefix == ""
            prefix = "data_"
        end

        new(dir, prefix, software, software_version, tinfo)
    end
end
function Output(dir::String, prefix::String, tinfo::TimeInfo{T}) where {T<:Real}
    software         = "Jecco"
    software_version = "0.8.0"
    Output{T}(dir, prefix, software, software_version, tinfo)
end

# make Output a callable struct

function (out::Output)(fields::Union{Vector, Tuple}; params::Union{Tuple, NamedTuple, Dict}=())
    it   = out.tinfo.it
    time = out.tinfo.t
    dt   = out.tinfo.dt

    filename = "$(out.prefix)$(lpad(string(it), 8, string(0))).h5"
    fullpath = abspath(out.dir, filename)

    # check if file already exists
    if isfile(fullpath)
        mode = "r+"
        firsttime = false
    else
        mode = "cw"
        firsttime = true
    end

    # open file
    fid = h5open(fullpath, mode)

    # create openPMD structure
    if firsttime
        setup_openpmd_file(out, fid)
    end

    if haskey(fid, "data/$it")
        grp = fid["data/$it/fields"]
    else
        grp = create_group(fid, "data/$it")
        attributes(grp)["time"] = time
        attributes(grp)["dt"] = dt
        grp = create_group(grp, "fields")
    end

    # write given parameters as attributes
    if firsttime
        for key in keys(params)
            name = String(key)
            val  = params[key]
            attributes(grp)[name] = val
        end
    end

    for field in fields
        write_hdf5(out, grp, field)
    end

    # close file
    close(fid)

    nothing
end

(out::Output)(fields::Vararg{Field,N}; params=()) where {N} =
    (out::Output)((fields...,), params=params)


write_hdf5(out::Output, grp::HDF5.Group, field::Field) =
    write_hdf5(out, grp, field.name, field.data, field.chart)

function write_hdf5(out::Output, grp::HDF5.Group, fieldname::String, data::AbstractArray,
                    chart::Chart)
    dset, dtype = create_dataset(grp, fieldname, data)
    # write actual data
    write_dataset(dset, dtype, data)
    setup_openpmd_mesh(dset, chart)

    nothing
end

function setup_openpmd_file(out::Output, fid::HDF5.File)
    timenow = now()
    date    = Dates.format(timenow, "yyyy-mm-dd HH:MM:SS")

    attributes(fid)["software"] = out.software
    attributes(fid)["softwareVersion"] = out.software_version
    attributes(fid)["openPMD"] = "1.1.0"
    attributes(fid)["openPMDextension"] = 0
    attributes(fid)["author"] = try
        ENV["USER"]
    catch e
        if isa(e, KeyError) # probably on windows
            splitdir(homedir())[end]
        else
            throw(e)
        end
    end
    attributes(fid)["date"] = string(date)
    attributes(fid)["iterationEncoding"] = "fileBased"
    attributes(fid)["iterationFormat"] = "$(out.prefix)%T.h5"

    attributes(fid)["basePath"] = "/data/%T/"
    attributes(fid)["meshesPath"] = "fields/"
end

openpmd_geometry(coord::CartesianCoord) = "cartesian"
openpmd_geometry(coord::GaussLobattoCoord) = "other"

function openpmd_geometry(chart::Chart)
    geometries = [openpmd_geometry(coord) for coord in chart.coords]
    if all(geometries .== "cartesian")
        geometry = "cartesian"
    else
        geometry = "other"
    end
    geometry
end

function setup_openpmd_mesh(dset::HDF5.Dataset, coord::AbstractCoord)
    attributes(dset)["geometry"]         = openpmd_geometry(coord)
    attributes(dset)["gridGlobalOffset"] = coord.min
    attributes(dset)["gridSpacing"]      = Jecco.delta(coord)
    attributes(dset)["gridMax"]          = coord.max
    attributes(dset)["gridType"]         = Jecco.coord_type(coord)
    attributes(dset)["axisLabels"]       = coord.name
    nothing
end
# there is no method to write BigFloat attributes in HDF5, so convert them to Float64
function setup_openpmd_mesh(dset::HDF5.Dataset, coord::AbstractCoord{N,BigFloat}) where {N}
    T = Float64
    attributes(dset)["geometry"]         = openpmd_geometry(coord)
    attributes(dset)["gridGlobalOffset"] = T(coord.min)
    attributes(dset)["gridSpacing"]      = T(Jecco.delta(coord))
    attributes(dset)["gridMax"]          = T(coord.max)
    attributes(dset)["gridType"]         = T(Jecco.coord_type(coord))
    attributes(dset)["axisLabels"]       = coord.name
    nothing
end

function setup_openpmd_mesh(dset::HDF5.Dataset, chart::Chart)
    mins      = Jecco.min(chart)
    deltas    = Jecco.delta(chart)
    maxs      = Jecco.max(chart)
    gridtypes = Jecco.coord_type(chart)
    names     = Jecco.name(chart)

    attributes(dset)["geometry"]         = openpmd_geometry(chart)

    # Julia, like Fortran and Matlab, stores arrays in column-major order. HDF5
    # uses C's row-major order, and consequently every array's dimensions are
    # inverted compared to what is seen with tools like h5dump. This is the same
    # convention as for the Fortran and Matlab HDF5 interfaces. For consistency,
    # then, we here flip the order of the grid arrays (and remember to flip back
    # when reading data in). The advantage is that no data rearrangement takes
    # place when reading or writing the data itself, which is more intensive.
    attributes(dset)["gridGlobalOffset"] = mins[end:-1:1]
    attributes(dset)["gridSpacing"]      = deltas[end:-1:1]
    attributes(dset)["gridMax"]          = maxs[end:-1:1]
    attributes(dset)["gridType"]         = gridtypes[end:-1:1]
    attributes(dset)["axisLabels"]       = names[end:-1:1]
    nothing
end
function setup_openpmd_mesh(dset::HDF5.Dataset, chart::Chart{N,BigFloat}) where{N}
    T = Float64
    mins      = T.(Jecco.min(chart))
    deltas    = T.(Jecco.delta(chart))
    maxs      = T.(Jecco.max(chart))
    gridtypes = Jecco.coord_type(chart)
    names     = Jecco.name(chart)

    attributes(dset)["geometry"]         = openpmd_geometry(chart)
    attributes(dset)["gridGlobalOffset"] = mins[end:-1:1]
    attributes(dset)["gridSpacing"]      = deltas[end:-1:1]
    attributes(dset)["gridMax"]          = maxs[end:-1:1]
    attributes(dset)["gridType"]         = gridtypes[end:-1:1]
    attributes(dset)["axisLabels"]       = names[end:-1:1]
    nothing
end

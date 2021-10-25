using RecipesBase

@recipe function plot(f_ts::TimeSeries, ii::Vararg)
    tmp = get_coords(f_ts,ii...)
    f   = f_ts[ii...]

    # we want to select only those that are *not* Numbers (ie, Arrays or Ranges)
    idx = .!isa.(tmp, Number)
    idx = [idx...]
    coord = [tmp...]

    xx = coord[idx]

    xx_f = [xx..., f]

    match_dimensions := true

    (xx_f...,)
end

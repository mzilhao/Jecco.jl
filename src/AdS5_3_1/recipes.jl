using RecipesBase

@recipe function plot(f_ts::TimeSeries{2}, i0::Union{Colon, UnitRange}, i1::Int, i2::Int)
    xguide := "t"
    tt, xx, yy = get_coords(f_ts,i0,i1,i2)
    f = f_ts[i0,i1,i2]
    tt, f
end

@recipe function plot(f_ts::TimeSeries{2}, i0::Int, i1::Union{Colon, UnitRange}, i2::Int)
    xguide := "x"
    tt, xx, yy = get_coords(f_ts,i0,i1,i2)
    f = f_ts[i0,i1,i2]
    xx, f
end

@recipe function plot(f_ts::TimeSeries{2}, i0::Int, i1::Int, i2::Union{Colon, UnitRange})
    xguide := "y"
    tt, xx, yy = get_coords(f_ts,i0,i1,i2)
    f = f_ts[i0,i1,i2]
    yy, f
end

@recipe function plot(f_ts::TimeSeries{2}, i0::Int, i1::Union{Colon, UnitRange}, i2::Union{Colon, UnitRange})
    xguide := "y"
    yguide := "x"
    tt, xx, yy = get_coords(f_ts,i0,i1,i2)
    f = f_ts[i0,i1,i2]
    yy, xx, f
end

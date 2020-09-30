using RecipesBase
using Printf

@recipe function plot(f_ts::TimeSeries{3}, i0, i1, i2, i3)
    f = f_ts[i0,i1,i2,i3]
    coord = get_coords(f_ts, i0,i1,i2,i3)

    ix_arr = isa.([i0,i1,i2,i3], Union{Colon,UnitRange})

    if sum(ix_arr) != 1
        error("Need to have a 1-D array...")
    end

    # direction along which we're plotting
    ix = argmin(abs.(ix_arr .- 1))
    xx = coord[ix]

    if ix == 1
        xguide := "t"

        xc = "x =  $(@sprintf("%.2f", coord[3]))"
        yc = ", y =  $(@sprintf("%.2f", coord[4]))"

        if !isnan(coord[2])
            uc = ", u =  $(@sprintf("%.2f", coord[2]))"
        else
            uc = ""
        end
        label --> xc * yc * uc

    elseif ix == 2
        xguide := "u"

        tc = "t =  $(@sprintf("%.2f", coord[1]))"
        xc = ", x =  $(@sprintf("%.2f", coord[3]))"
        yc = ", y =  $(@sprintf("%.2f", coord[4]))"

        label --> tc * xc * yc

    elseif ix == 3
        xguide := "x"

        tc = "t =  $(@sprintf("%.2f", coord[1]))"
        yc = ", y =  $(@sprintf("%.2f", coord[4]))"

        if !isnan(coord[2])
            uc = ", u =  $(@sprintf("%.2f", coord[2]))"
        else
            uc = ""
        end
        label --> tc * yc * uc

    elseif ix == 4
        xguide := "y"

        tc = "t =  $(@sprintf("%.2f", coord[1]))"
        xc = ", x =  $(@sprintf("%.2f", coord[3]))"

        if !isnan(coord[2])
            uc = ", u =  $(@sprintf("%.2f", coord[2]))"
        else
            uc = ""
        end
        label --> tc * yc * uc
    end

    xx, f
end

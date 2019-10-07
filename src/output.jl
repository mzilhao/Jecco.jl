
using Printf

function out_info(it::Integer, t::Real, f, label::String, info_every::Integer,
                  header_every::Integer)

    if it % header_every == 0
        println("--------------------------------------------------------------------------------")
        println(" it \t   Time     \t  min($label) \t max($label)")
        println("--------------------------------------------------------------------------------")
    end

    if it % info_every == 0
        @printf " %d \t %6.3f \t %g \t %g\n" it t minimum(f) maximum(f)
    end

    nothing
end

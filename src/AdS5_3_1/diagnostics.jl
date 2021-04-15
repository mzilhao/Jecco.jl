
struct NoDiag <: Diagnostics end

"""
Look for the Apparent Horizon during the evolution
"""
Base.@kwdef struct DiagAH <: Diagnostics
    # negative values turns off
    find_AH_every       :: Int  = -1
    find_AH_every_t     :: Float32  = -1.0
end


function (::NoDiag)(bulkevols, bulkconstrains, bulkderivs, boundary::Boundary,
                    gauge::Gauge, horizoncache::HorizonCache, systems::SystemPartition,
                    tinfo::Jecco.TimeInfo, evoleq::EvolutionEquations, io::InOut)
    # no-op
    () -> nothing
end


# initial guess for AH surface
function sigma_guess!(sigma, gaugecondition::ConstantAH)
    uAH = gaugecondition.u_AH
    fill!(sigma, 1/uAH)
end

function sigma_guess!(sigma, gaugecondition::GaugeCondition)
    # for other gauge conditions, just initialize sigma to 1, for lack of better
    # alternative
    fill!(sigma, 1)
end

function (diagnostics::DiagAH)(
    bulkevols, bulkconstrains, bulkderivs, boundary::Boundary, gauge::Gauge,
    horizoncache::HorizonCache, systems::SystemPartition, tinfo::Jecco.TimeInfo,
    evoleq::EvolutionEquations, io::InOut)

    # Apparent Horizon surface
    sigma = similar(gauge.xi)

    empty   = Cartesian{1}("u", systems[1].ucoord[1], systems[1].ucoord[1], 1)
    chart2D = Chart(empty, systems[end].xcoord, systems[end].ycoord)

    out_diag  = Jecco.Output(io.out_dir, "diagnostics_", tinfo)

    fields = Jecco.Field("sigma", sigma, chart2D)

    last_output_AH = -diagnostics.find_AH_every_t

    function ()
        it = tinfo.it
        tt = tinfo.t

        do_find_AH = false

        if (diagnostics.find_AH_every > 0 && it % diagnostics.find_AH_every == 0)
            do_find_AH = true
        end

        if (diagnostics.find_AH_every_t > 0 &&
            tt >= last_output_AH + diagnostics.find_AH_every_t - 1e-12)
            do_find_AH     = true
            last_output_AH = tt
        end

        if do_find_AH
            # it's useful to reset the initial guess, as the finder may sometimes fail
            # to find the AH, and the sigma surface would then be very displaced.
            sigma_guess!(sigma, evoleq.gaugecondition)
            find_AH!(sigma, bulkconstrains[end], bulkevols[end], bulkderivs[end], gauge,
                     horizoncache, systems[end], evoleq.ahf)

            # write sigma output
            out_diag(fields)
        end

        nothing
    end
end

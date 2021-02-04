

function (filters::Filters)(bulkevol::BulkEvolved)
    B1  = getB1(bulkevol)
    B2  = getB2(bulkevol)
    G   = getG(bulkevol)
    phi = getphi(bulkevol)

    @sync begin
        @spawn filters.exp_filter(B1)
        @spawn filters.exp_filter(B2)
        @spawn filters.exp_filter(G)
        @spawn filters.exp_filter(phi)
    end
    nothing
end

function setup_rhs(bulkconstrains::BulkPartition{Nsys}, bulkderivs::BulkPartition{Nsys},
                   cache::HorizonCache, systems::SystemPartition,
                   integration::Integration) where {Nsys}
    # function to solve the nested system
    nested = Nested(systems, bulkconstrains, bulkderivs)

    function (ff_t::EvolVars, ff::EvolVars, evoleq::EvolutionEquations, t)
        bulkevols_t = getbulkevolvedpartition(ff_t)
        boundary_t  = getboundary(ff_t)
        gauge_t     = getgauge(ff_t)

        bulkevols   = getbulkevolvedpartition(ff)
        boundary    = getboundary(ff)
        gauge       = getgauge(ff)

        compute_boundary_t!(boundary_t, bulkevols[1], boundary, gauge, systems[1], evoleq)
        apply_dissipation!(boundary_t, boundary,  systems[1])

        # solve nested system for the constrained variables
        nested(bulkevols, boundary, gauge, evoleq)

        compute_xi_t!(gauge_t, bulkconstrains[Nsys], bulkevols[Nsys], bulkderivs[Nsys],
                      gauge, cache, systems[Nsys], evoleq.gaugecondition)
        apply_dissipation!(gauge_t, gauge,  systems[Nsys])

        @inbounds @threads for aa in 1:Nsys
            sys           = systems[aa]
            bulkevol_t    = bulkevols_t[aa]
            bulkevol      = bulkevols[aa]
            bulkconstrain = bulkconstrains[aa]

            compute_bulkevolved_t!(bulkevol_t, bulkconstrain, gauge_t, bulkevol,
                                   boundary, gauge, sys, evoleq)
            apply_dissipation!(bulkevol_t, bulkevol, sys)
            # exponential filter
            sys.filters(bulkevol_t)
        end
        sync_bulkevolved!(bulkevols_t, bulkconstrains, gauge_t, systems, evoleq)

        nothing
    end
end

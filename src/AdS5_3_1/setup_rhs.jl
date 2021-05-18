
function (filters::Filters)(boundary::Boundary)
    filters.ko_filter2D_x(boundary.a4)
    filters.ko_filter2D_x(boundary.fx2)
    filters.ko_filter2D_x(boundary.fy2)
    filters.ko_filter2D_y(boundary.a4)
    filters.ko_filter2D_y(boundary.fx2)
    filters.ko_filter2D_y(boundary.fy2)
    nothing
end

function (filters::Filters)(gauge::Gauge)
    filters.ko_filter2D_x(gauge.xi)
    filters.ko_filter2D_y(gauge.xi)
    nothing
end

function (filters::Filters)(bulkevol::BulkEvolved)
    @sync begin
        @spawn filters.ko_filter_x(bulkevol.B1)
        @spawn filters.ko_filter_x(bulkevol.B2)
        @spawn filters.ko_filter_x(bulkevol.G)
        @spawn filters.ko_filter_x(bulkevol.phi)
    end
    @sync begin
        @spawn filters.ko_filter_y(bulkevol.B1)
        @spawn filters.ko_filter_y(bulkevol.B2)
        @spawn filters.ko_filter_y(bulkevol.G)
        @spawn filters.ko_filter_y(bulkevol.phi)
    end
    @sync begin
        @spawn filters.exp_filter(bulkevol.B1)
        @spawn filters.exp_filter(bulkevol.B2)
        @spawn filters.exp_filter(bulkevol.G)
        @spawn filters.exp_filter(bulkevol.phi)
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

        # filter after each integration (sub)step
        if t > 0 && integration.filter_poststep
            @inbounds for aa in 1:Nsys
                sys = systems[aa]
                sys.filters(bulkevols[aa])
            end
            systems[1].filters(boundary)
            systems[Nsys].filters(gauge)
        end

        compute_boundary_t!(boundary_t, bulkevols[1], boundary, gauge, systems[1], evoleq)

        # solve nested system for the constrained variables
        nested(bulkevols, boundary, gauge, evoleq)

        compute_xi_t!(gauge_t, bulkconstrains[Nsys], bulkevols[Nsys], bulkderivs[Nsys],
                      gauge, cache, systems[Nsys], evoleq.gaugecondition)

        # TODO: check if this loop is thread-safe
        # @inbounds @threads for aa in 1:Nsys
        @inbounds for aa in 1:Nsys
            sys           = systems[aa]
            bulkevol_t    = bulkevols_t[aa]
            bulkevol      = bulkevols[aa]
            bulkconstrain = bulkconstrains[aa]

            compute_bulkevolved_t!(bulkevol_t, bulkconstrain, gauge_t, bulkevol,
                                   boundary, gauge, sys, evoleq)
        end
        sync_bulkevolved!(bulkevols_t, bulkconstrains, gauge_t, systems, evoleq)

        nothing
    end
end

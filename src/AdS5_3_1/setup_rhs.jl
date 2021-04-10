
function (filters::Filters)(boundary::Boundary)
    a4  = geta4(boundary)
    fx2 = getfx2(boundary)
    fy2 = getfy2(boundary)

    filters.ko_filter2D_x(a4)
    filters.ko_filter2D_x(fx2)
    filters.ko_filter2D_x(fy2)
    filters.ko_filter2D_y(a4)
    filters.ko_filter2D_y(fx2)
    filters.ko_filter2D_y(fy2)
    nothing
end

function (filters::Filters)(gauge::Gauge)
    xi  = getxi(gauge)

    filters.ko_filter2D_x(xi)
    filters.ko_filter2D_y(xi)
    nothing
end


function (filters::Filters)(bulkevol::BulkEvolved)
    B1  = getB1(bulkevol)
    B2  = getB2(bulkevol)
    G   = getG(bulkevol)
    phi = getphi(bulkevol)

    @sync begin
        @spawn filters.ko_filter_x(B1)
        @spawn filters.ko_filter_x(B2)
        @spawn filters.ko_filter_x(G)
        @spawn filters.ko_filter_x(phi)
    end
    @sync begin
        @spawn filters.ko_filter_y(B1)
        @spawn filters.ko_filter_y(B2)
        @spawn filters.ko_filter_y(G)
        @spawn filters.ko_filter_y(phi)
    end
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

```
In what follows we removed the @threads for the parallelization due to
suspicion for possible problem with race conditions. The motivation
comes from a similar problem in the KG branch, that has similar
structure as below.

The parallelization was implemented via putting @threads after the
@inbounds in the two loops below.

TODO: fix the parallelization below
```

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


struct EvolTest0 <: AbstractEvolEq end

Base.@kwdef struct EvolEq{T,TP<:Potential} <: AbstractEvolEq
    phi0          :: T   = 0.0
    potential     :: TP  = ZeroPotential()
end

function setup_rhs(bulkconstrains::BulkPartition{Nsys}, systems::SystemPartition,
                   evoleq::AbstractEvolEq) where {Nsys}

    # function to solve the nested system
    solve_nested! = nested_solver(systems, evoleq)

    function (ff_t::EvolVars, ff::EvolVars, t)
        bulkevols_t = getbulkevolvedpartition(ff_t)
        boundary_t  = getboundary(ff_t)
        gauge_t     = getgauge(ff_t)

        bulkevols   = getbulkevolvedpartition(ff)
        boundary    = getboundary(ff)
        gauge       = getgauge(ff)

        # TODO
        # compute_boundary_t!(boundary_t, bulkevols[1], boundary, gauge, systems[1], evoleq)

        # solve nested system for the constrained variables
        solve_nested!(bulkconstrains, bulkevols, boundary, gauge)

        # TODO
        # compute_xi_t!(gauge_t, bulkconstrains[end], bulkevols[end], boundary, gauge, systems[end], evoleq)

        for aa in 1:Nsys
            sys = systems[aa]
            bulkevol_t = bulkevols_t[aa]

            # TODO
            # compute_bulkevolved_t!(bulkevol_t, bulkconstrains, bulkevols, boundary, gauge, systems, evoleq)
        end

        nothing
    end
end


struct EvolTest0 <: AbstractEvolEq end

Base.@kwdef struct EvolEq{T,TP<:Potential} <: AbstractEvolEq
    phi0          :: T   = 0.0
    potential     :: TP  = ZeroPotential()
end


function setup_rhs(bulkconstrains, systems::SystemPartition, evoleq::AbstractEvolEq)
    Nsys = length(systems)

    # function to solve the nested system
    solve_nested = nested_solver(systems, evoleq)

    function (ff_t::EvolPartition, ff::EvolPartition, t)
        bulkevols_t = getbulkevols(ff_t)
        boundary_t  = getboundary(ff_t)
        gauge_t     = getgauge(ff_t)

        bulkevols   = getbulkevols(ff)
        boundary    = getboundary(ff)
        gauge       = getgauge(ff)


        # TODO
        # compute_boundary_t!(boundary_t, bulkevols[1], boundary, gauge, systems[1], evoleq)

        # all bulk variables
        bulks = Bulk.(bulkevols, bulkconstrains)

        @show typeof(bulks)

        # solve nested system
        solve_nested(bulks, boundary, gauge)

        # TODO
        # compute_xi_t!(gauge_t, bulks[end], boundary, gauge, systems[end], evoleq)

        for aa in 1:Nsys
            sys = systems[aa]
            bulk       = bulks[aa]
            bulkevol_t = bulkevols_t[aa]

            # TODO
            # compute_bulk_t!(bulkevol_t, bulks, boundary, gauge, systems, evoleq)
        end

        nothing
    end
end



# TODO
function get_f_t!(evol_t::EvolPartition, bulks, boundary::Boundary, gauge::Gauge,
                  systems::SystemPartition, ::EvolTest0)
    Nsys = length(systems)

    bulkevols_t = getbulkevols(evol_t)
    boundary_t  = getboundary(evol_t)
    gauge_t     = getgauge(evol_t)


    fill!(boundary_t.a4, 0)
    fill!(boundary_t.fx2, 0)
    fill!(boundary_t.fy2, 0)
    fill!(gauge_t.xi, 0)

    for aa in 1:Nsys
        sys = systems[aa]
        bulk       = bulks[aa]
        bulkevol_t = bulkevols_t[aa]

        fill!(bulkevol_t.B1,  0)
        fill!(bulkevol_t.B2,  0)
        fill!(bulkevol_t.G,   0)
        fill!(bulkevol_t.phi, 0)
    end
    nothing
end


# TODO
function get_f_t!(evol_t::EvolPartition, bulks, boundary::Boundary, gauge::Gauge,
                  systems::SystemPartition, evoleq::EvolEq)
    Nsys = length(systems)

    bulkevols_t = getbulkevols(evol_t)
    boundary_t  = getboundary(evol_t)
    gauge_t     = getgauge(evol_t)


    fill!(boundary_t.a4, 0)
    fill!(boundary_t.fx2, 0)
    fill!(boundary_t.fy2, 0)
    fill!(gauge_t.xi, 0)

    for aa in 1:Nsys
        sys = systems[aa]
        bulk       = bulks[aa]
        bulkevol_t = bulkevols_t[aa]

        fill!(bulkevol_t.B1,  0)
        fill!(bulkevol_t.B2,  0)
        fill!(bulkevol_t.G,   0)
        fill!(bulkevol_t.phi, 0)
    end
    nothing
end

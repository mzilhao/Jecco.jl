
struct EvolTest0 <: AbstractEvolEq end

Base.@kwdef struct EvolEq{T,TP<:Potential} <: AbstractEvolEq
    phi0          :: T   = 0.0
    potential     :: TP  = ZeroPotential()
end


function setup_rhs(bulks, systems::SystemPartition, evoleq::AbstractEvolEq)
    Nsys = length(systems)

    # function to solve the nested system
    solve_nested = nested_solver(systems, evoleq)

    function (ff_t::EvolPartition, ff::EvolPartition, t)
        boundary    = getboundary(ff)
        gauge       = getgauge(ff)

        # solve nested system
        solve_nested(bulks, boundary, gauge)

        get_f_t!(ff_t, bulks, boundary, gauge, systems, evoleq)

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

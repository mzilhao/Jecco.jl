
struct EvolTest0 <: AbstractEvolEq end

function get_f_t!(ff_t, ff, systems, evoleq::EvolTest0)
    Nsys = length(systems)

    boundary    = getboundary(ff)
    gauge       = getgauge(ff)
    bulkevols   = getbulkevols(ff)

    boundary_t  = getboundary(ff_t)
    gauge_t     = getgauge(ff_t)
    bulkevols_t = getbulkevols(ff_t)

    fill!(boundary_t.a4, 0)
    fill!(boundary_t.fx2, 0)
    fill!(boundary_t.fy2, 0)
    fill!(gauge_t.xi, 0)

    for aa in 1:Nsys
        sys = systems[aa]
        bulkevol   = bulkevols[aa]
        bulkevol_t = bulkevols_t[aa]

        fill!(bulkevol_t.B1,  0)
        fill!(bulkevol_t.B2,  0)
        fill!(bulkevol_t.G,   0)
        fill!(bulkevol_t.phi, 0)
    end
    nothing
end

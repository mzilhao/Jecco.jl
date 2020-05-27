
function compute_bulkevolved_t!(bulkevol_t::BulkEvolved,
                                bulkconstrain::BulkConstrained, gauge_t::Gauge,
                                bulkevol::BulkEvolved, boundary::Boundary,
                                gauge::Gauge, sys::System, ::EvolTest0)

    B1_t, B2_t, G_t, phi_t = unpack(bulkevol_t)
    # B1  , B2  , G  , phi   = unpack(bulkevol)

    fill!(B1_t,  0)
    fill!(B2_t,  0)
    fill!(G_t,   0)
    fill!(phi_t, 0)

    nothing
end

# TODO
function compute_bulkevolved_t!(bulkevol_t::BulkEvolved,
                                bulkconstrain::BulkConstrained, gauge_t::Gauge,
                                bulkevol::BulkEvolved, boundary::Boundary,
                                gauge::Gauge, sys::System, evoleq::AffineNull)

    B1_t, B2_t, G_t, phi_t = unpack(bulkevol_t)
    # B1  , B2  , G  , phi   = unpack(bulkevol)

    fill!(B1_t,  0)
    fill!(B2_t,  0)
    fill!(G_t,   0)
    fill!(phi_t, 0)

    nothing
end


# TODO
function compute_xi_t!(gauge_t::Gauge, bulk::BulkConstrained, gauge::Gauge,
                       sys::System{Outer}, ::EvolutionEquations)
    xi_t = getxi(gauge_t)
    fill!(xi_t, 0)
    nothing
end

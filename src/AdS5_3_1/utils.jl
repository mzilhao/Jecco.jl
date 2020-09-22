
"""
    compute_energy(a4, phi2, phi0, phiM2)

Returns the normalized energy density ``ϵ / ϕ₀⁴``.
"""
function compute_energy(a4, phi2, phi0, oophiM2)
    if phi0 < 1.e-9
        return -3//4 * a4
    else
        return -3//4 * a4 / phi0^4 .- phi2 / phi0^3 .- 1//4 * oophiM2 .+ 7//36
    end
end

function get_energy(ts::OpenPMDTimeSeries; it::Int, phi0, oophiM2, verbose::Bool=false)
    a4,   chart = get_field(ts, it=it, field="a4", verbose=verbose)
    phi2, chart = get_field(ts, it=it, field="phi2", verbose=verbose)

    en = compute_energy(a4, phi2, phi0, oophiM2)

    en, chart
end

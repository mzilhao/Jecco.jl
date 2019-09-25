
# TODO: does inline work with the if statement??
@inline function dphifdt(u::Real, Af::Real, phif::Real, phif_du::Real, phifd::Real,
                         phifd_du::Real)
    if u > 1.e-8
        return 1.5 * (phif - phifd) / u + 0.5 * phif_du + 1.5 * u * u * u * Af * phif +
            0.5 * Af * u * u * u * u * phif_du
    else
        return 1.5 * (phif_du - phifd_du) + 0.5 * phif_du
    end
end


function rhs!(dphif::Array, bulk::BulkVars, sys)
    coords = sys.coords
    derivs = sys.derivs

    uu, xx, yy = Vivi.xx(coords)

    Du_phi  = Vivi.D(bulk.phi, derivs[1], 1)
    Du_phid = Vivi.D(bulk.phid, derivs[1], 1)

    for j in eachindex(yy)
        for i in eachindex(xx)
            for a in eachindex(uu)
                u        = uu[a]
                phif     = bulk.phi[a,i,j]
                phifd    = bulk.phid[a,i,j]
                Af       = bulk.A[a,i,j]
                phif_du  = Du_phi[a,i,j]
                phifd_du = Du_phid[a,i,j]

                dphif[a,i,j] = dphifdt(u, Af, phif, phif_du, phifd, phifd_du)
            end
        end
    end

    nothing
end

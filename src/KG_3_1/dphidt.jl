

## TODO: separate


# TODO: does inline work with the if statement??
@inline function dphifdt(vars::AllVars)
    u         =    vars.u
    phi       =    vars.phi_d0
    phid      =    vars.phid_d0
    A         =    vars.A_d0
    phi_du    =    vars.phi_du
    phid_du   =    vars.phid_du

    if u > 1.e-8
        return 1.5 * (phi - phid) / u + 0.5 * phi_du + 1.5 * u * u * u * A * phi +
            0.5 * A * u * u * u * u * phi_du
    else
        return 1.5 * (phi_du - phid_du) + 0.5 * phi_du
    end
end


function rhs!(dphif::Array, bulk::BulkVars, sys)
    coords = sys.coords
    derivs = sys.derivs

    uu, xx, yy = Vivi.xx(coords)

    Du_phi  = Vivi.D(bulk.phi, derivs[1], 1)
    Du_phid = Vivi.D(bulk.phid, derivs[1], 1)

    vars  = AllVars{Float64}()

    for j in eachindex(yy)
        for i in eachindex(xx)
            for a in eachindex(uu)
                vars.u       = uu[a]

                vars.phi_d0  = bulk.phi[a,i,j]
                vars.phid_d0 = bulk.phid[a,i,j]
                vars.A_d0    = bulk.A[a,i,j]

                vars.phi_du  = Du_phi[a,i,j]
                vars.phid_du = Du_phid[a,i,j]

                dphif[a,i,j] = dphifdt(vars)
            end
        end
    end

    nothing
end

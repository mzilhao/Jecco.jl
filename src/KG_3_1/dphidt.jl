
@inline function dphig1dt(vars::AllVars)
    u         =    vars.u
    phi       =    vars.phi_d0
    phid      =    vars.phid_d0
    A         =    vars.A_d0
    phi_du    =    vars.phi_du

    1.5 * (phi - phid) / u + 0.5 * phi_du + 1.5 * u * u * u * A * phi +
        0.5 * A * u * u * u * u * phi_du
end

@inline function dphig1dt_u0(vars::AllVars)
    phi_du    =    vars.phi_du
    phid_du   =    vars.phid_du

    1.5 * (phi_du - phid_du) + 0.5 * phi_du
end


function rhs_g1!(dphif::Array, bulk::BulkVars, sys)
    coords = sys.coords
    derivs = sys.derivs

    uu, xx, yy = Vivi.xx(coords)
    Nu = length(uu)

    Du_phi  = Vivi.D(bulk.phi, derivs[1], 1)
    Du_phid = Vivi.D(bulk.phid, derivs[1], 1)

    vars  = AllVars{Float64}()

    # u = 0
    # TODO: parallelize here
    @fastmath @inbounds for j in eachindex(yy)
        @inbounds for i in eachindex(xx)
            vars.u       = uu[1]

            vars.phi_d0  = bulk.phi[1,i,j]
            vars.phid_d0 = bulk.phid[1,i,j]
            vars.A_d0    = bulk.A[1,i,j]

            vars.phi_du  = Du_phi[1,i,j]
            vars.phid_du = Du_phid[1,i,j]

            dphif[a,i,j] = dphig1dt_u0(vars)
        end
    end

    # TODO: parallelize here
    @fastmath @inbounds for j in eachindex(yy)
        @inbounds for i in eachindex(xx)
            @inbounds @simd for a in 2:Nu
                vars.u       = uu[a]

                vars.phi_d0  = bulk.phi[a,i,j]
                vars.phid_d0 = bulk.phid[a,i,j]
                vars.A_d0    = bulk.A[a,i,j]

                vars.phi_du  = Du_phi[a,i,j]
                vars.phid_du = Du_phid[a,i,j]

                dphif[a,i,j] = dphig1dt(vars)
            end
        end
    end

    nothing
end

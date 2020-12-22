
# only for the first u-domain, since we need to separate the u=0 point
function compute_bulkevolved_t_1st!(bulkevol_t::BulkEvolved,
                                    bulkconstrain::BulkConstrained,
                                    bulkevol::BulkEvolved, deriv::BulkDeriv,
                                    sys::System, evoleq::AffineNull)
    phiGF   = getphi(bulkevol)
    phiGF_t = getphi(bulkevol_t)

    phidGF = getphid(bulkconstrain)
    # SdGF   = getSd(bulkconstrain)
    AGF    = getA(bulkconstrain)

    Du_phi  = deriv.Du_phi
    Du_phid = deriv.Du_phid

    uu  = sys.ucoord
    # Du  = sys.Du
    # Dx  = sys.Dx
    # Dy  = sys.Dy

    Nu, Nx, Ny = size(sys)

    # u = 0
    @fastmath @inbounds for j in 1:Ny
        @inbounds @simd for i in 1:Nx
            phi_u  = Du_phi[1,i,j]
            phid_u = Du_phid[1,i,j]

            phiGF_t[1,i,j] = 3//2 * (phi_u - phid_u) + phi_u / 2
        end
    end

    # remaining inner grid points
    @fastmath @inbounds for j in 1:Ny
        @inbounds for i in 1:Nx
            @inbounds @simd for a in 2:Nu
                u      = uu[a]
                u3     = u * u * u
                u4     = u * u3

                phi    = phiGF[a,i,j]
                phid   = phidGF[a,i,j]
                A      = AGF[a,i,j]

                phi_u  = Du_phi[a,i,j]
                # phid_u = Du_phid[a,i,j]

                phiGF_t[a,i,j] =  3//2 * (phi - phid) / u + phi_u / 2 + 3//2 * u3 * A * phi +
                    A * u4 * phi_u / 2
            end
        end
    end

    nothing
end

# remaining u-domains
function compute_bulkevolved_t!(bulkevol_t::BulkEvolved,
                                bulkconstrain::BulkConstrained,
                                bulkevol::BulkEvolved, deriv::BulkDeriv,
                                sys::System, evoleq::AffineNull)
    phiGF   = getphi(bulkevol)
    phiGF_t = getphi(bulkevol_t)

    phidGF = getphid(bulkconstrain)
    # SdGF   = getSd(bulkconstrain)
    AGF    = getA(bulkconstrain)

    Du_phi  = deriv.Du_phi
    Du_phid = deriv.Du_phid

    uu  = sys.ucoord
    # Du  = sys.Du
    # Dx  = sys.Dx
    # Dy  = sys.Dy

    Nu, Nx, Ny = size(sys)

    @fastmath @inbounds for j in 1:Ny
        @inbounds for i in 1:Nx
            @inbounds @simd for a in 1:Nu
                u      = uu[a]
                u3     = u * u * u
                u4     = u * u3

                phi    = phiGF[a,i,j]
                phid   = phidGF[a,i,j]
                A      = AGF[a,i,j]

                phi_u  = Du_phi[a,i,j]
                # phid_u = Du_phid[a,i,j]

                phiGF_t[a,i,j] =  3//2 * (phi - phid) / u + phi_u / 2 + 3//2 * u3 * A * phi +
                    A * u4 * phi_u / 2
            end
        end
    end

    nothing
end


function sync_bulkevolved!(bulkevols_t, bulkconstrains,
                           systems::SystemPartition, evoleq)
    Nsys = length(systems)

    @inbounds for i in 1:Nsys-1
        sync_bulkevolved!(bulkevols_t[i], bulkevols_t[i+1], bulkconstrains[i+1],
                          systems[i], systems[i+1], evoleq)
    end

    nothing
end

function sync_bulkevolved!(bulkevol1_t::BulkEvolved, bulkevol2_t::BulkEvolved,
                           bulkconstrain2::BulkConstrained,
                           sys1::System, sys2::System, evoleq::AffineNull)
    u0 = sys2.ucoord[1]

    _, Nx, Ny = size(sys2)

    phiGF1_t = getphi(bulkevol1_t)
    phiGF2_t = getphi(bulkevol2_t)
    AGF      = getA(bulkconstrain2)

    @fastmath @inbounds for j in 1:Ny
        @inbounds for i in 1:Nx
            A = AGF[1,i,j]

            # characteristic speed
            c = u0 * u0 * A/2

            # if c > 0, mode is entering grid1 from grid2;
            # if c < 0, mode is entering grid2 from grid1.
            # we assume here that grids merely touch at the interface
            if c > 0
                phiGF1_t[end,i,j] = phiGF2_t[1,i,j]
            elseif c < 0
                phiGF2_t[1,i,j] = phiGF1_t[end,i,j]
            end
        end
    end

    nothing
end

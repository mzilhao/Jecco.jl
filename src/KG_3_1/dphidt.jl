
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

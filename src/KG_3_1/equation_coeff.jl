
## TODO
# Need signature Vf(phi::Real)
function Vf end
function Vfp end

# assuming
# (A d_uu + B d_u + C Id) f = -S

function phig1_eq_coeff!(ABCS::Vector, vars::AllVars)
    u  = vars.u
    u2 = u * u
    u3 = u2 * u
    u4 = u2 * u2
    u5 = u3 * u2
    u6 = u4 * u2
    u9 = u6 * u3

    phi  = vars.phi_d0
    phi3 = phi * phi * phi
    phi4 = phi * phi3

    ABCS[1] = 0
    ABCS[2] = 9.0 * u
    ABCS[3] = 4.5

    ABCS[4] = -4.5 * phi - 27 * u4 * phi * vars.Sd_d0 -
        4 * u6 * phi3 * Vf(phi) - u9 * phi4 * Vfp(phi) +
        3 * u2 * (vars.phi_dxx + vars.phi_dyy) -
        4.5 * u * vars.phi_du - 9 * u5 * vars.Sd_d0 * vars.phi_du

    nothing
end

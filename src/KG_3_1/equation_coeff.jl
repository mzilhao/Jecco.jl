
# assuming
# (A d_uu + B d_u + C Id) f = -Source

function phid_eq_coeff!(ABCS::Vector, vars)
    (potential, u, Sd, phi, phi_u, phi_xx, phi_yy) = vars

    u2 = u * u
    u3 = u2 * u
    u4 = u2 * u2
    u5 = u3 * u2
    u6 = u4 * u2
    u9 = u6 * u3

    phi3 = phi * phi * phi
    phi4 = phi * phi3

    ABCS[1] = 0
    ABCS[2] = 9 * u
    ABCS[3] = 9//2

    ABCS[4] = -9//2 * phi - 27 * u4 * phi * Sd -
        4 * u6 * phi3 * UU(phi,potential) - u9 * phi4 * UUp(phi,potential) +
        3 * u2 * (phi_xx + phi_yy) -
        9//2 * u * phi_u - 9 * u5 * Sd * phi_u

    nothing
end

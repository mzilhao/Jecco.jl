
@inline S_inner_to_outer(S_in, u, xi, phi0) =
    1/u + xi - u * phi0*phi0/3 * (1 - u * xi) + u*u*u * S_in

@inline S_u_inner_to_outer(S_u_in, S_in, u, xi, phi0) =
    -1/(u*u) - phi0*phi0/3 * (1 - u * xi) + u * xi * phi0*phi0/3 +
    3 * u*u * S_in + u*u*u * S_u_in

@inline F_inner_to_outer(F_in, u) = u*u * F_in

@inline F_u_inner_to_outer(F_u_in, F_in, u) = 2 * u * F_in + u*u * F_u_in

@inline Sd_inner_to_outer(Sd_in, u, xi, phi0) =
    1/(2*u*u) + xi/u + xi*xi/2 - phi0*phi0/6 + u*u * Sd_in

@inline Bd_inner_to_outer(Bd_in, u) = u*u*u * Bd_in

@inline phid_inner_to_outer(phid_in, u, phi0) =
    -phi0/2 + u*u * phi0*phi0*phi0 * phid_in

@inline A_inner_to_outer(A_in, u, xi, phi0) =
    1/(u*u) + 2*xi/u + xi*xi - 2/3 * phi0*phi0 + u*u * A_in

@inline A_u_inner_to_outer(A_u_in, A_in, u, xi, phi0) =
    -2/(u*u*u) - 2*xi/(u*u) + 2*u * A_in + u*u * A_u_in

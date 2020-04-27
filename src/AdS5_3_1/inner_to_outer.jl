
@inline S_inner_to_outer(S_in, u, xi, phi0) =
    1/u + xi - u * phi0 * phi0 / 3 * ( 1 - u * xi) + u * u * u * S_in

@inline S_u_inner_to_outer(S_u_in, S_in, u, xi, phi0) =
    -1/(u*u) + xi - phi0 * phi0 / 3 * ( 1 - u * xi) + u * xi * phi0 * phi0 / 3 +
    3 * u * u * S_in + u * u * u * S_u_in

@inline F_inner_to_outer(F_in, u) = u * u * F_in

@inline F_u_inner_to_outer(F_u_in, F_in, u) = 2 * u * F_in + u * u * F_u_in

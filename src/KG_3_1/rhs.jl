
function setup_rhs(phi::Array, sys)

    a4 = -ones2D(sys)
    # boundary = BoundaryVars(a4)

    bulk = BulkVars(phi)
    boundary = bulk[1,:,:]

    solve_nested_g1! = nested_g1(sys)

    function (df, f, sys, t)
        bulk.phi .= f

        # FIXME
        boundary.Sd   .= 0.5 * a4
        boundary.phid .= bulk.phi[1,:,:] # phi2
        boundary.A    .= a4

        solve_nested_g1!(bulk, boundary)

        rhs_g1!(df, bulk, sys)
        nothing
    end
end

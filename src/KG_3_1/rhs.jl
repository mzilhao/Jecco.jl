
function setup_rhs(phi::Array, sys)

    a4 = -ones2D(sys)
    boundary = BoundaryVars(a4)

    bulk = BulkVars(phi)

    solve_nested_g1! = nested_g1(sys)

    function (df, f, sys, t)
        # TODO: check if it's better to use "=" with a mutable struct
        bulk.phi .= f
        solve_nested_g1!(bulk, boundary)

        rhs_g1!(df, bulk, sys)
        nothing
    end
end

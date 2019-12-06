
function setup_rhs(phi::Array{<:Number,N}, sys::System) where {N}

    a4 = -ones2D(sys)
    boundary = BoundaryVars(a4)

    bulk = BulkVars(phi)
    BC = bulk[1,:,:]

    nested = Nested(sys)

    function (df, f, sys, t)
        bulk.phi .= f

        solve_nested_g1!(bulk, BC, boundary, nested)

        df .= bulk.dphidt
        nothing
    end
end

function setup_rhs(phis::Vector, systems::Vector)

    a4 = -ones2D(systems[1])
    boundary = BoundaryVars(a4)

    bulks = BulkVars(phis)
    phis_slice  = [phi[1,:,:] for phi in phis]
    BCs  = BulkVars(phis_slice)

    Nsys    = length(systems)
    nesteds = Nested(systems)

    function (df::ArrayPartition, f::ArrayPartition, systems, t)
        for i in 1:Nsys
            bulks[i].phi .= f.x[i]
        end

        solve_nested_g1!(bulks, BCs, boundary, nesteds)

        for i in 1:Nsys
            df.x[i] .= bulks[i].dphidt
        end

        nothing
    end
end

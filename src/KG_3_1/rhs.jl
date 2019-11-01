
function setup_rhs(phi::Array{<:Number,N}, sys::System) where {N}

    a4 = -ones2D(sys)
    # boundary = BoundaryVars(a4)

    bulk = BulkVars(phi)
    boundary = bulk[1,:,:]

    nested = Nested(sys)

    function (df, f, sys, t)
        bulk.phi .= f

        # u=0 boundary
        boundary.Sd   .= 0.5 * a4
        boundary.phid .= bulk.phi[1,:,:] # phi2
        boundary.A    .= a4

        nested_g1!(nested, bulk, boundary)

        df .= bulk.dphidt
        nothing
    end
end

function setup_rhs(phis::Vector, systems::Vector)

    a4 = -ones2D(systems[1])

    bulks = BulkVars(phis)
    phis_slice  = [phi[1,:,:] for phi in phis]
    boundaries  = BulkVars(phis_slice)

    Nsys    = length(systems)
    nesteds = Nested(systems)

    function (df::ArrayPartition, f::ArrayPartition, sys, t)
        for i in 1:Nsys
            bulks[i].phi .= f.x[i]
        end

        # u=0 boundary
        boundaries[1].Sd   .= 0.5 * a4
        boundaries[1].phid .= bulks[1].phi[1,:,:] # phi2
        boundaries[1].A    .= a4

        for i in 1:Nsys-1
            nested_g1!(nesteds[i], bulks[i], boundaries[i])
            boundaries[i+1] = bulks[i][end,:,:]
        end
        nested_g1!(nesteds[Nsys], bulks[Nsys], boundaries[Nsys])

        # FIXME: sync boundary points

        for i in 1:Nsys
            df.x[i] .= bulks[i].dphidt
        end

        nothing
    end
end

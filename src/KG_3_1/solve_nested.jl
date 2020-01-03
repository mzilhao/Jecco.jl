
import Base.Threads.@threads
import Base.Threads.@spawn
using LinearAlgebra

function solve_lin_system!(sol, A_mat, b_vec)
    A_fact = lu!(A_mat)
    ldiv!(A_fact, b_vec)
    sol .= b_vec
    nothing
end

struct Aux{T<:Real}
    A_mat   :: Matrix{T}
    b_vec   :: Vector{T}
    ABCS    :: Vector{T}
    vars    :: AllVars{T}

    function Aux{T}(N::Int) where {T<:Real}
        A_mat  = zeros(T, N, N)
        b_vec  = zeros(T, N)
        ABCS   = zeros(T, 4)
        vars   = AllVars{T}()

        new(A_mat, b_vec, ABCS, vars)
    end
end

struct Nested{S,D,T<:Real}
    sys     :: S
    uu      :: Vector{T}
    xx      :: Vector{T}
    yy      :: Vector{T}
    Du_phi  :: D
    Dxx_phi :: D
    Dyy_phi :: D
    Du_phid :: D
    aux_acc :: Vector{Aux{T}}
end
function Nested(sys::System)
    coords = sys.coords

    uu, xx, yy = Vivi.xx(coords)
    Nu = length(uu)
    Nx = length(xx)
    Ny = length(yy)

    Du_phi    = zeros(Nu, Nx, Ny)
    Dxx_phi   = zeros(Nu, Nx, Ny)
    Dyy_phi   = zeros(Nu, Nx, Ny)
    Du_phid   = zeros(Nu, Nx, Ny)

    nt = Threads.nthreads()
    # pre-allocate thread-local aux quantities
    aux_acc = [Aux{eltype(uu)}(Nu) for _ in 1:nt]

    Nested{typeof(sys),typeof(Du_phi), eltype(uu)}(sys,uu, xx, yy, Du_phi, Dxx_phi,
                                                      Dyy_phi, Du_phid, aux_acc)
end

Nested(systems::Vector) = [Nested(sys) for sys in systems]


function solve_nested_g1!(bulk::BulkVars, BC::BulkVars, nested::Nested)
    sys  = nested.sys
    uu   = nested.uu
    xx   = nested.xx
    yy   = nested.yy

    Du_phi  = nested.Du_phi
    Dxx_phi = nested.Dxx_phi
    Dyy_phi = nested.Dyy_phi
    Du_phid = nested.Du_phid

    aux_acc = nested.aux_acc

    Du  = sys.Du
    Duu = sys.Duu
    # Dx  = sys.Dx
    Dxx = sys.Dxx
    # Dy  = sys.Dy
    Dyy = sys.Dyy

    t1 = @spawn mul!(Du_phi, Du, bulk.phi)
    t2 = @spawn mul!(Dxx_phi, Dxx, bulk.phi)
    t3 = @spawn mul!(Dyy_phi, Dyy, bulk.phi)

    # set Sd
    @fastmath @inbounds for j in eachindex(yy)
        @inbounds for i in eachindex(xx)
            @inbounds @simd for a in eachindex(uu)
                bulk.Sd[a,i,j] = BC.Sd[i,j]
            end
        end
    end

    wait(t1)
    wait(t2)
    wait(t3)

    # solve for phidg1

    @fastmath @inbounds @threads for j in eachindex(yy)
        @inbounds for i in eachindex(xx)
            id  = Threads.threadid()
            aux = aux_acc[id]

            @inbounds @simd for a in eachindex(uu)
                aux.vars.u       = uu[a]
                aux.vars.Sd_d0   = bulk.Sd[a,i,j]
                aux.vars.phi_d0  = bulk.phi[a,i,j]
                aux.vars.phi_du  = Du_phi[a,i,j]
                aux.vars.phi_dxx = Dxx_phi[a,i,j]
                aux.vars.phi_dyy = Dyy_phi[a,i,j]

                phig1_eq_coeff!(aux.ABCS, aux.vars)

                aux.b_vec[a]     = -aux.ABCS[4]

                # TODO: replace Duu.D and Du.D with something more general and portable
                @inbounds @simd for aa in eachindex(uu)
                    aux.A_mat[a,aa] = aux.ABCS[1] * Duu.D[a,aa] + aux.ABCS[2] * Du.D[a,aa]
                end
                aux.A_mat[a,a] += aux.ABCS[3]
            end

            # boundary condition
            aux.b_vec[1]    = BC.phid[i,j]
            aux.A_mat[1,:] .= 0.0
            aux.A_mat[1,1]  = 1.0

            sol = view(bulk.phid, :, i, j)
            solve_lin_system!(sol, aux.A_mat, aux.b_vec)
        end
    end

    # set Ag1
    @fastmath @inbounds for j in eachindex(yy)
        @inbounds for i in eachindex(xx)
            @inbounds @simd for a in eachindex(uu)
                bulk.A[a,i,j] = BC.A[i,j]
            end
        end
    end


    # finally compute dphidt_g1

    mul!(Du_phid, Du, bulk.phid)

    @fastmath @inbounds @threads for j in eachindex(yy)
        @inbounds for i in eachindex(xx)
            id  = Threads.threadid()
            aux = aux_acc[id]
            @inbounds @simd for a in eachindex(uu)
                aux.vars.u       = uu[a]

                aux.vars.phi_d0  = bulk.phi[a,i,j]
                aux.vars.phid_d0 = bulk.phid[a,i,j]
                aux.vars.A_d0    = bulk.A[a,i,j]

                aux.vars.phi_du  = Du_phi[a,i,j]
                aux.vars.phid_du = Du_phid[a,i,j]

                if aux.vars.u > 1.e-9
                    bulk.dphidt[a,i,j]  = dphig1dt(aux.vars)
                else
                    bulk.dphidt[a,i,j]  = dphig1dt_u0(aux.vars)
                end
            end
        end
    end

    nothing
end


function solve_nested_g1!(bulk::BulkVars, BC::BulkVars, boundary::BoundaryVars,
                          nested::Nested)
    # u=0 boundary
    BC.Sd   .= 0.5 * boundary.a4
    BC.phid .= bulk.phi[1,:,:] # phi2
    BC.A    .= boundary.a4

    solve_nested_g1!(bulk, BC, nested)

    nothing
end

function solve_nested_g1!(bulks::Vector, BCs::Vector, boundary::BoundaryVars,
                          nesteds::Vector)
    Nsys = length(nesteds)

    # u=0 boundary
    BCs[1].Sd   .= 0.5 * boundary.a4
    BCs[1].phid .= bulks[1].phi[1,:,:] # phi2
    BCs[1].A    .= boundary.a4

    for i in 1:Nsys-1
        solve_nested_g1!(bulks[i], BCs[i], nesteds[i])
        BCs[i+1] = bulks[i][end,:,:]
    end
    solve_nested_g1!(bulks[Nsys], BCs[Nsys], nesteds[Nsys])

    # sync boundary points. note: in a more general situation we may need to
    # check the characteristic speeds (in this case we just know where the
    # horizon is)
    for i in 1:Nsys-1
        bulks[i].dphidt[end,:,:] .= bulks[i+1].dphidt[1,:,:]
    end

    nothing
end


function solve_nested_g1(phi::Array{<:Number,N}, sys::System) where {N}
    a4 = -ones2D(sys)
    boundary = BoundaryVars(a4)

    bulk = BulkVars(phi)
    BC = bulk[1,:,:]

    nested = Nested(sys)

    solve_nested_g1!(bulk, BC, boundary, nested)
    bulk
end

function solve_nested_g1(phis::Vector, systems::Vector)
    a4 = -ones2D(systems[1])
    boundary = BoundaryVars(a4)

    bulks = BulkVars(phis)
    phis_slice  = [phi[1,:,:] for phi in phis]
    BCs  = BulkVars(phis_slice)

    Nsys    = length(systems)
    nesteds = Nested(systems)

    solve_nested_g1!(bulks, BCs, boundary, nesteds)
    bulks
end

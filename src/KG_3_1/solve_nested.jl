
function solve_lin_system!(A_mat, b_vec)
    # passing Val(false) to the second argument turns off pivoting. it seems to
    # improve speed for the small matrices that we typically consider. we can
    # revisit this (or make it a parameter) if needed.
    A_fact = lu!(A_mat, Val(false))
    ldiv!(A_fact, b_vec)        # b_vec is overwritten to store the result
    nothing
end

struct Aux{T<:Real}
    A_mat   :: Matrix{T}
    b_vec   :: Vector{T}
    ABCS    :: Vector{T}
    function Aux{T}(N::Int) where {T<:Real}
        A_mat  = zeros(T, N, N)
        b_vec  = zeros(T, N)
        ABCS   = zeros(T, 4)
        new(A_mat, b_vec, ABCS)
    end
end
function Aux{T}(sys::System) where {T}
    Nu, Nx, Ny = size(sys)
    nt = Threads.nthreads()
    # pre-allocate thread-local aux quantities
    [Aux{T}(Nu) for _ in 1:nt]
end

struct BC{T}
    phid :: Array{T,2}
    Sd   :: Array{T,2}
    A    :: Array{T,2}
end
function BC{T}(Nx::Int, Ny::Int) where {T<:Real}
    phid = Array{T}(undef, Nx, Ny)
    Sd   = Array{T}(undef, Nx, Ny)
    A    = Array{T}(undef, Nx, Ny)
    BC{T}(phid, Sd, A)
end

function BCs(systems::SystemPartition)
    sys1  = systems[1]
    T     = Jecco.coord_eltype(sys1.ucoord)
    _, Nx, Ny = size(sys1)

    [BC{T}(Nx, Ny) for sys in systems]
end


function solve_phid!(bulkconstrain::BulkConstrained, bulkevol::BulkEvolved, bc::BC,
                     deriv::BulkDeriv, aux_acc, sys::System, evoleq::AffineNull)
    phiGF  = getphi(bulkevol)
    phidGF = getphid(bulkconstrain)
    SdGF   = getSd(bulkconstrain)
    AGF    = getA(bulkconstrain)

    Du_phi  = deriv.Du_phi

    Nu, Nx, Ny = size(sys)

    Du  = sys.Du
    Duu = sys.Duu
    Dx  = sys.Dx
    Dxx = sys.Dxx
    Dy  = sys.Dy
    Dyy = sys.Dyy

    potential = evoleq.potential

    @fastmath @inbounds @threads for j in 1:Ny
        @inbounds for i in 1:Nx
            id  = Threads.threadid()
            aux = aux_acc[id]

            @inbounds @simd for a in 1:Nu
                u      = sys.ucoord[a]

                Sd     = SdGF[a,i,j]
                phi    = phiGF[a,i,j]
                phi_u  = Du_phi[a,i,j]
                phi_xx = Dxx(phiGF,a,i,j)
                phi_yy = Dyy(phiGF,a,i,j)

                vars = (potential, u, Sd, phi, phi_u, phi_xx, phi_yy)

                phid_eq_coeff!(aux.ABCS, vars)

                aux.b_vec[a]   = -aux.ABCS[4]
                @inbounds @simd for aa in 1:Nu
                    aux.A_mat[a,aa] = aux.ABCS[1] * Duu[a,aa] + aux.ABCS[2] * Du[a,aa]
                end
                aux.A_mat[a,a] += aux.ABCS[3]
            end

            # BC (first order equation)

            aux.b_vec[1]    = bc.phid[i,j]
            aux.A_mat[1,:] .= 0.0
            aux.A_mat[1,1]  = 1.0

            solve_lin_system!(aux.A_mat, aux.b_vec)

            @inbounds @simd for aa in 1:Nu
                phidGF[aa,i,j] = aux.b_vec[aa]
            end

        end
    end

    nothing
end


function solve_nested!(bulkconstrain::BulkConstrained, bulkevol::BulkEvolved, bc::BC,
                       deriv::BulkDeriv, aux_acc, sys::System, evoleq::AffineNull)
    phiGF  = getphi(bulkevol)
    phidGF = getphid(bulkconstrain)
    SdGF   = getSd(bulkconstrain)
    AGF    = getA(bulkconstrain)

    Du_phi  = deriv.Du_phi
    Du_phid = deriv.Du_phid

    Nu, Nx, Ny = size(sys)

    Du  = sys.Du
    Duu = sys.Duu
    Dx  = sys.Dx
    Dxx = sys.Dxx
    Dy  = sys.Dy
    Dyy = sys.Dyy


    # set Sd and A
    @fastmath @inbounds for j in 1:Ny
        @inbounds for i in 1:Nx
            @inbounds @simd for a in 1:Nu
                SdGF[a,i,j] = bc.Sd[i,j]
                AGF[a,i,j]  = bc.A[i,j]
            end
        end
    end

    # solve for phid
    solve_phid!(bulkconstrain, bulkevol, bc, deriv, aux_acc, sys, evoleq)

    # take u-derivatives of phid
    mul!(Du_phid,  Du,  phidGF)

    nothing
end


function syncBCs!(bc::BC, bulk::BulkConstrained, deriv::BulkDeriv, sys::System)
    _, Nx, Ny = size(sys)

    phidGF = getphid(bulk)
    SdGF   = getSd(bulk)
    AGF    = getA(bulk)

    # we are here assuming that the inner and outer grids merely touch at the
    # interface, so we pass the values at this point without any interpolation
    @fastmath @inbounds @threads for j in 1:Ny
        @inbounds @simd for i in 1:Nx
            bc.Sd[i,j]   = SdGF[end,i,j]
            bc.phid[i,j] = phidGF[end,i,j]
            bc.A[i,j]    = AGF[end,i,j]
        end
    end

    nothing
end

function set_innerBCs!(bc::BC, bulk::BulkEvolved, boundary::Boundary,
                       deriv::BulkDeriv, sys::System)
    _, Nx, Ny = size(sys)

    phiGF = getphi(bulk)
    a4GF  = geta4(boundary)

    @fastmath @inbounds @threads for j in 1:Ny
        @inbounds @simd for i in 1:Nx
            bc.phid[i,j] = phiGF[1,i,j]
            bc.Sd[i,j]   = 0.5 * a4GF[1,i,j]
            bc.A[i,j]    = a4GF[1,i,j]
        end
    end

    nothing
end


function solve_nesteds!(bulkconstrains, bulkevols, boundary::Boundary,
                        bcs, derivs, aux_accs,
                        systems::SystemPartition, evoleq::AffineNull)
    Nsys = length(systems)

    # take u-derivatives of the bulkevols functions
    @inbounds for i in 1:Nsys
        Du_phi = derivs[i].Du_phi
        Du     = systems[i].Du
        phiGF  = getphi(bulkevols[i])
        mul!(Du_phi, Du, phiGF)
    end

    set_innerBCs!(bcs[1], bulkevols[1], boundary, derivs[1], systems[1])

    @inbounds for i in 1:Nsys-1
        solve_nested!(bulkconstrains[i], bulkevols[i], bcs[i],
                      derivs[i], aux_accs[i], systems[i], evoleq)
        syncBCs!(bcs[i+1], bulkconstrains[i], derivs[i], systems[i+1])
    end
    solve_nested!(bulkconstrains[Nsys], bulkevols[Nsys], bcs[Nsys],
                  derivs[Nsys], aux_accs[Nsys], systems[Nsys], evoleq)

    nothing
end


struct Nested{S,C,D,B,A}
    systems        :: S
    bulkconstrains :: C
    derivs         :: D
    bcs            :: B
    aux_accs       :: A
end

function Nested(systems::SystemPartition, bulkconstrains::BulkPartition,
                derivs::BulkPartition)
    T        = Jecco.coord_eltype(systems[1].ucoord)
    bcs      = BCs(systems)
    aux_accs = [Aux{T}(sys) for sys in systems]

    Nested{typeof(systems),typeof(bulkconstrains),typeof(derivs),
           typeof(bcs),typeof(aux_accs)}(systems, bulkconstrains, derivs, bcs,
                                         aux_accs)
end

function Nested(systems::SystemPartition, bulkconstrains::BulkPartition)
    T        = Jecco.coord_eltype(systems[1].ucoord)
    bcs      = BCs(systems)
    aux_accs = [Aux{T}(sys) for sys in systems]
    derivs   = BulkDerivPartition(systems)

    Nested(systems, bulkconstrains, derivs)
end

function Nested(systems::SystemPartition)
    T        = Jecco.coord_eltype(systems[1].ucoord)
    bcs      = BCs(systems)
    aux_accs = [Aux{T}(sys) for sys in systems]

    bulkconstrains = BulkConstrainedPartition(systems)
    derivs         = BulkDerivPartition(systems)

    Nested(systems, bulkconstrains, derivs)
end

function (nested::Nested)(bulkevols::BulkPartition, boundary::Boundary,
                          evoleq::EvolutionEquations)
    solve_nesteds!(nested.bulkconstrains, bulkevols, boundary, nested.bcs,
                   nested.derivs, nested.aux_accs, nested.systems, evoleq)
end

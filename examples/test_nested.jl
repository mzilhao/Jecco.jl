
using Jecco
using Jecco.AdS5_3_1

import Base.Threads.@threads
import Base.Threads.@spawn
using LinearAlgebra

# par_base = ParamBase(
#     which_potential = "square",
# )

function IDtest0(sys::AbstractSystem{Outer})
    Nu, Nx, Ny = size(sys)
    ucoord = sys.ucoord
    xcoord = sys.xcoord
    ycoord = sys.ycoord

    B1  = zeros(Nu, Nx, Ny)
    B2  = zeros(Nu, Nx, Ny)
    G   = zeros(Nu, Nx, Ny)
    phi = zeros(Nu, Nx, Ny)

    b14 = 0.01
    b24 = 0.02

    phi0 = 1.0
    phi2 = 0.01

    for j in 1:Ny
        for i in 1:Nx
            for a in 1:Nu
                u = ucoord[a]
                x = xcoord[i]
                y = ycoord[j]
                B1[a,i,j]  = u^4 * b14
                B2[a,i,j]  = u^4 * b24
                phi[a,i,j] = phi0 * u + phi2 * u^3
                G[a,i,j]   = 0.0
            end
        end
    end

    BulkVars(B1, B2, G, phi)
end

function IDtest0(sys::AbstractSystem{Inner})
    Nu, Nx, Ny = size(sys)
    ucoord = sys.ucoord
    xcoord = sys.xcoord
    ycoord = sys.ycoord

    B1  = zeros(Nu, Nx, Ny)
    B2  = zeros(Nu, Nx, Ny)
    G   = zeros(Nu, Nx, Ny)
    phi = zeros(Nu, Nx, Ny)

    b14 = 0.01
    b24 = 0.02

    phi0 = 1.0
    phi2 = 0.01

    for j in 1:Ny
        for i in 1:Nx
            for a in 1:Nu
                u = ucoord[a]
                x = xcoord[i]
                y = ycoord[j]
                B1[a,i,j]  = b14
                B2[a,i,j]  = b24
                phi[a,i,j] = phi2
                G[a,i,j]   = 0.0
            end
        end
    end

    BulkVars(B1, B2, G, phi)
end

IDtest0(systems::Vector{T}) where {T<:System} = [IDtest0(sys) for sys in systems]


par_grid = ParamGrid(
    x_min            = -5.0,
    x_max            =  5.0,
    x_nodes          =  128,
    y_min            = -5.0,
    y_max            =  5.0,
    y_nodes          =  128,
    u_outer_min      =  0.1,
    u_outer_max      =  1.0,
    # u_outer_domains  =  1,
    # u_outer_nodes    =  64,
    u_outer_domains  =  2,
    u_outer_nodes    =  16,
    u_inner_nodes    =  12,
)

# TODO: make parameter
kappa = 1.0

systems = Jecco.AdS5_3_1.create_systems(par_grid)

sys = systems[2]

Nu, Nx, Ny = size(sys)
u = sys.ucoord[:]
# x = sys.xcoord[:]
# y = sys.ycoord[:]


# Note: it's important that xi is defined on a 1*Nx*Ny grid, rather than a Nx*Ny
# one, so that the same Dx and Dy differential operators defined for the bulk
# quantities can also straightforwardly apply on xi. Remember that the axis
# along which the operator applies is defined on the operator itself. So, by
# defining things this way, the Dx operator (which acts along the 2nd index)
# will also do the correct thing when acting on xi.
xi = zeros(Float64, 1, Nx, Ny)
gauge = GaugeVars(xi, kappa)

nested = Jecco.AdS5_3_1.Nested(sys)

bulks = IDtest0(systems)

nesteds = [Jecco.AdS5_3_1.Nested(sys) for sys in systems]
BCs     = [BulkVars(Nx, Ny) for sys in systems]
dBCs    = [BulkVars(Nx, Ny) for sys in systems]

u0 = u[1]

fx2_0 = 0.02
fy2_0 = 0.1

BCs[2].S  .= 1.0/u0
dBCs[2].S .= -1.0/(u0*u0)

BCs[2].Fx .= fx2_0 * u0 * u0
BCs[2].Fy .= fy2_0 * u0 * u0

dBCs[2].Fx .= 2 * fx2_0 * u0
dBCs[2].Fy .= 2 * fy2_0 * u0

BCs[2].Sd .= 0.5/(u0*u0)

BCs[2].B2d .= -2.0 * u0*u0*u0 * 0.02
BCs[2].B1d .= -2.0 * u0*u0*u0 * 0.01

BCs[2].Gd   .= 0.0
BCs[2].phid .= -0.5 + u0*u0 * ( 1.0/3.0 - 1.5 * 0.01 )

BCs[2].A  .= 1.0/(u0*u0)
dBCs[2].A .= -2.0/(u0*u0*u0)

Jecco.AdS5_3_1.solve_nested!(bulks, BCs, dBCs, gauge, nesteds)

bulk = bulks[2]
BC   = BCs[2]
dBC  = dBCs[2]
nested = nesteds[2]

# Jecco.AdS5_3_1.solve_S_outer!(bulk, BC, dBC, nested)
# Jecco.AdS5_3_1.solve_A_outer!(bulk, BC, dBC, nested)

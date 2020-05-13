
using Jecco
using Jecco.AdS5_3_1

import Base.Threads.@threads
import Base.Threads.@spawn
using LinearAlgebra

function IDtest0(sys::System{Outer})
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

    BulkVars(sys.gridtype, B1, B2, G, phi)
end

function IDtest0(sys::System{Inner})
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

    BulkVars(sys.gridtype, B1, B2, G, phi)
end

IDtest0(systems::Vector{T}) where {T<:System} = [IDtest0(sys) for sys in systems]


par_grid = Grid3D(
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

# par_evol = ParamEvol(
#     ODE_method = "AB3",
#     # ODE_method = "RK4",
#     dt      = 0.008,
#     tmax    = 1.0,
# )



potential   = ZeroPotential()
phi0        = 0.0

base  = BaseVars(potential, phi0)


kappa = 1.0


systems = Jecco.AdS5_3_1.Systems(par_grid)

sys = systems[1]

_, Nx, Ny = size(sys)
# u = sys.ucoord[:]
# x = sys.xcoord[:]
# y = sys.ycoord[:]


# Note: it's important that xi, a4 and f2 are defined on a 1*Nx*Ny grid, rather
# than a Nx*Ny one, so that the same Dx and Dy differential operators defined
# for the bulk quantities can also straightforwardly apply on them. Remember that
# the axis along which the operator applies is defined on the operator itself.
# So, by defining things this way, the Dx operator (which acts along the 2nd
# index) will also do the correct thing when acting on xi, a4 or f2.
xi = zeros(Float64, 1, Nx, Ny)
gauge = GaugeVars(xi, kappa)


fx2_0 = 0.02
fy2_0 = 0.1


f2x = fx2_0 * ones(Float64, 1, Nx, Ny)
f2y = fy2_0 * ones(Float64, 1, Nx, Ny)
a4  = -ones(Float64, 1, Nx, Ny)
boundary = BoundaryVars(a4, f2x, f2y)

bulks = IDtest0(systems)

solve_nested = Jecco.AdS5_3_1.nested_solver(base, systems)
solve_nested(bulks, boundary, gauge)



# nesteds = [Jecco.AdS5_3_1.Nested(sys) for sys in systems]
# BCs     = [BulkVars(sys.gridtype, Float64, Nx, Ny) for sys in systems]
# dBCs    = [BulkVars(sys.gridtype, Float64, Nx, Ny) for sys in systems]


# Jecco.AdS5_3_1.set_innerBCs!(BCs[1], dBCs[1], bulks[1], boundary, gauge, base, nesteds[1])


# Jecco.AdS5_3_1.set_outerBCs!(BCs[2], dBCs[2], bulks[1], gauge, base, nesteds[1])


# u0 = u[1]


# BCs[2].Sd .= 0.5/(u0*u0)

# BCs[2].B2d .= -2.0 * u0*u0*u0 * 0.02
# BCs[2].B1d .= -2.0 * u0*u0*u0 * 0.01

# BCs[2].Gd   .= 0.0
# BCs[2].phid .= -0.5 + u0*u0 * ( 1.0/3.0 - 1.5 * 0.01 )

# BCs[2].A  .= 1.0/(u0*u0)
# dBCs[2].A .= -2.0/(u0*u0*u0)


i = 2

u = systems[i].ucoord[:]

bulk = bulks[i]
# BC   = BCs[i]
# dBC  = dBCs[i]
# nested = nesteds[i]

bulk.A[:,1,1]

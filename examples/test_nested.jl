
using Jecco
using Jecco.AdS5_3_1

import Base.Threads.@threads
import Base.Threads.@spawn
using LinearAlgebra

# par_base = ParamBase(
#     which_potential = "square",
# )

function IDtest0(sys::System)
    Nu, Nx, Ny = size(sys.grid)

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
                u,x,y = sys.grid[a,i,j]
                B1[a,i,j]  = u^4 * b14
                B2[a,i,j]  = u^4 * b24
                phi[a,i,j] = phi0 * u + phi2 * u^3
                G[a,i,j]   = 0.0
            end
        end
    end

    BulkVars(B1, B2, G, phi)
end


par_grid = ParamGrid(
    xmin        = -5.0,
    xmax        =  5.0,
    xnodes      =  128,
    ymin        = -5.0,
    ymax        =  5.0,
    ynodes      =  128,
    # umin        =  0.01,
    umin        =  0.1,
    umax        =  1.0,
    udomains    =  1,
    # unodes      =  16,
    unodes      =  32,
    # unodes      =  64,
)


systems = Jecco.AdS5_3_1.create_systems(par_grid)

sys = systems[1]

# Nu, Nx, Ny = size(sys.grid)
Nu_, Nx, Ny = size(sys.grid)
u, x, y = sys.grid[:]

nested = Jecco.AdS5_3_1.Nested(sys)

bulk = IDtest0(sys)

BC  = BulkVars(Nx, Ny)
dBC = BulkVars(Nx, Ny)

u0 = u[1]

fx2_0 = 0.02
fy2_0 = 0.1

BC.S  .= 1.0/u0
dBC.S .= -1.0/(u0*u0)

BC.Fx .= fx2_0 * u0 * u0
BC.Fy .= fy2_0 * u0 * u0

dBC.Fx .= 2 * fx2_0 * u0
dBC.Fy .= 2 * fy2_0 * u0

BC.Sd .= 0.5/(u0*u0)

BC.B2d .= -2.0 * u0*u0*u0 * 0.02
BC.B1d .= -2.0 * u0*u0*u0 * 0.01

BC.Gd   .= 0.0
BC.phid .= -0.5 + u0*u0 * ( 1.0/3.0 - 1.5 * 0.01 )

BC.A  .= 1.0/(u0*u0)
dBC.A .= -2.0/(u0*u0*u0)

Jecco.AdS5_3_1.solve_nested_outer!(bulk, BC, dBC, nested)


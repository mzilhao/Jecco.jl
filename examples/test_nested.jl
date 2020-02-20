
using Jecco
using Jecco.AdS_3_1

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

    for j in 1:Ny
        for i in 1:Nx
            for a in 1:Nu
                u,x,y = sys.grid[a,i,j]
                B1[a,i,j]  = u^4 * b14
                B2[a,i,j]  = u^4 * b24
                phi[a,i,j] = 0.0
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
    umin        =  0.1,
    umax        =  1.0,
    udomains    =  1,
    unodes      =  16,
)


systems = Jecco.AdS_3_1.create_systems(par_grid)

sys = systems[1]

Nu, Nx, Ny = size(sys.grid)

nested = Jecco.AdS_3_1.Nested(sys)

bulk = IDtest0(sys)

Jecco.AdS_3_1.solve_nested_outer!(bulk, bulk, nested)

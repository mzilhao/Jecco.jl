using Jecco
using HDF5
#using Plots
#gr()
include("./Fourier_modes.jl")


dirname = "/Users/apple/Documents/Jecco.jl/data/spinodal_e_0.9_L_10/"
#Specify the initial iteration, the final and the step in the jump.
#Remember to set same output frequency for both boundary and gauge files
initial_it = 0
step_it = 500
final_it = 89500

#Phi0 and inverse of phiM squared
phi0 = 1.0
phiM2 = -1.0

To_hdf5 = "no"
To_mathematica = "yes"
Fourier_modes = "yes"
order_x, order_y = (15,0)
integration_method = "trapezoidal"
#Read_tensor = "yes"

files = OpenPMDTimeSeries(dirname,prefix="boundary")
files_xi = OpenPMDTimeSeries(dirname,prefix="gauge")


#Computation of the energy stress tensor. Each component is going to be a 3D array, where the index are [time,x,y].
u, x, y = get_field(files,it=initial_it,field="a4")[2][:]
n_t = Int(floor((final_it-initial_it)/step_it))+1
#TODO: Add the last point
#=append!(x, x[end]+(x[2]-x[1]))
append!(y, y[end]+(y[2]-y[1]))
=#
n_x, n_y = length(x), length(y)


e, Jx, Jy = zeros(n_t,n_x,n_y), zeros(n_t,n_x,n_y), zeros(n_t,n_x,n_y)
px, py, pxy, pz = zeros(n_t,n_x,n_y), zeros(n_t,n_x,n_y), zeros(n_t,n_x,n_y), zeros(n_t,n_x,n_y)
O_phi = zeros(n_t,n_x,n_y)
t = zeros(n_t)
#Matrix of coefficients of fourier modes of the energy density.
if Fourier_modes == "yes"
    a, b = zeros(n_t,order_x+1,order_y+1), zeros(n_t,order_x+1,order_y+1)
end


for it in initial_it:step_it:final_it
    #We first read all the files
    a4 = get_field(files,it=it,field="a4")[1][1,:,:]
    b14 = get_field(files,it=it,field="b14")[1][1,:,:]
    b24 = get_field(files,it=it,field="b24")[1][1,:,:]
    fx2 = get_field(files,it=it,field="fx2")[1][1,:,:]
    fy2 = get_field(files,it=it,field="fy2")[1][1,:,:]
    g4 = get_field(files,it=it,field="g4")[1][1,:,:]
    phi2 = get_field(files,it=it,field="phi2")[1][1,:,:]
    #xi = get_field(files_xi,it=it,field="xi")[1][1,:,:]

   #TODO :Add the last points
    #=append!(a4, a4[1])
    append!(b14, b14[1])
    append!(b24, b24[1])
    append!(fx2, fx2[1])
    append!(fy2, fy2[1])
    append!(g4, g4[1])
    append!(phi2, phi2[1])=#

    i = Int(floor((it-initial_it)/step_it))+1
    #phi2 = phi0^3 .*phi3 .- phi0.*xi.^2

    e[i,:,:] = -3/4 .*a4./phi0^4 - phi2./phi0.^3 .- (9-7*phiM2)/(36*phiM2)
    Jx[i,:,:] = fx2./phi0^4
    Jy[i,:,:] = fy2./phi0^4
    px[i,:,:] = -a4./(4*phi0^4) - b14./phi0^4 - b24./phi0^4 + phi2./(3*phi0^3) .+ (-5/108 + 1/(4*phiM2))
    py[i,:,:] = -a4./(4*phi0^4) + b14./phi0^4 - b24./phi0^4 + phi2./(3*phi0^3) .+ (-5/108 + 1/(4*phiM2))
    pxy[i,:,:] = -g4./phi0^4
    pz[i,:,:] = -a4./(4*phi0^4) + 2 .*b24./phi0^4 + phi2./(3*phi0^3) .+ (-5/108 + 1/(4*phiM2))
    O_phi[i,:,:] = -2 .*phi2./phi0^3 .+ (1/3 - 1/phiM2)

    t[i] = files.current_t

    if Fourier_modes == "yes"
        a[i,:,:], b[i,:,:] = Fourier_modes_2D(e[i,:,:], x, y, order_x, order_y, integration_method)
    end

end

#Writing the data into a h5df file.
if To_hdf5 == "yes"

    output = string(dirname,"Stress-Tensor.h5df")
    fid = h5open(output,"w")

    group_coord = g_create(fid, "Coordinates")
    group_st = g_create(fid, "Stress-Tensor")
    group_svev = g_create(fid, "Scalar_VEV")
    group_fourier = g_create(fid, "Fourier modes")

    Jecco.write_dataset(group_coord, "t",t)
    Jecco.write_dataset(group_coord, "x",x)
    Jecco.write_dataset(group_coord, "y",y)

    Jecco.write_dataset(group_st, "e",e)
    Jecco.write_dataset(group_st, "Jx",Jx)
    Jecco.write_dataset(group_st, "Jy",Jy)
    Jecco.write_dataset(group_st, "px",px)
    Jecco.write_dataset(group_st, "py",py)
    Jecco.write_dataset(group_st, "pxy",pxy)
    Jecco.write_dataset(group_st, "pz",pz)

    Jecco.write_dataset(group_svev, "O_phi",O_phi)

    Jecco.write_dataset(group_fourier, "a",a)
    Jecco.write_dataset(group_fourier, "b",b)

    close(fid)
end

#This writes into a h5df file that has the mathematica array format directly
if To_mathematica == "yes"
    T_m = zeros(11,n_x*n_y*n_t)
    if Fourier_modes == "yes"
        f_m = zeros(2*(order_x+1)*(order_y+1)+1,n_t)
    end
    n = 0
    for i in 1:n_t
        for j in 1:n_x
            for k in 1:n_y
                global n += 1
                T_m[:,n] = [t[i] x[j] y[k] e[i,j,k] Jx[i,j,k] Jy[i,j,k] px[i,j,k] py[i,j,k] pxy[i,j,k] pz[i,j,k] O_phi[i,j,k]]
            end
        end
        if Fourier_modes == "yes"
            ap, bp = reshape(transpose(a[i,:,:]),1,:), reshape(transpose(b[i,:,:]),1,:)
            f_m[:,i] = hcat(t[i], ap, bp)
        end
    end

    output = string(dirname,"data_mathematica.h5df")
    fid = h5open(output,"w")

    group_st = g_create(fid, "data")

    Jecco.write_dataset(group_st, "VEVs",T_m)
    Jecco.write_dataset(group_st, "Fourier modes",f_m)

    close(fid)

end

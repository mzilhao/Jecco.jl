using Jecco
using HDF5
#using Plots
#gr()
include("./Fourier_modes.jl")


dirname = "/Users/apple/Documents/Jecco.jl/data/old_potential_e_10/"
#Specify the initial iteration, the final and the desired step.
initial_it = 0
step_it = 500
final_it = 100000

#Phi0 and inverse of phiM squared
phi0 = 1.0
phiM2 = -(2.3^2)

To_hdf5 = "no"
To_mathematica = "yes"
Fourier_modes = "no"
order_x, order_y = (4,0)
integration_method = "trapezoidal"
#Read_tensor = "yes"

files = OpenPMDTimeSeries(dirname,prefix="boundary")
files_xi = OpenPMDTimeSeries(dirname,prefix="gauge")



u, x, y = get_field(files,it=initial_it,field="a4")[2][:]
nt = Int(floor((final_it-initial_it)/step_it))+1
# Add the last point
append!(x, x[end]+(x[2]-x[1]))
append!(y, y[end]+(y[2]-y[1]))
nx, ny = length(x), length(y)


e, Jx, Jy = zeros(nt,nx,ny), zeros(nt,nx,ny), zeros(nt,nx,ny)
px, py, pxy, pz = zeros(nt,nx,ny), zeros(nt,nx,ny), zeros(nt,nx,ny), zeros(nt,nx,ny)
O_phi = zeros(nt,nx,ny)
t = zeros(nt)
#Matrix of coefficients of fourier modes of the energy density.
if Fourier_modes == "yes"
    a, b = zeros(nt,order_x+1,order_y+1), zeros(nt,order_x+1,order_y+1)
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

    i = Int(floor((it-initial_it)/step_it))+1
    #phi2 = phi0^3 .*phi3 .- phi0.*xi.^2

    e[i,1:end-1,1:end-1] = -3/4 .*a4./phi0^4 - phi2./phi0.^3 .- (9-7*phiM2)/(36*phiM2)
    Jx[i,1:end-1,1:end-1] = fx2./phi0^4
    Jy[i,1:end-1,1:end-1] = fy2./phi0^4
    px[i,1:end-1,1:end-1] = -a4./(4*phi0^4) - b14./phi0^4 - b24./phi0^4 + phi2./(3*phi0^3) .+ (-5/108 + 1/(4*phiM2))
    py[i,1:end-1,1:end-1] = -a4./(4*phi0^4) + b14./phi0^4 - b24./phi0^4 + phi2./(3*phi0^3) .+ (-5/108 + 1/(4*phiM2))
    pxy[i,1:end-1,1:end-1] = -g4./phi0^4
    pz[i,1:end-1,1:end-1] = -a4./(4*phi0^4) + 2 .*b24./phi0^4 + phi2./(3*phi0^3) .+ (-5/108 + 1/(4*phiM2))
    O_phi[i,1:end-1,1:end-1] = -2 .*phi2./phi0^3 .+ (1/3 - 1/phiM2)

    t[i] = files.current_t

    #Add the last points
    e[i,end,:], e[i,:,end], e[i,end,end] = e[i,1,:], e[i,:,1], e[i,1,1]
    Jx[i,end,:], Jx[i,:,end], Jx[i,end,end] = Jx[i,1,:], Jx[i,:,1], Jx[i,1,1]
    Jy[i,end,:], Jy[i,:,end], Jy[i,end,end] = Jy[i,1,:], Jy[i,:,1], Jy[i,1,1]
    px[i,end,:], px[i,:,end], px[i,end,end] = px[i,1,:], px[i,:,1], px[i,1,1]
    py[i,end,:], py[i,:,end], px[i,end,end] = py[i,1,:], py[i,:,1], py[i,1,1]
    pxy[i,end,:], pxy[i,:,end], pxy[i,end,end] = pxy[i,1,:], pxy[i,:,1], pxy[i,1,1]
    pz[i,end,:], pz[i,:,end], pz[i,end,end] = pz[i,1,:], pz[i,:,1], pz[i,1,1]
    O_phi[i,end,:], O_phi[i,:,end], O_phi[i,end,end] = O_phi[i,1,:], O_phi[i,:,1], O_phi[i,1,1]

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
    T_m = zeros(11,nx*ny*nt)
    if Fourier_modes == "yes"
        f_m = zeros(2*(order_x+1)*(order_y+1)+1,nt)
    end
    n = 0
    for i in 1:nt
        for j in 1:nx
            for k in 1:ny
                global n += 1
                T_m[:,n] = [t[i] x[j] y[k] e[i,j,k] Jx[i,j,k] Jy[i,j,k] px[i,j,k] pxy[i,j,k] py[i,j,k] pz[i,j,k] O_phi[i,j,k]]
            end
        end
        if Fourier_modes == "yes"
            f_m[:,i] = hcat(t[i], reshape(transpose(a[i,:,:]),1,:), reshape(transpose(b[i,:,:]),1,:))
        end
    end

    output = string(dirname,"data_mathematica.h5df")
    fid = h5open(output,"w")

    group_st = g_create(fid, "data")

    Jecco.write_dataset(group_st, "VEVs",T_m)
    if Fourier_modes == "yes"
        Jecco.write_dataset(group_st, "Fourier modes",f_m)
    end

    close(fid)

end

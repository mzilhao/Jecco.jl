using Jecco
using HDF5
using Plots
gr()

dirname = "/home/mikel/Documentos/Jecco.jl/BB_AH_01/"

#Phi0 and inverse of phiM squared
phi0 = 0.0
IphiM2 = -1.0

To_hdf5 = "no"
To_mathematica = "yes"
Read_tensor = "yes"

files = OpenPMDTimeSeries(dirname,prefix="boundary")
#Specify the initial iteration, the final and the step in the jump.
initial_it = 0
step_it = 30
final_it = 1830





u, x, y = get_field(files,it=initial_it,field="a4")[2][:]
n_t = Int(floor((final_it-initial_it)/step_it))+1
n_x, n_y = length(x), length(y)

e, Jx, Jy = zeros(n_t,n_x,n_y), zeros(n_t,n_x,n_y), zeros(n_t,n_x,n_y)
px, py, pxy, pz = zeros(n_t,n_x,n_y), zeros(n_t,n_x,n_y), zeros(n_t,n_x,n_y), zeros(n_t,n_x,n_y)
O_phi = zeros(n_t,n_x,n_y)
t = zeros(n_t)

for it in initial_it:step_it:final_it
    #We first read all the files
    a4 = get_field(files,it=it,field="a4")[1][1,:,:]
    b14 = get_field(files,it=it,field="b14")[1][1,:,:]
    b24 = get_field(files,it=it,field="b24")[1][1,:,:]
    fx2 = get_field(files,it=it,field="fx2")[1][1,:,:]
    fy2 = get_field(files,it=it,field="fy2")[1][1,:,:]
    g4 = get_field(files,it=it,field="g4")[1][1,:,:]
    phi3 = get_field(files,it=it,field="phi3")[1][1,:,:]

    i = Int(floor((it-initial_it)/step_it))+1

    e[i,:,:] = -3/4 .*a4 - phi0.*phi3 .- (IphiM2/4 - 7/36)*phi0^4
    Jx[i,:,:] = fx2
    Jy[i,:,:] = fy2
    px[i,:,:] = -a4./4 - b14 - b24 + phi0.*phi3./3 .+ (-5/108 .+ IphiM2/4)*phi0^4
    py[i,:,:] = -a4./4 + b14 - b24 + phi0.*phi3./3 .+ (-5/108 .+ IphiM2/4)*phi0^4
    pxy[i,:,:] = -g4
    pz[i,:,:] = -a4./4 + 2 .*b24 + phi0.*phi3./3 .+ (-5/108 .+ IphiM2/4)*phi0^4
    O_phi[i,:,:] = -2 .*phi3 .+ (1/3 - IphiM2)*phi0^3

    t[i] = files.current_t

end


if To_hdf5 == "yes"

    output = string(dirname,"Stress-Tensor.h5df")
    fid = h5open(output,"w")

    group_coord = g_create(fid, "Coordinates")
    group_st = g_create(fid, "Stress-Tensor")
    group_svev = g_create(fid, "Scalar_VEV")

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

    close(fid)
end

if To_mathematica == "yes"
    file_e = open(string(dirname,"e.m"),"w")
    
    write(file_e,"ee={")
    for i in 1:n_t
        write(file_e,"{")
        for j in 1:n_x
            for k in 1:n_y
                if j == n_x && k == n_y
                    write(file_e, string("{",x[j],",",y[k],",",e[i,j,k],"}"))
                else
                    write(file_e, string("{",x[j],",",y[k],",",e[i,j,k],"},"))
                end
            end
        end
        if i == n_t
            write(file_e,"}")
        else
            write(file_e,"},")
        end
    end
    write(file_e,"}")
    close(file_e)




end





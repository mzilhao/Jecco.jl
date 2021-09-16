using Jecco, Jecco.AdS5_3_1
using HDF5

dirname_16 = "/gpfs/projects/ub48/ub48946/jecco_test/convergence_test_16/"
dirname_32 = "/gpfs/projects/ub48/ub48946/jecco_test/convergence_test_32/"
dirname_64 = "/gpfs/projects/ub48/ub48946/jecco_test/convergence_test_64/"
output     = "/gpfs/projects/ub48/ub48946/jecco_test/convergence_test_result/"


A_16_c1 = ConstrainedTimeSeries(dirname_16, :A, 1)
A_32_c1 = ConstrainedTimeSeries(dirname_32, :A, 1)
A_64_c1 = ConstrainedTimeSeries(dirname_64, :A, 1)

A_16_c2 = ConstrainedTimeSeries(dirname_16, :A, 2)
A_32_c2 = ConstrainedTimeSeries(dirname_32, :A, 2)
A_64_c2 = ConstrainedTimeSeries(dirname_64, :A, 2)

t_16, u_16_c1, x_16, y_16 = get_coords(A_16,:,:,:,:)
t_32, u_32_c1, x_32, y_32 = get_coords(A_32,:,:,:,:)
t_64, u_64_c1, x_64, y_64 = get_coords(A_64,:,:,:,:)

t_16, u_16_c2, x_16, y_16 = get_coords(A_16,:,:,:,:)
t_32, u_32_c2, x_32, y_32 = get_coords(A_32,:,:,:,:)
t_64, u_64_c2, x_64, y_64 = get_coords(A_64,:,:,:,:)

Nu = length(u_16_c2)

h_16 = x_16[2]-x_16[1]
h_32 = x_32[2]-x_32[1]
h_64 = x_64[2]-x_64[1]
n    = 2
Q    = (h_16^n-h_32^n)/(h_32^n-h_64^n)

∆A_16_32_c1 = similar(A_16_c1)
∆A_32_64_c1 = similar(A_32_c1)
∆A_16_32_c2 = similar(A_16_c2)
∆A_32_64_c2 = similar(A_32_c2)

max_index_t = length(t_64)

∆A_16_32_c1 .= abs.(A_16_c1[1:max_index_t,:,:,:] - A_32_c1[1:max_index_t,:,1:2:32,1:2:32])
∆A_32_64_c1 .= abs.(A_32_c1[1:max_index_t,:,1:2:32,1:2:32] - A_64_c1[1:max_index_t,:,1:4:64,1:4:64])
∆A_16_32_c2 .= abs.(A_16_c2[1:max_index_t,:,:,:] - A_32_c2[1:max_index_t,:,1:2:32,1:2:32])
∆A_32_64_c2 .= abs.(A_32_c2[1:max_index_t,:,1:2:32,1:2:32] - A_64_c2[1:max_index_t,:,1:4:64,1:4:64])

Q_c1 = ∆A_1

norm_convergence = zeros(max_index_t)


for i in 1:max_index_t
    s_16_32_c1 = sum(∆A_16_32_c1[i,:,:,:].^2)
    s_32_64_c1 = sum(∆A_32_64_c1[i,:,:,:].^2)

    s_16_32_c2 = sum(∆A_16_32_c2[i,2:Nu,:,:].^2)
    s_32_64_c2 = sum(∆A_32_64_c2[i,2:Nu,:,:].^2)

    norm_convergence[i] = sqrt((s_16_32_c1+s_16_32_c2)/(s_32_64_c1+s_32_64_c2))
end

fid = h5open(output,"w")
group_st = g_create(fid, "data")
Jecco.write_dataset(group_st, "Q", Q)
Jecco.write_dataset(group_st, "16_32_c1", ∆A_16_32_c1)
Jecco.write_dataset(group_st, "16_32_c2", ∆A_16_32_c2)
Jecco.write_dataset(group_st, "32_64_c1", ∆A_32_64_c1)
Jecco.write_dataset(group_st, "32_64_c2", ∆A_32_64_c2)
Jecco.write_dataset(group_st, "Norm", norm_convergence)

using HDF5
using Plots



file = "/home/mikel/Documents/Jecco.jl/data/convergence_test/convergence_test.h5"

t        = h5read(file, "data/t")
x        = h5read(file, "data/x")
y        = h5read(file, "data/y")
dA_16_32_c1 = h5read(file, "data/16_32_c1")
dA_32_64_c1 = h5read(file, "data/32_64_c1")
norm     = h5read(file, "data/Norm")

#plot(t, norm, lw=3)
tt=1200
yy=6
uu=6
plot(x, dA_16_32_c1[tt,uu,:,yy]./dA_32_64_c1[tt,uu,:,yy], lw=3)
#plot!(x, 8 .*dA_32_64_c1[tt,uu,:,yy], lw=3)

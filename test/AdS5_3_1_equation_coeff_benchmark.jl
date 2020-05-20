using BenchmarkTools

using Jecco
using Jecco.AdS5_3_1

ABCS  = zeros(4)
AA    = zeros(2,2)
BB    = zeros(2,2)
CC    = zeros(2,2)
SS    = zeros(2)

svars    = ones(11)
fvars    = ones(37)
sdvars   = tuple([ZeroPotential(); ones(70)]...)
bdgvars  = ones(71)
phidvars = tuple([ZeroPotential(); ones(72)]...)
avars    = tuple([ZeroPotential(); ones(76)]...)


# inner grid

println("S equation inner grid coefficients:")
@btime Jecco.AdS5_3_1.S_eq_coeff!($ABCS, $svars, Inner())

println("Fxy equations inner grid coefficients:")
@btime Jecco.AdS5_3_1.Fxy_eq_coeff!($AA, $BB, $CC, $SS, $fvars, Inner())

println("Sd equation inner grid coefficients:")
@btime Jecco.AdS5_3_1.Sd_eq_coeff!($ABCS, $sdvars, Inner())

println("B2d equation inner grid coefficients:")
@btime Jecco.AdS5_3_1.B2d_eq_coeff!($ABCS, $bdgvars, Inner())

println("B1dGd equations inner grid coefficients:")
@btime Jecco.AdS5_3_1.B1dGd_eq_coeff!($AA, $BB, $CC, $SS, $bdgvars, Inner())

println("phid equations inner grid coefficients:")
@btime Jecco.AdS5_3_1.phid_eq_coeff!($ABCS, $phidvars, Inner())

println("A equations inner grid coefficients:")
@btime Jecco.AdS5_3_1.A_eq_coeff!($ABCS, $avars, Inner())


# outer grid

println("S equation outer grid coefficients:")
@btime Jecco.AdS5_3_1.S_eq_coeff!($ABCS, $svars, Outer())

println("Fxy equations outer grid coefficients:")
@btime Jecco.AdS5_3_1.Fxy_eq_coeff!($AA, $BB, $CC, $SS, $fvars, Outer())

println("Sd equation outer grid coefficients:")
@btime Jecco.AdS5_3_1.Sd_eq_coeff!($ABCS, $sdvars, Outer())

println("B2d equation outer grid coefficients:")
@btime Jecco.AdS5_3_1.B2d_eq_coeff!($ABCS, $bdgvars, Outer())

println("B1dGd equations outer grid coefficients:")
@btime Jecco.AdS5_3_1.B1dGd_eq_coeff!($AA, $BB, $CC, $SS, $bdgvars, Outer())

println("phid equations outer grid coefficients:")
@btime Jecco.AdS5_3_1.phid_eq_coeff!($ABCS, $phidvars, Outer())

println("A equations outer grid coefficients:")
@btime Jecco.AdS5_3_1.A_eq_coeff!($ABCS, $avars, Outer())

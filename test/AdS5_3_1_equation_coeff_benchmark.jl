using BenchmarkTools

using Jecco
using Jecco.AdS5_3_1

ABCS  = zeros(4)
AA    = zeros(2,2)
BB    = zeros(2,2)
CC    = zeros(2,2)
SS    = zeros(2)
abccS = 0

svars    = tuple(ones(11)...)
fvars    = tuple(ones(37)...)
sdvars   = tuple([ZeroPotential(); ones(70)]...)
bdgvars  = tuple(ones(71)...)
phidvars = tuple([ZeroPotential(); ones(72)]...)
avars    = tuple([ZeroPotential(); ones(76)]...)
xivars   = tuple(ones(84)...)

# inner grid

println("S equation inner grid coefficients:")
@btime Jecco.AdS5_3_1.S_eq_coeff!($ABCS, $svars, AdS5_3_1.Inner())

println("Fxy equations inner grid coefficients:")
@btime Jecco.AdS5_3_1.Fxy_eq_coeff!($AA, $BB, $CC, $SS, $fvars, AdS5_3_1.Inner())

println("Sd equation inner grid coefficients:")
@btime Jecco.AdS5_3_1.Sd_eq_coeff!($ABCS, $sdvars, AdS5_3_1.Inner())

println("B2d equation inner grid coefficients:")
@btime Jecco.AdS5_3_1.B2d_eq_coeff!($ABCS, $bdgvars, AdS5_3_1.Inner())

println("B1dGd equations inner grid coefficients:")
@btime Jecco.AdS5_3_1.B1dGd_eq_coeff!($AA, $BB, $CC, $SS, $bdgvars, AdS5_3_1.Inner())

println("phid equations inner grid coefficients:")
@btime Jecco.AdS5_3_1.phid_eq_coeff!($ABCS, $phidvars, AdS5_3_1.Inner())

println("A equations inner grid coefficients:")
@btime Jecco.AdS5_3_1.A_eq_coeff!($ABCS, $avars, AdS5_3_1.Inner())


# outer grid

println("S equation outer grid coefficients:")
@btime Jecco.AdS5_3_1.S_eq_coeff!($ABCS, $svars, AdS5_3_1.Outer())

println("Fxy equations outer grid coefficients:")
@btime Jecco.AdS5_3_1.Fxy_eq_coeff!($AA, $BB, $CC, $SS, $fvars, AdS5_3_1.Outer())

println("Sd equation outer grid coefficients:")
@btime Jecco.AdS5_3_1.Sd_eq_coeff!($ABCS, $sdvars, AdS5_3_1.Outer())

println("B2d equation outer grid coefficients:")
@btime Jecco.AdS5_3_1.B2d_eq_coeff!($ABCS, $bdgvars, AdS5_3_1.Outer())

println("B1dGd equations outer grid coefficients:")
@btime Jecco.AdS5_3_1.B1dGd_eq_coeff!($AA, $BB, $CC, $SS, $bdgvars, AdS5_3_1.Outer())

println("phid equations outer grid coefficients:")
@btime Jecco.AdS5_3_1.phid_eq_coeff!($ABCS, $phidvars, AdS5_3_1.Outer())

println("A equations outer grid coefficients:")
@btime Jecco.AdS5_3_1.A_eq_coeff!($ABCS, $avars, AdS5_3_1.Outer())

#xi equation

println("xi_t equations coefficients:")
@btime $abccS = Jecco.AdS5_3_1.xi_t_eq_coeff($xivars, AdS5_3_1.Outer())


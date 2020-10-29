using Test
using Jecco
using Jecco.AdS5_3_1


sigmavars   = tuple(ones(58)...)

axx, ayy, axy, bx, by, cc= AdS5_3_1.AH_eq_coeff(sigmavars, AdS5_3_1.Outer())

S = AdS5_3_1.AH_eq_res(sigmavars, AdS5_3_1.Outer())

print("(axx, ayy, axy, bx, by, cc, S) =",(axx, ayy, axy, bx, by, cc, S),"\n")


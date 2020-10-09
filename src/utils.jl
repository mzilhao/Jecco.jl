
#= iterative improvement to a solution of a linear system A x = b

cf. Numerical Recipes (2007), sec 2.5

this is the same procedure as gsl_linalg_LU_refine:
https://github.com/ampl/gsl/blob/master/linalg/lu.c

the algorithm is:

1. compute the residual = A * x - b

2. solve for the correction delta, A * delta = residual

3. subtract this from the old solution for the improved solution
   xnew = x - delta

Note: in Numerical Recipes it is suggested to use higher precision to compute
the residual, since this can be very small. from my (very limited) testing, this
didn't seem crucial, and i'm also not sure what would be the best way to do that
here. the issue is that we just want to compute the residual in higher
precision, not performing the left division. so, if we were to expect the caller
to convert to higher precision upon calling this function, we'd still have to
convert back to lower precision before calling ldiv! (for efficiency reasons and
also since most ldiv methods don't work for higher precision types anyway). so
for now it stays like this.
=#
function refine_solution!(sol::Vector, A_fact, A_mat::AbstractMatrix, b_vec::Vector)
    mul!(b_vec, A_mat, sol, 1, -1)
    ldiv!(A_fact, b_vec)
    sol .-= b_vec
end

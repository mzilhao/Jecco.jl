# Convergence tests

*Instructions to reproduce the convergence tests of the section 4.3 of
arxiv xxxx.xxxxx.*

If you want to obtain the results of sectio 4.3 you have to first run
the three julia scripts in the directory, these are t10_Nx16.jl,
t10_N32.jl and t10_N64.jl. Here, 16, 32 and 64 refer to the number of
nodes on each trasverse domain and 10 is the total time in code
units.

To run the example of choice move to the directory where the example
is located and type in bash
``` julia example.jl ```

Multi-threading is supported. To enable it type in bash
```
export JULIA_NUM_THREADS=n
julia example.jl
```
where n is the number of threads in multi-threading.

Data are saved in "examples/convergence_tests" and their analysis can
be performed with the jupyter notebook
"examples/convergence_tests/convergence_plots.ipynb".


function run_model(par_base, par_grid, par_id, par_evol, par_io)
    Jecco.startup()

    # define potential
    Jecco.AdS_3_1.setup(par_base)

    Jecco.AdS_3_1.ibvp(par_grid, par_id, par_evol, par_io)

    println("-------------------------------------------------------------")
    println("Done.")
end

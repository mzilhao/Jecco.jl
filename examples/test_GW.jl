using Jecco, Jecco.AdS5_3_1


dirname = "/home/mikel/Dropbox/Spinodal_GWs/mathematica/data/spinodal_L_63.6_N_100_same_a4/"
# eos     = "/home/mikel/Dropbox/Spinodal_GWs/mathematica/data/eos.h5"
outdir  = dirname
dt = 0.05
δk = 0.1
AdS5_3_1.solve_GW(outdir, dirname, dt=dt, alg=AdS5_3_1.AB3())
# AdS5_3_1.get_TT_stress_tensor(dirname, dirname)
convert_to_mathematica(dirname)
AdS5_3_1.modes_to_mathematica(dirname, :energy, "E_modes.h5")
# AdS5_3_1.modes_to_mathematica(dirname, :px, "Txx_modes.h5")
# AdS5_3_1.modes_to_mathematica(dirname, :pxy, "Txy_modes.h5")
# AdS5_3_1.modes_to_mathematica(dirname, :py, "Tyy_modes.h5")
# AdS5_3_1.modes_to_mathematica(dirname, :pz, "Tzz_modes.h5")
# AdS5_3_1.modes_to_mathematica(dirname, :Ideal_Txx, "Ideal_Txx_modes.h5", path_to_eos=eos)
# AdS5_3_1.modes_to_mathematica(dirname, :Ideal_Txy, "Ideal_Txy_modes.h5", path_to_eos=eos)
# AdS5_3_1.modes_to_mathematica(dirname, :Ideal_Tyy, "Ideal_Tyy_modes.h5", path_to_eos=eos)
# AdS5_3_1.modes_to_mathematica(dirname, :Ideal_Tzz, "Ideal_Tzz_modes.h5", path_to_eos=eos)
# AdS5_3_1.modes_to_mathematica(dirname, :vivj, "vivj_modes.h5")
# AdS5_3_1.modes_to_mathematica(dirname, :Txx, "TTxx_modes.h5")
# AdS5_3_1.modes_to_mathematica(dirname, :hdxx, "hdxx_modes.h5")
# AdS5_3_1.modes_to_mathematica(dirname, :hdxy, "hdxy_modes.h5")
# AdS5_3_1.modes_to_mathematica(dirname, :hdyy, "hdyy_modes.h5")
# AdS5_3_1.modes_to_mathematica(dirname, :hdzz, "hdzz_modes.h5")
AdS5_3_1.modes_to_mathematica(dirname, :drhod2k, "drhod2k_modes.h5")
AdS5_3_1.drhodlogk_to_mathematica_1(dirname, δk, outfile="drhoGWdlogk_mathematica.h5")
# AdS5_3_1.ρGW(dirname)
nothing

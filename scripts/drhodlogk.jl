using Jecco, Jecco.AdS5_3_1

dir    = "/home/mikel/Dropbox/Spinodal_GWs/mathematica/data/spinodal_4/"
δk = 0.1
@time AdS5_3_1.drhodlogk_to_mathematica_1(dir, δk, outfile="drhoGWdlogk_mathematica_"*string(δk)*".h5")

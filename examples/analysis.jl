using Vivi
using Plots
gr()

ts = OpenPMDTimeSeries("./data")

phi, info=get_field(ts, it=2, field="phi")

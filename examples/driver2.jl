

# This is a Julia port of SunlightHNC 
# Copyright (C) 2025 Hiroharu Sugawara, Tokyo Metropolitan University <hsugawa@tmu.ac.jp>

# Based on an original code copyright (c) 2007 Lucian Anton.
# Modifications copyright (c) 2008, 2009 Andrey Vlasov.  Additional
# modifications copyright (c) 2009-2017 Unilever UK Central Resources
# Ltd (Registered in England & Wales, Company No 29140; Registered
# Office: Unilever House, Blackfriars, London, EC4P 4BQ, UK).

using SunlightHNC
using Printf

ncomp = 2
ng = 2^12
w2 = Wizard(; ncomp=ncomp, ng=ng)
w2.sigma = 0.5e0
w2.lb = 50.0e0
w2.arep .= 25.0e0
w2.arep[1, 2] = 30.0e0
w2.z[1] = 1.0e0
w2.z[2] = -1.0e0
dpd_potential!(w2)

rhotot = 3.0e0
w2.rho .= 0.5e0 * rhotot
write_params(w2)
w2.verbose = true
w2.maxsteps = 1000
hnc_solve!(w2)
@show w2.ERROR

ddsf = calc_ddsf(w2)
ccsf = calc_ccsf(w2)

for j in 1:20 #  w2.ng-1
    @printf "%5d %15.7e %15.7e %15.7e\n" j w2.k[j] ddsf[j] ccsf[j]
end

using Plots
plot()
plot!(Series_vs_k(), w2, ddsf; label="ddsf")
plot!(Series_vs_k(), w2, ccsf; label="ccsf")
plot!(title="Structure factors")
plot!(xlims=(0,20))
savefig("driver2_snn.png")

plot(Series_hr(), w2; rcut=10)
plot!(xlims=(0,4))
savefig("driver2_hr.png")

write_thermodynamics(w2)

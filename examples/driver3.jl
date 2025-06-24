# This is a Julia port of SunlightHNC 
# Copyright (C) 2025 Hiroharu Sugawara, Tokyo Metropolitan University <hsugawa@tmu.ac.jp>

# Based on an original code copyright (c) 2007 Lucian Anton.
# Modifications copyright (c) 2008, 2009 Andrey Vlasov.  Additional
# modifications copyright (c) 2009-2017 Unilever UK Central Resources
# Ltd (Registered in England & Wales, Company No 29140; Registered
# Office: Unilever House, Blackfriars, London, EC4P 4BQ, UK).
# Additional modifications copyright (c) 2020-2021 Patrick B Warren
# <patrick.warren@stfc.ac.uk> and STFC.

using SunlightHNC
using Printf

ng = 4096
ncomp = 3
ncomp = 3
w3 = Wizard(; ng=ng, ncomp=ncomp)
w3.sigma = 0.5e0
w3.lb = 20.0e0
w3.arep .= 25.0e0
w3.arep[1, 2] = 30.0e0
w3.arep[1, 3] = 27.0e0
w3.arep[2, 3] = 20.0e0
w3.z[1] = 1.0e0
w3.z[2] = -1.0e0

dpd_potential!(w3)

rhotot = 3.0e0
mfcharge = 0.2e0

w3.rho[1] = 0.5e0 * rhotot * mfcharge
w3.rho[2] = w3.rho[1]
w3.rho[3] = rhotot * (1.0e0 - mfcharge)

dpd_potential!(w3) # , charge_type=OZ.DPD_GAUSSIAN_CHARGES)
write_params(w3)
w3.verbose = true
w3.maxsteps = 1000
hnc_solve!(w3)
ddsf = calc_ddsf(w3)
ccsf = calc_ccsf(w3)

for j in 1:20 #  w2.ng-1
    @printf "%5d %15.7e %15.7e %15.7e\n" j w3.k[j] ddsf[j] ccsf[j]
end

using Plots
plot()
plot!(Series_vs_k(), w3, ddsf; label="ddsf")
plot!(Series_vs_k(), w3, ccsf; label="ccsf")
plot!(title="Structure factors")
plot!(xlims=(0,20))
savefig("driver2_snn.png")

plot(Series_hr(), w3; rcut=10)
plot!(xlims=(0,4))
savefig("driver3_hr.png")

write_thermodynamics(w3)
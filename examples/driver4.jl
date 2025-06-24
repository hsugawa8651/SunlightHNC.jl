# This is a Julia port of SunlightHNC 
# Copyright (C) 2025 Hiroharu Sugawara, Tokyo Metropolitan University <hsugawa@tmu.ac.jp>

# Based on an original code copyright (c) 2007 Lucian Anton.
# Modifications copyright (c) 2008, 2009 Andrey Vlasov.  Additional
# modifications copyright (c) 2009-2017 Unilever UK Central Resources
# Ltd (Registered in England & Wales, Company No 29140; Registered
# Office: Unilever House, Blackfriars, London, EC4P 4BQ, UK).

using SunlightHNC
using Printf

ng = 65536
ng = 4096
ncomp = 4
w4 = Wizard(; ng=ng, ncomp=ncomp)
w4.sigma = 0.5e0
w4.lb = 25.0e0
w4.arep .= 25.0e0
w4.arep[1, 2] = 28.0e0
w4.arep[1, 3] = 22.0e0
w4.arep[1, 4] = 15.0e0
w4.arep[2, 3] = 30.0e0
w4.arep[2, 4] = 26.0e0
w4.arep[3, 4] = 20.0e0

w4.z[2] = 1.0e0
w4.z[3] = -1.0e0
w4.z[4] = -1.0e0

dpd_potential!(w4)

rhotot = 3.0e0
mfcharge = 0.2e0

w4.rho[1] = rhotot * (1.0e0 - mfcharge)
w4.rho[2] = 0.5e0 * rhotot * mfcharge
w4.rho[3] = w4.rho[2] * 0.2e0
w4.rho[4] = w4.rho[2] - w4.rho[3]

write_params(w4)
# w4.verbose = true
# w4.maxsteps = 1000
hnc_solve!(w4)
write_thermodynamics(w4)

ddsf = calc_ddsf(w4)
ccsf = calc_ccsf(w4)

for j in 1:20 #  w2.ng-1
  @printf "%5d %15.7e %15.7e %15.7e\n" j w4.k[j] ddsf[j] ccsf[j]
end

using Plots
plot()
plot!(Series_vs_k(), w4, ddsf; label="ddsf")
plot!(Series_vs_k(), w4, ccsf; label="ccsf")
plot!(title="Structure factors")
plot!(xlims=(0,20))
savefig("driver4_snn.png")

plot(Series_hr(), w4; rcut=10)
plot!(xlims=(0,4))
savefig("driver4_hr.png")

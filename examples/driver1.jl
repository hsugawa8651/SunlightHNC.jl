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

using Plots

ng = 2^14
w1 = Wizard(;ng=ng)
w1.arep .= 25.0e0
dpd_potential!(w1)
w1.rho .= 3.0e0
write_params(w1)

w1.verbose = true
hnc_solve!(w1)
@show w1.ERROR

for j in 1:200
   @printf "%12.5f %12.5f\n" w1.k[j] w1.sk[j,1,1]
end 

using Plots
plot(Series_sk(), w1; kcut=30)
plot!(xlims=(0,20))
savefig("driver1_sk.png")

plot(Series_hr(), w1; rcut=10)
plot!(xlims=(0,4))
savefig("driver1_hr.png")

write_thermodynamics(w1);


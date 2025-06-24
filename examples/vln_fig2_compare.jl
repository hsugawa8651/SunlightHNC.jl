# This is a Julia port of SunlightHNC 
# Copyright (C) 2025 Hiroharu Sugawara, Tokyo Metropolitan University <hsugawa@tmu.ac.jp>

# Based on an original code copyright (c) 2007 Lucian Anton.
# Modifications copyright (c) 2008, 2009 Andrey Vlasov.  
# Additional modifications copyright (c) 2009-2017 Unilever UK Central
# Resources Ltd (Registered in England & Wales, Company No 29140;
# Registered Office: Unilever House, Blackfriars, London, EC4P 4BQ,
# UK).

using SunlightHNC
using Printf

ncomp = 2
w = Wizard(;ncomp=ncomp)

w.rho[1] = 3.0
w.rho[2] = 0.0 # zero density (still works!)

Amin = 25.0
Amax = 106.5
npt = 50
As = LinRange(Amin, Amax, npt)

for i in 1:npt
    A = As[i]
    w.arep[1,1] = w.arep[1,2] = w.arep[2,2] = A
    dpd_potential!(w)
    hnc_solve!(w)
    @show i, A, w.muex[1], w.muex[2], w.ERROR
end

A = As[npt]

dAmin = -5.0
dAmax = 20.0
npt = 11
dAs = LinRange(dAmin, dAmax, npt)
ys = zeros(npt)
for i in 1:npt
    dA=dAs[i]
    w.arep[1,2] = A + dA
    dpd_potential!(w)
    hnc_solve!(w)
    ys[i] = (w.muex[2] - w.muex[1]) / log(10.0) # ie log_10(gamma)
    @show i, A, w.muex[1], w.muex[2], w.ERROR
end

xx = [dAmin, dAmax] # first and last elements
yy = [0.144*x for x in xx] # straight line slope 0.144

using Plots
using LaTeXStrings
plot()
plot!(title="Fig.2, Vishnyakov, Lee and Neimark (2013)")
plot!(xlabel=L"$\Delta{a}$", ylabel=L"$\log_{10}{\gamma}$")
plot!(dAs, ys, label="HNC")
plot!(xx, yy, label=L"$0.144\Delta{a}$")
savefig("vln_fig2.png")

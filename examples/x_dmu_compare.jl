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


ncomp = 2
w = Wizard(; ncomp=ncomp)
rho = 3.0
A = 25.0
dA = 5.0

w.arep[1, 1] = w.arep[2, 2] = A
w.arep[1, 2] = A + dA
dpd_potential!(w)
npt = 41
xs = LinRange(0, 1, npt)
ys = zeros(npt)
for i in 1:npt
    x = xs[i]
    w.rho[1] = rho * (1 - x)
    w.rho[2] = rho * x
    hnc_solve!(w)
    ys[i] = w.muex[2] - w.muex[1]
end


#  x   d(mu)  std-error
data_text = """ 
 0.0,  1.702,  0.004
 0.1,  1.281,  0.003
 0.2,  0.915,  0.003
 0.3,  0.597,  0.004
 0.4,  0.282,  0.004
 0.5, -0.002,  0.003
 0.6, -0.290,  0.003
 0.7, -0.596,  0.003
 0.8, -0.920,  0.004
 0.9, -1.284,  0.003
 1.0, -1.702,  0.004"""

using DelimitedFiles
df = readdlm(IOBuffer(data_text), ',', skipstart=1)
@show size(df)

using Plots
using LaTeXStrings
plot()
plot!(xlabel=L"$x$", ylabel=L"\Delta\mu")
plot!(xs, ys, label="HNC")
scatter!(df[:, 1], df[:, 2],
    label="Monte-Carlo",
    markersize=8, markerstrokewidth=0, dpi=150)
savefig("x_dmu-plot.png")

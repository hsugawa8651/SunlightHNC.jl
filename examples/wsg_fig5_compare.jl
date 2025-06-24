# This is a Julia port of SunlightHNC 
# Copyright (C) 2025 Hiroharu Sugawara, Tokyo Metropolitan University <hsugawa@tmu.ac.jp>

# Based on an original code copyright (c) 2007 Lucian Anton.
# Modifications copyright (c) 2008, 2009 Andrey Vlasov.  Additional
# modifications copyright (c) 2009-2017 Unilever UK Central Resources
# Ltd (Registered in England & Wales, Company No 29140; Registered
# Office: Unilever House, Blackfriars, London, EC4P 4BQ, UK).
# Additional modifications copyright (c) 2020-2021 Patrick B Warren
# <patrick.warren@stfc.ac.uk> and STFC.

using DelimitedFiles

# dA dF
data1_text = """
0, 0
2, 0.699745
4, 1.370293
6, 2.0262583
8, 2.6384122
10, 3.2359834
12, 3.7897591
14, 4.3435348
16, 4.853515
18, 5.3488966
20, 5.8296797
22, 6.2812658
24, 6.7182534
26, 7.1406424
28, 7.5338344
30, 7.9270106
"""
df1 = readdlm(IOBuffer(data1_text), ',', skipstart=1)
# @show size(df1)
wsgx = view(df1[:, 1], :)
wsgy = view(df1[:, 2], :)

# dA   d(mu)   std-err
data2_text = """
0, 0.0007, 0.00299772
1, 0.3571, 0.00355517
2, 0.7052, 0.004213
3, 1.0451, 0.00498353
4, 1.3771, 0.00587967
5, 1.7014, 0.00690962
6, 2.0181, 0.00808128
7, 2.3276, 0.00940074
8, 2.6299, 0.0108687
9, 2.9254, 0.0124841
10, 3.2143, 0.014242
11, 3.4969, 0.0161377
12, 3.7734, 0.0181624
13, 4.0441, 0.0203056
14, 4.3093, 0.0225574
15, 4.5692, 0.0249076
16, 4.824, 0.0273447
17, 5.0741, 0.0298606
18, 5.3197, 0.0324429
19, 5.5609, 0.0350847
20, 5.7982, 0.0377775
21, 6.0317, 0.0405129
22, 6.2615, 0.0432836
23, 6.488, 0.0460852
24, 6.7113, 0.0489117
25, 6.9317, 0.0517575
26, 7.1492, 0.0546188
27, 7.3641, 0.0574919
28, 7.5765, 0.0603722
29, 7.7866, 0.0632587
30, 7.9945, 0.0661484
"""

df2 = readdlm(IOBuffer(data2_text), ',', skipstart=1)
# @show size(df2)
mcx = view(df2[:, 1], :)
mcy = view(df2[:, 2], :)
mce = view(df2[:, 3], :)

using Plots
using LaTeXStrings
plot()
plot!(xlabel=L"$\Delta{A}$", ylabel=L"$\Delta{F}$")
scatter!(mcx, mcy,
    label="Widom insertion",
    markersize=8, markerstrokewidth=0, dpi=150)
scatter!(wsgx, wsgy,
    label="Wijmans et al. (2001)",
    markersize=8, markerstrokewidth=0, dpi=150)

using SunlightHNC
using Printf

ncomp = 2
w = Wizard(; ncomp=ncomp)

rho = 3.0
Amin = 25.0
A = Amin
Amax = 55.0

w.arep[1, 1] = A
w.arep[2, 2] = A
w.rho[1] = rho
w.rho[2] = 0.0
npt = 41
As = range(Amin, Amax, length=npt)
ys = zeros(Float64, npt)
for i in 1:npt
    global A = As[i]
    dA = A - Amin
    w.arep[1, 2] = A
    dpd_potential!(w)
    hnc_solve!(w)
    ys[i] = w.muex[2] - w.muex[1]
end

plot!(As .- Amin, ys, label="HNC")
savefig("wsg_fig5-plot.png")

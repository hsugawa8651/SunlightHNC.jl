# This is a Julia port of SunlightHNC 
# Copyright (C) 2025 Hiroharu Sugawara, Tokyo Metropolitan University <hsugawa@tmu.ac.jp>

# Based on an original code copyright (c) 2007 Lucian Anton.
# Modifications copyright (c) 2008, 2009 Andrey Vlasov.  
# Additional modifications copyright (c) 2009-2017 Unilever UK Central
# Resources Ltd (Registered in England & Wales, Company No 29140;
# Registered Office: Unilever House, Blackfriars, London, EC4P 4BQ,
# UK).


using SunlightHNC

ncomp = 1
ng = 65536
ng = 4096
deltar = 0.001
w = Wizard(; ncomp=ncomp, ng=ng, deltar=deltar)
w.sigma = 1.0
eta = 0.49e0
w.rho[1] = eta * 6.0 / (pi * w.sigma^3)
hs_potential!(w)
write_params(w)
# npt = 150
hnc_solve!(w)

# w.verbose = true
# w.maxsteps = 1000

using Plots
using LaTeXStrings
plot()
plot!(Series_gr(), w;
    title="Fig.4.3, Hansen and McDonald, 4th ed. (2013)")
scatter!([], []; markershape=:x, markersize=0, markerstrokewidth=0, color=:white,
    label="Hard-sphere liquid")
scatter!([], []; markershape=:x, markersize=0, markerstrokewidth=0, color=:white,
    label=L"$\eta=%$(eta)$")
plot!(xlims=(0.5, 3), ylims=(-0.1, 6))

savefig("hm4_fig4-3_plot.png")


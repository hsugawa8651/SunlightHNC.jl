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
hs_potential!(w)
write_params(w)
npt = 150
rho_hi = 0.43 * 6.0 / (pi * w.sigma^3)
drho = rho_hi / npt
rho_s = ((1:npt) .+ 1) * drho
eta_s = rho_s * w.sigma^3 * pi / 6.0
carnahan_starling = (1 .+ eta_s .+ eta_s .^ 2 .- eta_s .^ 3) ./ (1 .- eta_s) .^ 3
percus_yevick_virial = (1 .+ 2 * eta_s .+ 3 * eta_s .^ 2) ./ (1 .- eta_s) .^ 2
percus_yevick_comp = (1 .+ eta_s .+ eta_s .^ 2) ./ (1 .- eta_s) .^ 3

# w.verbose = true
# w.maxsteps = 1000

hnc_vs = zeros(npt)
hnc_cs = zeros(npt)
msa_vs = zeros(npt)
msa_cs = zeros(npt)
hnc_p_xc = hnc_prev = 0.0e0
msa_p_xc = msa_prev = 0.0e0

for i in 1:npt
    rho = rho_s[i]
    @show i, rho
    w.rho[1] = rho
    hnc_solve!(w)
    hnc_vs[i] = w.press / rho
    global hnc_p_xc += 0.5 * drho * (hnc_prev + w.comp_xc)
    global hnc_prev = w.comp_xc
    hnc_cs[i] = hnc_p_xc / rho + 1
    #
    msa_solve!(w)
    msa_vs[i] = w.press / rho
    global msa_p_xc += 0.5 * drho * (msa_prev + w.comp_xc)
    global msa_prev = w.comp_xc[]
    msa_cs[i] = msa_p_xc / rho + 1
end

using Plots
using LaTeXStrings
plot()
plot!(xlabel=L"\eta", ylabel=L"$\beta{P} /\rho$")
plot!(eta_s, carnahan_starling, label="Carnahan-Starling")
plot!(eta_s, percus_yevick_virial, label="PY(v) exact")
plot!(eta_s, percus_yevick_comp, label="PY(c) exact")
plot!(eta_s, hnc_vs, label="HNC(v)")
plot!(eta_s, hnc_cs, label="HNC(c)")
plot!(eta_s, msa_vs, label="MSA(v)")
plot!(eta_s, msa_cs, label="MSA(c)")
vline!([1], label="")
plot!(title="Fig.4.2, Hansen and McDonald, 4th ed. (2013)")
plot!(legend=:topleft)
plot!(xlims=(0.05, 0.5), ylims=(0, 10))

savefig("hm4_fig4-2_plot.png")


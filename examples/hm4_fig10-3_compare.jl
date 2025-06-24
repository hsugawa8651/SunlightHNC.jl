# This is a Julia port of SunlightHNC 
# Copyright (C) 2025 Hiroharu Sugawara, Tokyo Metropolitan University <hsugawa@tmu.ac.jp>

# Copyright (c) 2009-2019 Unilever UK Central Resources Ltd
# (Registered in England & Wales, Company No 29140; Registered Office:
# Unilever House, Blackfriars, London, EC4P 4BQ, UK).  Additional
# modifications copyright (c) 2020-2024 Patrick B Warren
# <patrick.warren@stfc.ac.uk> and STFC.



# c_s/M      -U/NkT (MC)      contact (MC)     p/(rho kT)

data_text = """
0.00911,  0.1029, 0.0013,  0.0044, 0.0007,  0.9701, 0.0008
0.10376,  0.2739, 0.0014,  0.0359, 0.0011,  0.9445, 0.0012
0.42502,  0.4341, 0.0017,  0.1217, 0.0045,  0.9774, 0.0046
1.0001,   0.5516, 0.0016,  0.2777, 0.0045,  1.094,  0.005
1.9676,   0.6511, 0.0020,  0.5625, 0.0088,  1.346,  0.009
"""

using DelimitedFiles
df = readdlm(IOBuffer(data_text), ',', skipstart=1)
# @show size(df)

d = 0.425
lb = 0.71

ndata = size(df, 2)
sqrdat = sqrt.(d^3 * 1.204 * df[:, 1])
nundat = view(df[:, 2], :)
ctcdat = view(df[:, 4], :)
cmpdat = view(df[:, 6], :)

using Plots
using LaTeXStrings

function ep_plot()
    plot()
    plot!(xlabel=L"\sqrt{\rho {\sigma}^{3}}")
    plot!(ylabel=L"$-\beta{U_{\mathrm{ex}}}/N$ and $\beta{P}/\rho$")
    scatter!(sqrdat, nundat,
        label="Rasaiah et al. Energy", 
        markersize=8, markerstrokewidth=0, dpi=150)
    scatter!(sqrdat, cmpdat,
        label="Rasaiah et al. Pressure", 
        markersize=8, markerstrokewidth=0, dpi=150)
    plot!()
end

npt = 41
x = sqrt.(logrange(5e-4, 0.2, length=npt))

using SunlightHNC
ncomp = 2;
ng = 8192;
deltar = 0.01;
w = Wizard(; ncomp=ncomp, ng=ng, deltar=deltar)
model = restricted_primitive_model(w, lb / d);
# @show model

u = Float64[0.5, 0.5]
y = zeros(Float64, npt, 6)
for i in 1:npt
    rho = x[i]^2
    soln = solve!(model, rho .* u; closure="HNC")
    y[i, 1] = -soln.uex / rho
    y[i, 2] = soln.press / rho
end

p = ep_plot()
plot!(title="Fig.4.10, HNC, Hansen and McDonald, 4th ed. (2013)")
plot!(x, y[:, 1]; label="HNC Energy")
plot!(x, y[:, 2]; label="HNC Pressure")
plot!(ylims=(0.1, 1.3))
savefig("hm4_fig10-3_HNC-plot.png")
####################################

for i in 1:npt
    rho = x[i]^2
    soln = solve!(model, rho .* u; closure="MSA")
    y[i, 3] = -soln.uex / rho
    y[i, 4] = soln.press / rho
end

p = ep_plot()
plot!(title="Fig.4.10, MSA, Hansen and McDonald, 4th ed. (2013)")
#
plot!(x, y[:, 3]; label="MSA Energy")
plot!(x, y[:, 4]; label="MSA Pressure")
plot!(ylims=(0.1, 1.3))
savefig("hm4_fig10-3_MSA-plot.png")
####################################

p = ep_plot()
plot!(title="Fig.4.10, EXP, Hansen and McDonald, 4th ed. (2013)")

for i in 1:npt
    rho = x[i]^2
    w.lb = lb / d
    w.sigma = 1.0
    w.kappa = -1.0
    w.rho = rho * u
    w.cold_start = true
    hs_potential!(w)
    msa_solve!(w)
    save_reference(w)
    rpm_potential!(w)
    msa_solve!(w)
    exp_refine!(w)
    y[i, 5] = -w.uex / rho
    y[i, 6] = w.press / rho
end
plot!(x, y[:, 5]; label="EXP Energy")
plot!(x, y[:, 6]; label="EXP Pressure")
plot!(ylims=(0.1, 1.3))
savefig("hm4_fig10-3_EXP-plot.png")

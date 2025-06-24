# This is a Julia port of SunlightHNC 
# Copyright (C) 2025 Hiroharu Sugawara, Tokyo Metropolitan University <hsugawa@tmu.ac.jp>

# Based on an original code copyright (c) 2007 Lucian Anton.
# Modifications copyright (c) 2008, 2009 Andrey Vlasov.  
# Additional modifications copyright (c) 2009-2017 Unilever UK Central
# Resources Ltd (Registered in England & Wales, Company No 29140;
# Registered Office: Unilever House, Blackfriars, London, EC4P 4BQ,
# UK).

using SunlightHNC

w = Wizard()
w.arep[1, 1] = A = 25.0
dpd_potential!(w)

npt = 41
rho_max = 10.0

rho_s = LinRange(0, rho_max, npt)
y = zeros(npt)
for i in 1:npt
    w.rho[1] = rho = rho_s[i]
    hnc_solve!(w)
    y[i] = (w.press[] - rho) / (A * rho * rho)
end

# rho   (p-rho)/(A*rho^2)
data_text = """
  0.0, 0.0379935086163
  1.5, 0.0751786298043
  2.5, 0.0886823425022
  3.0, 0.0924251622846
  3.5, 0.0946639891655
  4.0, 0.0965259421847
  5.0, 0.0987451548125
  6.0, 0.0998358473824
  7.0, 0.1005510671090
  8.0,  0.102017933031"""

using DelimitedFiles
df = readdlm(IOBuffer(data_text), ',', skipstart=1)

using Plots
using LaTeXStrings
plot()
plot!(xlabel=L"$\rho$", ylabel=L"$(p-\rho)/(A \rho^{2})$")
plot!(rho_s, y, label="HNC")
scatter!(df[:, 1], df[:, 2],
    label="Groot & Warren (1997)",
    markersize=8, markerstrokewidth=0, dpi=150)
plot!(title="Fig 4, Groot and Warren (1997).")
savefig("gw_p-plot.png")

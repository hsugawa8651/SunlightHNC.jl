# This is a Julia port of SunlightHNC 
# Copyright (C) 2025 Hiroharu Sugawara, Tokyo Metropolitan University <hsugawa@tmu.ac.jp>

# Copyright (c) 2009-2019 Unilever UK Central Resources Ltd
# (Registered in England & Wales, Company No 29140; Registered Office:
# Unilever House, Blackfriars, London, EC4P 4BQ, UK).  Additional
# modifications copyright (c) 2020-2024 Patrick B Warren
# <patrick.warren@stfc.ac.uk> and STFC.


using Unitful
T = 300u"K" # temperature
εr = 78.5 # relative permittivity

lb = (1u"q")^2 / (4π * εr * 1u"ε0" * 1u"k" * T) |> u"Å"
@show lb

d = 5u"Å" # ion diameter
# concs = [0.5, 2, 5] * u"mol/L" # Molar units
# Molar units
concs = Quantity.([0.5, 2, 5], u"mol/L")

using Printf

println("Ion diameter= $(d)")
println("Bjerrum length = $(lb) = $(lb/d |> NoUnits) d")
println("T* = $(d/lb |> NoUnits)")

using SunlightHNC

ncomp = 2
# ng=2^14; deltar=1e-3
# ng=2^16; deltar=5e-4

w = Wizard(; ncomp=ncomp)
# w.verbose=true
model = restricted_primitive_model(w, lb / d)

# In Attard's figure the line for 0.5 M is surely h_{+-}
# imin = floor(Int64, 1.0 / w.deltar)
# imax = min(floor(Int64, 4.0 / w.deltar), length(w.r))
r_d = w.r * d .|> u"Å"
IR = 5.0u"Å" .<= r_d .<= 20u"Å"

using Plots
using LaTeXStrings
ax1 = plot(xlabel=L"$r / \mathrm{angstrom}$")
ax2 = plot(xlabel=L"$r / \mathrm{angstrom}$")

for (i, conc) in enumerate(concs)
    ρd3 = conc * 1u"Na" * d^3 |> NoUnits
    sol = solve!(model, Float64[ρd3, ρd3])

    for j in 1:2
        plot!(ax1, r_d[IR],
            sol.wizard.hr[IR, 1, j],
            label="hr$(j), " * string(conc),
            dpi=150, palette=:tab20)
    end

    h = sol.wizard.hr[IR, 1, (i == 1) ? 2 : 1]
    h[h.<0] .= 1e-20
    plot!(ax2, r_d[IR],
        log.(h), label="h, " * string(conc),
        dpi=150, palette=:tab10)

end
plot!(ax1, ylabel=L"hr")
plot!(ax2, ylabel=L"$\log{h}$")

plot!(ax1, ylims=(-0.5, 0.5))
savefig(ax1, "attard_fig1-plot1.png")
plot!(ax2, ylims=(-8, 1))
savefig(ax2, "attard_fig1-plot2.png")

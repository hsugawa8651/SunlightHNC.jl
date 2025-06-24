# This is a Julia port of SunlightHNC 
# Copyright (C) 2025 Hiroharu Sugawara, Tokyo Metropolitan University <hsugawa@tmu.ac.jp>

# Based on an original code copyright (c) 2007 Lucian Anton.
# Modifications copyright (c) 2008, 2009 Andrey Vlasov.  Additional
# modifications copyright (c) 2009-2017 Unilever UK Central Resources
# Ltd (Registered in England & Wales, Company No 29140; Registered
# Office: Unilever House, Blackfriars, London, EC4P 4BQ, UK).

using SunlightHNC
using Printf

function cr_press(w, drho)
    @show w.rho[1]
    p_xc = prev = 0.0e0
    n = floor(Int64, w.rho[1] / drho)
    w.cold_start = true
    for i in 1:n
        w.rho[1] = drho * i
        hnc_solve!(w)
        curr = w.comp_xc
        p_xc += 0.5 * drho * (prev + curr)
        prev = curr
    end
    return w.rho[1] + p_xc
end

function energy_aex(w, dA)
    aex_xc = prev = 0.0
    n = floor(Int64, w.arep[1, 1] / dA)
    w.cold_start = true
    for i in 1:n
        w.arep[1, 1] = dA * i
        dpd_potential!(w)
        hnc_solve!(w)
        curr = w.uex_xc / w.arep[1, 1]
        aex_xc += 0.5 * dA * (prev + curr)
        prev = curr
    end
    return w.uex_mf[] + aex_xc
end

function mu_aex(w, drho)
    aex = prev = 0.0e0
    n = floor(Int64, w.rho[1] / drho)
    w.cold_start = true
    for i in 1:n
        w.rho[1] = drho * i
        hnc_solve!(w)
        curr = w.muex[1]
        aex += 0.5 * drho * (prev + curr)
        prev = curr
    end
    return aex
end

ng = 8192
# ng = 65536
deltar = 0.02
w = Wizard(; ng=ng, deltar=deltar)

w.arep[1, 1] = A = 25.0
w.rho[1] = rho = 3.0
# w.verbose=true
dpd_potential!(w)
hnc_solve!(w)

write_thermodynamics(w)
write_params(w)

println("\n*** Example 6.1 ***")
println("rho = $(rho), A = $(A)")
@printf "pressure = %0.5f\n" w.press
@printf "energy density = %0.5f\n" w.uex
@printf "mu_excess = %0.5f\n" w.muex[1]

using Plots
using LaTeXStrings
plot(Series_gr(), w; rcut=3)
savefig("examples_6_1-plot1.png")

plot(Series_sk(), w; kcut=25.0)
savefig("examples_6_1-plot2.png")

println("\n*** Example 6.2 ***")
cr_pres = cr_press(w, 0.05)
@printf "CR pressure = %0.5f\n" cr_pres
@printf "VR pressure = %0.5f\n" w.press
@show abs(cr_pres - w.press) / w.press

println("\n*** Example 6.3 ***")
@printf "energy aex = %0.5f\n" energy_aex(w, 0.1)
@printf "mu aex     = %0.5f\n" mu_aex(w, 0.05)
@printf "direct aex = %0.5f\n" w.aex


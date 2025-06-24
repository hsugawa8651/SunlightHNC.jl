
# This is a Julia port of SunlightHNC 
# Copyright (C) 2025 Hiroharu Sugawara, Tokyo Metropolitan University <hsugawa@tmu.ac.jp>

# Based on an original code copyright (c) 2007 Lucian Anton.
# Modifications copyright (c) 2008, 2009 Andrey Vlasov.  Additional
# modifications copyright (c) 2009-2017 Unilever UK Central Resources
# Ltd (Registered in England & Wales, Company No 29140; Registered
# Office: Unilever House, Blackfriars, London, EC4P 4BQ, UK).  Later
# modifications copyright (c) 2020-2024 Patrick B Warren
# <patrick.warren@stfc.ac.uk> and STFC.

# soften_rpm : not implemented yet

"""
    struct Model{NT}

a model built on top of the grid
"""
struct Model{NT}
    nt::NT
end

"""
    Model(wizard::Wizard, name="HNC")

Construct Model object.
"""
function Model(wizard::Wizard, name="HNC")
    Model((;
        name,
        wizard,
    ))
end

write_params(m::Model) = write_params(m.nt.wizard)

# Python's ROUND(x) returns x with the fractional portion of its magnitude eliminated by rounding to the nearest whole number and with its sign preserved, converted to an INTEGER of the default kind.
ROUND(x) = copysign(floor(Int64, abs(x) + dp(0.5)), x)

"""
    additive_primitive_model

Construct additive primitive model.

# Arguments
* `wizard::Wizard`
* `lb::Number`
* `diam`
* `z` 
# Keyword arguments
* `name="additive PM"`
"""
function additive_primitive_model(wizard::Wizard, lb::Number, diam, z; name="additive PM")
    @assert 2 <= wizard.ncomp <= 4
    (; ncomp, nfnc, deltar, r, k) = wizard
    wizard.lb = lb
    diam_ave = zeros(ncomp, ncomp)
    for (i, j, ij) in upperTriangularIndices(ncomp)
        ave = (diam[i] + diam[j]) / 2
        wizard.dd[ij] = ave
        diam_ave[i, j] = ave
    end
    wizard.sigma = sigma = minimum(wizard.dd)
    cut = ROUND(sigma / deltar)
    for (i, j, ij) in upperTriangularIndices(ncomp)
        zzlb = z[i] * z[j] * lb
        wizard.ulong[:, ij] .= zzlb ./ r
        wizard.ulong[1:cut, ij] .= zzlb / sigma # cut-off inside min hard core (see docs)
        wizard.dulong[:, ij] .= -zzlb ./ r .^ 2
        wizard.dulong[1:cut, ij] .= 0.0e0 # -- ditto --
        wizard.ulongk[:, ij] .= 4 * π * zzlb * sin.(k * sigma) ./ (sigma * k .^ 3)
    end
    wizard.ushort .= 0.0e0
    wizard.dushort .= 0.0e0
    wizard.expnegus .= 1.0e0
    for ij in 1:nfnc
        cut = ROUND(wizard.dd[ij] / deltar)
        wizard.expnegus[1:cut, ij] .= dp(0.0)
    end
    wizard.u0 .= wizard.z .^ 2 * lb / sigma
    wizard.tp .= dp(0.0)
    wizard.tu .= dp(0.0)
    wizard.tl .= dp(0.0)

    Model((; name, wizard, lb, z, diam, diam_ave))
end

"""
    restricted_primitive_model

Construct restricted primitive model.

# Arguments
* `wizard::Wizard`
* `lb::Number`
"""
function restricted_primitive_model(wizard::Wizard, lb::Number)
    @assert wizard.ncomp == 2
    diam = Float64[1.0e0, 1.0e0]
    z = Float64[1.0e0, -1.0e0]
    additive_primitive_model(wizard, lb, diam, z; name="RPM")
end

"""
    function solve!(model::Model, rho; kwarg...)

Solve a `model`

# Arguments
* `model::Model`
* `rho`
# Keyword arguments
* `closure="HNC"`
* `alpha=0.2`
* `maxsteps=100`
* `cold_start=false`
"""
function solve!(model::Model, rho;
    closure="HNC", alpha=0.2, maxsteps=100, cold_start=false, # npic=6
)
    @show rho
    model.nt.wizard.rho = rho
    # model.wizard.npic = npic
    model.nt.wizard.alpha = alpha
    model.nt.wizard.maxsteps = maxsteps

    if !cold_start
        model.nt.wizard.cold_start = cold_start
    end
    if occursin("hnc", lowercase(closure))
        hnc_solve!(model.nt.wizard)
    elseif occursin("msa", lowercase(closure))
        msa_solve!(model.nt.wizard)
    else
        throw(ErrorException(
            "unrecognised closure $(closure) in solve; use HNC or MSA"))
    end
    return (;
        wizard=model.nt.wizard, # give access to underlying FORTRAN sector
        ERROR=model.nt.wizard.ERROR,
        deficit=model.nt.wizard.deficit,
        aex=model.nt.wizard.aex,
        press=model.nt.wizard.press,
        uex=model.nt.wizard.uex,
        muex=model.nt.wizard.muex,
        hc=model.nt.wizard.hc,
        hr=model.nt.wizard.hr,
        closure=model.nt.wizard.closure_name,
        nsteps=model.nt.wizard.nsteps)
end



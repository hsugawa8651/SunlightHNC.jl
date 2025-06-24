
# This is a Julia port of SunlightHNC 
# Copyright (C) 2025 Hiroharu Sugawara, Tokyo Metropolitan University <hsugawa@tmu.ac.jp>

# This file is part of SunlightDPD - a home for open source software
# related to the dissipative particle dynamics (DPD) simulation
# method.

# Based on an original code copyright (c) 2007 Lucian Anton.
# Modifications copyright (c) 2008, 2009 Andrey Vlasov.  Additional
# modifications copyright (c) 2009-2018 Unilever UK Central Resources
# Ltd (Registered in England & Wales, Company No 29140; Registered
# Office: Unilever House, Blackfriars, London, EC4P 4BQ, UK).  Later
# modifications copyright (c) 2020-2024 Patrick B Warren
# <patrick.warren@stfc.ac.uk> and STFC.

# SunlightDPD is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# SunlightDPD is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with SunlightDPD.  If not, see <http://www.gnu.org/licenses/>.

# Full documentation is found in the accompanying LaTeX document.

# dp = sizeof(1.0e0)
dp(x) = convert(Float64, x)
NINT(x) = floor(Int64, x + 0.5)
ispositive(x) = x > zero(x)

# Enumeration of current closures
@doc "NO\\_CLOSURE - constant for `closure_type`"
const NO_CLOSURE = 0
@doc "HNC\\_CLOSURE - constant for `closure_type`"
const HNC_CLOSURE = 1
@doc "RPA\\_CLOSURE - constant for `closure_type`"
const RPA_CLOSURE = 2
@doc "MSA\\_CLOSURE - constant for `closure_type`"
const MSA_CLOSURE = 3
@doc "EXP\\_CLOSURE - constant for `closure_type`"
const EXP_CLOSURE = 4

# Enumeration of current potential models
@doc "NO\\_MODEL\\_TYPE - constant for `charge_type`"
const NO_MODEL_TYPE = 0
@doc "DPD\\_GAUSSIAN\\_CHARGES - constant for `charge_type`"
const DPD_GAUSSIAN_CHARGES = 1
@doc "DPD\\_BESSEL\\_CHARGES - constant for `charge_type`"
const DPD_BESSEL_CHARGES = 2
@doc "DPD\\_LINEAR\\_CHARGES - constant for `charge_type`"
const DPD_LINEAR_CHARGES = 3
@doc "DPD\\_SLATER_APPROX\\_CHARGES - constant for `charge_type`"
const DPD_SLATER_APPROX_CHARGES = 4
@doc "DPD\\_SLATER_EXACT\\_CHARGES - constant for `charge_type`"
const DPD_SLATER_EXACT_CHARGES = 5
@doc "URPM\\_WITHOUT\\_USHORT - constant for `charge_type`"
const URPM_WITHOUT_USHORT = 6
@doc "URPM\\_WITH\\_USHORT - constant for `charge_type`"
const URPM_WITH_USHORT = 7
@doc "RPM\\_WITHOUT\\_USHORT - constant for `charge_type`"
const RPM_WITHOUT_USHORT = 8
@doc "RPM\\_WITH\\_USHORT - constant for `charge_type`"
const RPM_WITH_USHORT = 9
@doc "HARD\\_SPHERES - constant for `charge_type`"
const HARD_SPHERES = 10

# Enumeration of error codes
@doc "NO\\_ERROR -- constant for `return_code`"
const NO_ERROR = 0
@doc "CONVERGENCE\\_ERROR - constant for `return_code`"
const CONVERGENCE_ERROR = 1
@doc "AXEQB\\_ERROR - constant for `return_code`"
const AXEQB_ERROR = 2
@doc "DSYSV\\_ERROR - constant for `return_code`"
const DSYSV_ERROR = 3
@doc "DGEEV\\_ERROR - constant for `return_code`"
const DGEEV_ERROR = 4
@doc "MISMATCH\\_ERROR - constant for `return_code`"
const MISMATCH_ERROR = 5

# const (π) = dp(4.0) * atan(dp(1.0))
# const (2π) = dp(2.0) * (π)
# const (4π) = dp(4.0) * (π)
# const (√π) = sqrt((π))

@doc """
     mutable struct Wizard

holds lots of HNC data.

# Important fields for input
* Also refer to immutable arguments in [`Wizard()`](@ref)
# Immutable fields set by [`Wizard()`](@ref)
* `ncomp`: number of chemical components
* `nfnc` : number of functions = `ncomp (ncomp + 1) / 2`
* `r` : real space grid
* `k` : reciprocal space grid (computed)
* `deltak` : reciprocal space grid spacing (computed)
# Mutable fields, which may be modified before invoking a solver.
* `rho` : density array
* `z`: valence array
* `arep` :repulsion amplitude array
* `dd`: hard core diameter array
* `diam` : hard core diameter array, pairwise
* `sigma` : used both for the long-range Coulomb smearing length for the soft potentials, and for hard core diameter for RPM models.
* `sigmap` : +- long-range Coulomb smearing length (URPM)
* `kappa` : +- long-range Coulomb smoothing parameter (RPM)
* `rgroot` : linear charge smearing range (Groot)
* `lbda` : Slater charge smearing range (exact)
* `beta` : Slater charge smearing range (approx)
* `alpha`: Picard method, fraction of new solution
* `rc` : short-range DPD repulsion range
* `lb` : long-range Coulomb coupling length
# Computed and allocated fields, calculated by solvers.
* `c` : direct correlation functions (dcfs)
* `e` : indirect correlation functions (icfs)
* `ck` : transform of dcfs
* `ek` : transform of icfs
* `sk` : partial structure factors
* `h0` : reference total correlation functions
* `hr` : current total correlation functions
* `hc` : current contact values
* `ushort` : short range potential in real space
* `dushort` : derivative of `ushort`
* `ulong` : long range potential in real space
* `dulong` : derivative of `ulong`
* `ulongk` : long range potential in reciprocal space
* `expnegus` : `exp(-ushort)` (includes hard cores)
* `u0` : r = 0 of long range potential (diagonal)
* `muex`: chemical potential array
* `tp` : mean field pressure contribution
* `tu` : mean field energy contribution
* `tl` : mean field long range potential
# Important fields for output
* `ERROR` : difference between current and previous solutions
* `return_code` : return (error) code, any of [`return_code`](@ref return_codes)
* `closure_type` : any of [`closure_type`](@ref closure_types)
* `model_type` : which potential was (last) chosen
* `nsteps` : actual number of steps taken
# Working fields for FFT
* `planFFT`
* `fftw`
# Some fields are calculated by following methods
* [`make_pair_functions!`](@ref) 
* [`make_structure_factors!`](@ref) 
* [`make_thermodynamics!`](@ref) 
"""
mutable struct Wizard{J<:Integer,R<:Number}
    # Grid
    const ncomp::J
    const nfnc::J 
    const ng
    const nps
    rho 
    z
    muex 
    arep 
    diam
    dd 
    tp
    tu 
    tl 
    c
    e
    h0 
    hr
    hc
    ck
    ek
    sk 
    ushort
    dushort
    ulong 
    dulong
    ulongk 
    expnegus
    u0
    deltar::R 
    deltak::R
    r::AbstractVector{R}
    k::AbstractVector{R} 
    # Potential
    model_type::J 
    model_name::String
    alpha::R
    tol::R
    rc::R
    lb::R 
    sigma::R
    sigmap::R
    kappa::R
    rgroot::R
    lbda::R
    beta::R
    # Solver
    npic::J
    cold_start::Bool
    start_type::J
    istep::J
    nsteps::J
    maxsteps::J
    error_msg::String
    return_code::J # error code
    ERROR::R
    auto_fns::Bool
    # messages
    verbose::Bool
    silent::Bool
    suppress_msgs::Bool
    # FFT
    const planFFT
    fftw
    # Closure
    closure_type::J
    closure_name::String
    # Thermodynamics
    cf_mf::R
    cf_xc::R
    cf_gc::R
    pex::R 
    press::R
    comp::R
    comp_xc::R
    uex::R
    uex_mf::R
    uex_xc::R
    aex_rl::R
    aex_rs::R
    aex_kl::R
    aex_kd::R
    aex_ks::R
    aex::R
    deficit::R
end

"""
    Wizard(; kwargs...)

Construct [`Wizard`](@ref) object.
# Immutable arguments, each used as the struct's field as it is
* `ncomp::Int64=1` >= 1: number of chemical components
* `ng::Int64=2^8` > 0: grid size
* `nps::Int64=6` >= 2: number of previous states used in Ng method
* `npic::Int64=6` >= 2: number of Picard steps
* `deltar::Float64=1e-2` > 0 : real space grid spacing
# Mutable arguments, used to control solver
* `tol::Float64=1e-12` > 0 : Error tolerance for claiming convergence
* `maxsteps::Int64=100` : max number of steps to take for convergence
* `verbose::Bool=false` : Print solver diagnostics
* `silent::Bool=false` : Prevent printing warning/error messages
* `suppress_msgs::Bool=false` : In case of severe numerical instability
* `auto_fns:Bool=true` : Calculate stuff at end
* `cold_start` :  Force a cold start
* `start_type` : how to initialise in a cold start
"""
function Wizard(; ncomp=1, ng=2^8, nps=6, npic=6,
    deltar=1e-2, tol=1e-12, maxsteps=100, start_type=3,
    verbose=false, silent=false, suppress_msgs=false, auto_fns=true)

    ### Grids
    (ncomp >= 1) || throw(ArgumentError("ncomp=$(ncomp), though should be >= 1"))
    rho = @MVector zeros(ncomp)
    z = zeros(size(rho)...)
    muex = @MMatrix zeros(ncomp, ncomp)
    arep = zeros(size(muex)...)
    diam = zeros(size(muex)...)

    nfnc = ncomp * (ncomp + 1) ÷ 2
    dd = @MVector zeros(nfnc)
    tp = zeros(size(dd)...)
    tu = zeros(size(dd)...)
    tl = zeros(size(dd)...)

    (ng >= 2) || throw(ArgumentError("ng=$(ng), though should be >= 2"))
    (nps >= 2) || throw(ArgumentError("nps=$(nps), though should be >= 2"))
    c = zeros(ng - 1, nfnc, nps)
    e = zeros(size(c)...)
    h0 = zeros(ng - 1, nfnc)
    hr = zeros(ng - 1, ncomp, ncomp)
    hc = @MMatrix zeros(ncomp, ncomp)
    ck = zeros(size(h0)...)
    ek = zeros(size(h0)...)
    sk = zeros(size(hr)...)
    ushort = zeros(size(h0)...)
    expnegus = zeros(size(h0)...)
    dushort = zeros(size(h0)...)
    ulong = zeros(size(h0)...)
    dulong = zeros(size(h0)...)
    ulongk = zeros(size(h0)...)
    u0 = zeros(ncomp)

    # Default values

    # Make grids
    ispositive(deltar) || throw(ArgumentError("deltar=$(deltar), though should be > 0"))
    ispositive(tol) || throw(ArgumentError("tol=$(tol), though should be > 0"))

    deltak = π / (ng * deltar)
    r = deltar * (1:ng-1)
    k = deltak * (1:ng-1)

    alpha = 2e-1
    rc = 1.0e0
    lb = 0.0e0
    sigma = 1.0e0

    sigmap = 1.0e0
    kappa = -1.0e0
    rgroot = 1.0e0
    lbda = 1.0e0
    beta = 1.0e0
    (maxsteps > 0) || throw(ArgumentError("maxsteps=$(maxsteps), though should be > 0"))

    cold_start::Bool = true
    (1 <= start_type <= 3) || throw(ArgumentError("start_type=$(start_type), though should 1 <= <= 3"))
    start_type = 3

    ### Potential
    model_name = ""
    model_type = NO_MODEL_TYPE

    ### Solver
    ERROR = 0.0e0
    error_msg = ""
    return_code = istep = nsteps = 0

    # Make the FFTW plan
    fftw = Vector{Float64}(undef, ng - 1)
    planFFT = FFTW.plan_r2r(fftw, FFTW.RODFT00)

    ## Closure
    closure_type = 0
    closure_name = ""

    ## Thermodynamics
    cf_mf = cf_xc = cf_gc = 0.0e0
    pex = press = comp = comp_xc = 0.0e0
    uex = uex_mf = uex_xc = 0.0e0
    aex_rl = aex_rs = aex_kl = aex_kd = aex_ks = aex = 0.0e0
    deficit = 0.0

    Wizard(ncomp, nfnc, ng, nps,
        rho, z, muex, arep, diam, dd,
        tp, tu, tl, c, e, h0, hr, hc, ck, ek, sk,
        ushort, dushort, ulong, dulong, ulongk, expnegus, u0,
        deltar, deltak, r, k,
        model_type, model_name,
        alpha, tol, rc, lb, sigma, sigmap, kappa, rgroot, lbda, beta, npic,
        cold_start, start_type, istep, nsteps, maxsteps,
        error_msg, return_code, ERROR, auto_fns, verbose, silent, suppress_msgs,
        planFFT, fftw,
        closure_type, closure_name,
        cf_mf, cf_xc, cf_gc, pex, press, comp, comp_xc, uex, uex_mf, uex_xc,
        aex_rl, aex_rs, aex_kl, aex_kd, aex_ks, aex, deficit)
end


"""
    upperTriangularIndices(n::Int)

Traverse every element `[i,j]` of the upper triangular matrix
"""
function upperTriangularIndices(n::Int)
    Channel{Tuple{Int64,Int64,Int64}}(16) do ch
        ij = 1
        for j in 1:n
            for i in 1:j
                put!(ch, (i, j, ij))
                ij += 1
            end
        end
    end
end

"""
    upperTriangularIndices(w::Wizard)

Equivalent to `upperTriangularIndices(w.ncomp)`.
"""
upperTriangularIndices(w::Wizard) = upperTriangularIndices(w.ncomp)


"""
    upperTriangularIndices(n::Int, nlimit::Int)

Traverse every element `[i,j]` of the upper triangular matrix.
Restrict `i` and `j` <= `nlimit`
"""
function upperTriangularIndices(n::Int, nlimit::Int)
    (0 <= nlimit <= n) || throw(ArgumentError("nlimit=$(nlimit) should be within 0 .. $(n)"))
    Channel{Tuple{Int64,Int64,Int64}}(16) do ch
        ij = 1
        for j in 1:n
            for i in 1:j
                if i <= nlimit && j <= nlimit
                    put!(ch, (i, j, ij))
                end
                ij += 1
            end
        end
    end
end

"""
    upperTriangularIndices(w::Wizard, nlimit::Int)

Equivalent to `upperTriangularIndices(w.ncomp, nlimit)`.
"""
upperTriangularIndices(w::Wizard, nlimit::Int) = upperTriangularIndices(w.ncomp, nlimit)


"""
    diagonalUpperTriangularIndices(n::Int)

Traverse diagonal element [i,i] of the upper triangular matrix
"""
function diagonalUpperTriangularIndices(n::Int)
    Channel{Tuple{Int64,Int64}}(16) do ch
        ij = 1
        for j in 1:n
            for i in 1:j
                if i == j
                    put!(ch, (i, ij))
                end
                ij += 1
            end
        end
    end
end

"""
    diagonalUpperTriangularIndices(w::Wizard)

Equivalen to `diagonalUpperTriangularIndices(w.ncomp)`.
"""
diagonalUpperTriangularIndices(w::Wizard) = diagonalUpperTriangularIndices(w.ncomp)

"""
    write_params(w::Wizard)

Write out parameters for system and potentials.
"""
function write_params(w::Wizard)
    println(stdout, "=====================================================")
    println(stdout, "GRID DETAILS")
    println(stdout, " ng = ", w.ng, " ncomp = ", w.ncomp, " nfnc = ", w.nfnc, " nps = ", w.nps)
    println(stdout, " deltar = ", w.deltar, " deltak = ", w.deltak)
    println(stdout, " deltar*deltak*ng/pi = ", w.deltar * w.deltak / (π) * w.ng)
    println(stdout, " r(ng-1) = ", w.r[end], " k(ng-1) = ", w.k[end])
    println(stdout, "POTENTIAL DETAILS :: model type =", w.model_type, " ", w.model_name)
    if w.model_type == NO_MODEL_TYPE
        println(stdout, "No potential has been selected / potential set externally")
    elseif w.model_type <= DPD_SLATER_EXACT_CHARGES
        print(stdout, "DPD potential was selected, matrix A = ")
        # println(stdout, w.arep)
        # println(stdout, replace(repr("text/plain", w.arep), "\n" => raw"\n"))
        display(w.arep)
        # for i = 1:w.ncomp
        #    println(stdout, " ", w.arep[i, :])
        # end
        println(stdout, " valencies, z = ", w.z)
        println(stdout, " rc = ", w.rc, " lb = ", w.lb)
        if w.model_type == DPD_GAUSSIAN_CHARGES
            println(stdout, " Gaussian smearing, sigma = ", w.sigma)
        end
        if w.model_type == DPD_BESSEL_CHARGES
            println(stdout, " Bessel smearing, sigma = ", w.sigma)
        end
        if w.model_type == DPD_LINEAR_CHARGES
            println(stdout, " linear (Groot) smearing, R = ", w.rgroot)
            println(stdout, " equivalent Gaussian sigma = ",
                sqrt(dp(2.0) / dp(15.0)) * w.rgroot)
        end
        if w.model_type == DPD_SLATER_APPROX_CHARGES
            println(stdout, " Slater smearing (approx), beta = ", beta)
            println(stdout,
                " 1 / beta = ", dp(1.0) / w.beta,
                ", 5 / (8 beta) = ", dp(0.625) / w.beta)
            println(stdout,
                " equivalent Gaussian sigma = ", sqrt(dp(0.5)) / w.beta)
        end
        if w.model_type == DPD_SLATER_EXACT_CHARGES
            println(stdout, " Slater smearing (exact), lambda = ", w.lbda)
            println(stdout, " equivalent Gaussian sigma = ", w.lbda)
        end
    elseif w.model_type <= URPM_WITH_USHORT
        if w.model_type == URPM_WITHOUT_USHORT
            println(stdout, "URPM potential was selected with ushort unused")
        else
            println(stdout, "URPM potential was selected with ushort used")
        end
        println(stdout, "w.lb = ", w.lb, " sigma = ", w.sigma, " sigmap = ", w.sigmap)
    elseif w.model_type <= RPM_WITH_USHORT
        if w.model_type == RPM_WITHOUT_USHORT
            println(stdout, "RPM potential was selected with ushort unused")
        else
            println(stdout, "RPM potential was selected with ushort used")
        end
        if w.kappa < dp(0.0)
            println(stdout, "w.lb = ", w.lb, " sigma = ", w.sigma, " kappa -> infinity")
        else
            println(stdout, "w.lb = ", w.lb, " sigma = ", w.sigma, " kappa = ", w.kappa)
        end
        println(stdout, " valencies, z = ", w.z)
        print(stdout, "hard core diameter array")
        #for i = 1:w.ncomp
        #    println(stdout, " ", w.diam[i, :])
        #end
        display(w.diam)
        println(stdout, " hard sphere diameters, dd = ", w.dd)
    elseif w.model_type == HARD_SPHERES
        println(stdout, "Hard sphere potential with sigma = ", w.sigma)
        println(stdout, "hard core diameter array")
        for i = 1:w.ncomp
            println(stdout, " ", w.diam[i, :])
        end
        println(stdout, " hard sphere diameters, dd = ", w.dd)
    else
        println(stdout, "Undefined potential")
    end
    println(stdout, "SYSTEM DETAILS")
    println(stdout, " rho = ", w.rho)
    sum_rho = sum(w.rho)
    if ispositive(sum_rho)
        println(stdout, " x = ", w.rho[:] / sum_rho)
    end
    println(stdout, " sum(rho) = ", sum_rho)
    println(stdout, " sum(rho*z) = ", sum(w.rho .* w.z))
    println(stdout, "=====================================================")
end


"""
    dpd_potential!(w::Wizard, charge_type=DPD_GAUSSIAN_CHARGES)

Build the potential arrays.

# Arguments
* `charge_type` : Use [the defined integer constants](@ref charge_types) for these.
  * `DPD_GAUSSIAN_CHARGES=1`, Gaussian (default)
  * `DPD_BESSEL_CHARGES=2`, Bessel
  * `DPD_LINEAR_CHARGES=3`, Groot (linear)
  * `DPD_SLATER_APPROX_CHARGES=4`, Slater (approx)
  * `DPD_SLATER_EXACT_CHARGES=5`, Slater (exact) 

# Referred fields in Wizard object
A factor `beta` = 1/kT is implicit in the following definitions.

## Fields referred by all `charge_type`
* `arep[:,:]` 
* `z`
"""
function dpd_potential!(w::Wizard, charge_type=DPD_GAUSSIAN_CHARGES)

    if ispositive(charge_type)
        (charge_type <= HARD_SPHERES) || throw(ArgumentError("charge_type=$(charge_type), though should be 0 <= <= $(HARD_SPHERES)"))
        w.model_type = charge_type
    end

    # w.model_type = ctype
    # Sort out some recoded potential parameters.  Also force the
    # amplitude matrix to be symmetric, set by upper triangular
    # entries.  This enforces the rule in the documentation and also
    # simplifies the printing.

    aa = zeros(w.nfnc)
    aa .= NaN
    zz = zeros(w.nfnc)
    zz .= NaN

    @views for (i, j, ij) in upperTriangularIndices(w)
        aa[ij] = w.arep[i, j]
        zz[ij] = w.z[i] * w.z[j]
        if (i < j)
            w.arep[j, i] = w.arep[i, j]
        end
    end

    irc = NINT(w.rc / w.deltar)

    # Leave out the amplitude, then the function can be re-used
    # (see below)

    @views begin
        w.ushort[:, w.nfnc] .= dp(0.0)
        w.ushort[1:irc, w.nfnc] .= ((dp(1.0) .- w.r[1:irc] / w.rc) .^ 2) / 2
        w.dushort[:, w.nfnc] .= dp(0.0)
        w.dushort[1:irc, w.nfnc] .= -(dp(1.0) .- w.r[1:irc] / w.rc) / w.rc
    end
    local ulong0
    # Gaussian charges

    if w.model_type == DPD_GAUSSIAN_CHARGES

        ulong0 = w.lb / (w.sigma * (√π))
        @views begin
            w.ulong[:, w.nfnc] .=
                (w.lb ./ w.r) .* erf.(dp(0.5) * w.r / w.sigma)
            w.ulongk[:, w.nfnc] .=
                ((4π) * w.lb ./ w.k .^ 2) .* exp.(-w.k .^ 2 * w.sigma^2)
            w.dulong[:, w.nfnc] .=
                w.lb * exp.(dp(-0.25) * w.r .^ 2 / w.sigma^2) ./
                ((√π) * w.r * w.sigma) .-
                w.lb * erf.(dp(0.5) * w.r / w.sigma) ./ w.r .^ 2
        end
        w.model_name = "DPD_with_Gaussian_charges"

    end

    # Bessel charges

    if w.model_type == DPD_BESSEL_CHARGES

        ulong0 = w.lb / w.sigma
        @views begin
            w.ulong[:, w.nfnc] .=
                (w.lb ./ w.r) .* (dp(1.0) .- exp.(-w.r / w.sigma))
            w.ulongk[:, w.nfnc] .=
                ((4π) * w.lb ./ w.k .^ 2) .*
                dp(1.0) ./ (dp(1.0) .+ (w.k .^ 2) * w.sigma^2)

            w.dulong[:, w.nfnc] .=
                -(w.lb ./ w.r .^ 2) .*
                (dp(1.0) .- exp.(-w.r / w.sigma) .* (dp(1.0) .+ w.r ./ w.sigma))
        end
        w.model_name = "DPD_with_Bessel_charges"

    end

    # Linear charge smearing as in Groot [JCP v118, 11265 (2003)].
    # Note we do not give the real space part here hence the
    # thermodynamic calculations will be wrong.  This could be fixed
    # up by doing a numerical FFT of the potential.  Best would be to
    # separate off the long range 4 piw.lb / k^2 first.  TO BE DONE!

    if w.model_type == DPD_LINEAR_CHARGES

        w.ulong[:, w.nfnc] .= dp(0.0)

        ulong0 = dp(0.0)

        w.ulongk[:, w.nfnc] .=
            ((4π) .* w.lb ./ (w.k .^ 2)) * dp(144.0) .*
            (dp(2.0) .-
             dp(2.0) * cos.(w.k * w.rgroot) .-
             w.k * w.rgroot .* sin.(w.k * w.rgroot)) .^ 2 ./
            ((w.k .^ 8) * (w.rgroot^8))

        w.dulong[:, w.nfnc] .= dp(0.0)

        w.model_name = "DPD_with_linear_charges"

    end

    # Slater charge smearing as in Gonzales-Melchor et al, [JCP v125,
    # 224107 (2006)] but here with exact expression for interaction.

    if w.model_type == DPD_SLATER_EXACT_CHARGES

        w.ulong[:, w.nfnc] .=
            (w.lb ./ w.r) .*
            (dp(1.0) .- exp.(dp(-2.0) * w.r / w.lbda) .*
                        (dp(1.0) .+ dp(1.375) * w.r / w.lbda .+
                         dp(0.75) * (w.r .^ 2) / w.lbda^2 .+
                         (w.r .^ 3) / (dp(6.0) * w.lbda^3)))

        ulong0 = dp(0.625) * w.lb / w.lbda

        w.ulongk[:, w.nfnc] .=
            ((4π) * w.lb ./ (w.k .^ 2)) * dp(1.0) ./
            (dp(1.0) .+ (w.k .^ 2) * w.lbda^2 / dp(4.0)) .^ 4

        w.dulong[:, w.nfnc] .=
            -(w.lb ./ w.r .^ 2) .*
            (dp(1.0) .- exp.(dp(-2.0) * (w.r) / w.sigma) .*
                        (dp(1.0) .+ dp(2.0) * w.r / w.lbda .+
                         dp(2.0) * (w.r .^ 2) / w.lbda^2 .+
                         dp(7.0) * (w.r .^ 3) / (dp(6.0) * w.lbda^3)
                         .+
                         (w.r .^ 4) / (dp(3.0) * w.lbda^4)))

        w.model_name = "DPD_with_Slater_(exact)_charges"

    end

    # Slater charge smearing as in Gonzales-Melchor et al, [JCP v125,
    # 224107 (2006)] with original approximate expression.

    if w.model_type == DPD_SLATER_APPROX_CHARGES

        w.ulong[:, w.nfnc] .=
            (w.lb ./ w.r) .*
            (dp(1.0) .- exp.(-2 * w.beta * w.r) .* (
                dp(1.0) .+ w.beta * w.r))

        ulong0 = w.beta * w.lb

        w.ulongk[:, w.nfnc] .=
            ((4π) * w.lb ./ (w.k .^ 2)) * dp(1.0) ./
            (dp(1.0) .+ (w.k .^ 2) ./
                        (dp(4.0) * w.beta^2)) .^ 2

        w.dulong[:, w.nfnc] .=
            -(w.lb ./ w.r .^ 2) .*
            (dp(1.0) .- exp.(dp(-2.0) * w.beta * w.r) .*
                        (dp(1.0) .+ dp(2.0) * w.beta * w.r .* (dp(1.0) .+ w.beta * w.r)))

        w.model_name = "DPD_with_Slater_(approx)_charges"

    end

    # Generate the pair potentials by walking through the functions.
    # In the final step we correctly normalise the final function.

    for i = 1:w.nfnc
        w.ushort[:, i] .= aa[i] .* w.ushort[:, w.nfnc]
        w.dushort[:, i] .= aa[i] .* w.dushort[:, w.nfnc]
        w.ulong[:, i] .= zz[i] .* w.ulong[:, w.nfnc]
        w.ulongk[:, i] .= zz[i] .* w.ulongk[:, w.nfnc]
        w.dulong[:, i] .= zz[i] .* w.dulong[:, w.nfnc]
    end

    # These individual species-pair contributions to the mean field
    # compressibility factor and the mean-field internal energy per
    # particle can be calculated analytically for the DPD potential.

    w.tp .= (π) * w.rc .^ 3 * aa / dp(30.0)
    w.tu .= w.tp
    w.tl .= dp(0.0)

    # The r = 0 diagonal parts of the long range potential

    w.u0 .= w.z .^ 2 * ulong0

    # Generate auxiliary function

    w.expnegus .= exp.(-w.ushort)
end

"""
    urpm_potential!(w::Wizard, use_ushort=false)

Build the potential arrays for the softened URPM (Gaussian charges),
This expects `ncomp = 2`, and will set `z[1] = 1`, `z[2] = -1`.  

# Arguments
* `use_ushort::Bool=false`: whether ushort is used or not.

# Referred fields for input
* `sigma`
* `sigmap`
* `lb`
"""
function urpm_potential!(w::Wizard, use_ushort=false)

    if w.ncomp == 2
    else
        w.return_code = MISMATCH_ERROR
        # w.error_msg = "mismatch ncomp <> 2 in urpm_potential"
        msg = "ncomp = $(w.ncomp), though urpm_potential requires 2."
        w.error_msg = msg
        if (!w.silent)
            println(stdout, "** error: ", w.error_msg)
        end
        # return nothing
        throw(ArgumentError(msg))

    end

    w.z[1] = dp(1.0)
    w.z[2] = dp(-1.0)

    w.ulong[:, 1] .= w.lb * erf.(dp(0.5) * w.r / w.sigma) ./ w.r

    w.ulongk[:, 1] .= (4π) * w.lb * exp.(-w.k .^ 2 * w.sigma^2) / w.k .^ 2

    w.dulong[:, 1] .=
        w.lb * exp.(dp(-0.25) * w.r .^ 2 / w.sigma^2) ./ ((√π) * w.r * w.sigma) .-
        w.lb * erf.(dp(0.5) * w.r / w.sigma) ./ r .^ 2

    w.ulong[:, 2] .= -w.lb * erf.(dp(0.5) * w.r / w.sigmap) ./ r

    w.ulongk[:, 2] .= -(4π) * w.lb * exp.(-w.k .^ 2 * w.sigmap^2) ./ w.k .^ 2

    w.dulong[:, 2] .=
        -lb * exp(dp(-0.25) * w.r .^ 2 / w.sigmap^2) ./
        ((√π) * w.r * w.sigmap) .+
        w.lb * erf(dp(0.5) * w.r / w.sigmap) / r .^ 2

    w.ulong[:, 3] .= w.ulong[:, 1]
    w.ulongk[:, 3] .= w.ulongk[:, 1]
    w.dulong[:, 3] .= w.dulong[:, 1]

    w.ushort[:, :] .= dp(0.0)
    w.dushort[:, :] .= dp(0.0)

    if use_ushort
        w.ushort[:, 2] .= w.ulong[:, 2] + w.ulong[:, 1]
        w.dushort[:, 2] .= w.dulong[:, 2] + w.dulong[:, 1]
        w.ulong[:, 2] .= -w.ulong[:, 1]
        w.ulongk[:, 2] .= -w.ulongk[:, 1]
        w.dulong[:, 2] .= -w.dulong[:, 1]
    end

    # These individual species-pair contributions to the mean field
    # compressibility factor and the mean-field internal energy per
    # particle can be calculated analytically for the URPM potential.
    # These are the same whether using ushort or not, as they are
    # defined in terms of the total potential.

    w.tp[1] = dp(0.0)
    w.tp[2] = (2π) * w.lb * (w.sigmap^2 - w.sigma^2)
    w.tp[3] = dp(0.0)

    w.tu .= w.tp

    # If not using ushort, we are off the symmetry point condition and
    # the contribution of the long range part should be incorporated
    # into the compressibility and chemical potential expressions.

    if use_ushort
        w.tl .= dp(0.0)
    else
        w.tl .= dp(2.0) * w.tp
    end

    # The r = 0 diagonal parts of the long range potential

    w.u0 .= w.z .^ 2 * w.lb / (w.sigma * (√π))

    # Generate auxiliary function

    w.expnegus .= exp.(-w.ushort)

    if use_ushort
        w.model_type = URPM_WITH_USHORT
        w.model_name = "URPM_with_U_short"
    else
        w.model_type = URPM_WITHOUT_USHORT
        w.model_name = "URPM_without_U_short"
    end
end

"""
    rpm_potential!(w::Wizard, use_ushort=false)

Build the potential arrays for the softened RPM (charged hard
spheres).

Using `kappa < 0` implies the pure RPM case (`kappa` -> infinity). 

The case `ncomp = 3` corresponds to the RPM in the presence of a neutral hard sphere solvent.

# Arguments
* `use_ushort::Bool=false`: whether ushort is used or not.

# Referred fields for input
* `diam`
* `sigma`
* `lb`
* `kappa`
"""
function rpm_potential!(w::Wizard, use_ushort=false)

    if 2 <= w.ncomp <= 3
        # do nothing
    else
        w.return_code = MISMATCH_ERROR
        # msg = "mismatch ncomp <> 2 or 3 in rpm_potential"
        msg = "ncomp = $(w.ncomp), though rpm_potential requires 2 or 3."
        w.error_msg = msg
        if (!w.silent)
            println(stdout, "** error: ", w.error_msg)
        end
        # return nothing
        throw(ArgumentError(msg))
    end

    w.z[1] = dp(1.0)
    w.z[2] = dp(-1.0)
    if (w.ncomp == 3)
        w.z[3] = dp(0.0)
    end

    # Set the diameters from an array, or sigma if none are preset.
    # The remaining logic means sigma, if not set, is set to the
    # smallest diameter

    if any(ispositive.(w.diam))
        @views for (i, j, ij) in upperTriangularIndices(w)
            w.dd[ij] = w.diam[i, j]
        end
        if w.sigma == dp(0.0)
            w.sigma = minimum(w.dd)
        end
    else
        w.dd .= w.sigma
    end

    # Zero everything to start with, as for hard spheres

    w.ulong .= dp(0.0)
    w.ulongk .= dp(0.0)
    w.dulong .= dp(0.0)
    w.ushort .= dp(0.0)
    w.dushort .= dp(0.0)

    irc = NINT(w.sigma / w.deltar)

    w.ulong[:, 1] .= w.lb ./ w.r
    w.ulong[1:irc, 1] .= w.lb / w.sigma
    w.ulongk[:, 1] .= (4π) * w.lb * sin.(w.k * w.sigma) ./ (w.sigma * w.k .^ 3)
    w.dulong[:, 1] .= -w.lb ./ w.r .^ 2
    w.dulong[1:irc, 1] .= dp(0.0)

    if ispositive(w.kappa)
        w.ulong[:, 2] = -w.lb * erf.(w.kappa * w.r) ./ w.r
        w.ulongk[:, 2] = -(4π) * w.lb * exp.(-w.k .^ 2 / (dp(4.0) * w.kappa^2)) / (w.k .^ 2)
        w.dulong[:, 2] = -dp(2.0) * w.kappa * w.lb * exp.(-w.kappa^2 .* w.r .^ 2) / ((√π) * w.r) +
                         w.lb * erf.(w.kappa * w.r) ./ (w.r .^ 2)
    else
        w.ulong[:, 2] .= -w.ulong[:, 1]
        w.ulongk[:, 2] .= -w.ulongk[:, 1]
        w.dulong[:, 2] .= -w.dulong[:, 1]
    end

    w.ulong[:, 3] .= w.ulong[:, 1]
    w.ulongk[:, 3] .= w.ulongk[:, 1]
    w.dulong[:, 3] .= w.dulong[:, 1]

    w.ushort[:, :] .= dp(0.0)
    w.dushort[:, :] .= dp(0.0)

    if use_ushort
        w.ushort[:, 2] .= w.ulong[:, 2] + w.ulong[:, 1]
        w.dushort[:, 2] .= w.dulong[:, 2] + w.dulong[:, 1]
        w.ulong[:, 2] .= -w.ulong[:, 1]
        w.ulongk[:, 2] .= -w.ulongk[:, 1]
        w.dulong[:, 2] .= -w.dulong[:, 1]
    end

    # Generate auxiliary function - this is where HS diam implemented

    w.expnegus = exp.(-w.ushort)

    for ij = 1:w.nfnc
        irc = NINT(w.dd[ij] / w.deltar)
        w.expnegus[1:irc, ij] .= dp(0.0)
    end

    # These are the analytic contributions to the thermodynamics.

    w.tp .= dp(0.0)
    w.tu .= dp(0.0)
    w.tl .= dp(0.0)

    if ispositive(w.kappa)
        w.tp[2] = (π) * w.lb *
                  (w.sigma * exp(-w.kappa^2 * w.sigma^2) /
                   (w.kappa * (√π)) +
                   (1 / (dp(2.0) * w.kappa^2) - w.sigma^2 / dp(3.0)) *
                   erfc(w.kappa * w.sigma))
        w.tu[2] = (π) * w.lb *
                  (w.sigma * exp(-w.kappa^2 * w.sigma^2) /
                   (w.kappa * (√π)) +
                   (1 / (dp(2.0) * w.kappa^2) - w.sigma^2) * erfc(w.kappa * w.sigma))
        if !use_ushort # off SYM condition;
            # missing last term added in v1.11
            w.tl[2] = (π) * w.lb / w.kappa^2 -
                      w.lb * w.sigma^2 * (2π) / dp(3.0)
        end
    end

    # The r = 0 diagonal parts of the long range potential

    w.u0 .= w.z .^ 2 * w.lb / w.sigma

    if use_ushort
        w.model_type = RPM_WITH_USHORT
        w.model_name = "RPM_with_U_short"
    else
        w.model_type = RPM_WITHOUT_USHORT
        w.model_name = "RPM_without_U_short"
    end
end

"""
    hs_potential!(w::Wizard)

Build the potential arrays for hard spheres with diameter `sigma`.
This works for arbitrary numbers of components.

# Referred fields for input
* `diam`
* `sigma`
"""
function hs_potential!(w::Wizard)
    # Set the hard core diameters from an array, or sigma if none are pre-set

    if any(ispositive.(w.diam))
        @views for (i, j, ij) in upperTriangularIndices(w)
            w.dd[ij] = w.diam[i, j]
        end
    else
        w.dd .= w.sigma
    end

    w.ulong .= dp(0.0)
    w.ulongk .= dp(0.0)
    w.dulong .= dp(0.0)
    w.ushort .= dp(0.0)
    w.dushort .= dp(0.0)

    # This is where the hard sphere diameters enter

    w.expnegus .= dp(1.0)

    for ij = 1:w.nfnc
        irc = NINT(w.dd[ij] / w.deltar)
        w.expnegus[1:irc, ij] .= dp(0.0)
    end

    # These are the analytic contributions to the thermodynamics.

    w.tp .= dp(0.0)
    w.tu .= dp(0.0)
    w.tl .= dp(0.0)
    w.u0 .= dp(0.0)

    w.model_type = HARD_SPHERES
    w.model_name = "hard_spheres"
end

"""
    oz_solve!(w::Wizard)

Solve the Ornstein-Zernicke equation 
to determine `e = h - c`, given `c`.  We re-partition the long range part of the
potential so that the routine actually calculates `c' = c + Ulong` and `e' = e - Ulong`.  
This is because the Fourier transform of `Ulong` can be computed in closed form.  
Note `h = e + c = e' + c'`.
"""
function oz_solve!(w::Wizard)

    i1 = mod(w.istep - 1, w.nps) + 1

    # Forward transform the real space functions c, to the reciprocal
    # space functions ck.

    fftwx = similar(w.fftw)
    fftwy = similar(w.fftw)
    for i = 1:w.nfnc
        fftwx .= w.r .* w.c[:, i, i1]
        fftwy .= w.planFFT * fftwx
        w.ck[:, i] .= (2π) * w.deltar * fftwy ./ w.k
    end

    if w.ncomp == 1

        # In the one component case OZ inversion is straightforward.
        # Note the implicit indexing on wavevector k.

        w.ek[:, 1] .=
            (w.ck[:, 1] .- w.ulongk[:, 1]) ./
            (dp(1.0) .- w.rho[1] * (w.ck[:, 1] .- w.ulongk[:, 1])) .- w.ck[:, 1]

    else # Multicomponent OZ inversion
        # irc = 0
        # perm = zeros(Int64, ncomp)
        a = zeros(w.ncomp, w.ncomp)
        # b = similar(a)
        # X = similar(a)

        cmat = similar(a)
        umat = similar(a)

        # First set up a unit matrix, and the diagonal R matrix -- see
        # the documentation for the math here.

        rhomat = Diagonal(w.rho)

        # Do the matrix calculations for each wavevector k.

        for ik = 1:w.ng-1

            # Unpack the reciprocal space functions into matrices.

            @views for (i, j, ij) in upperTriangularIndices(w)
                cmat[i, j] = w.ck[ik, ij]
                umat[i, j] = w.ulongk[ik, ij]
                if i < j
                    cmat[j, i] = cmat[i, j]
                    umat[j, i] = umat[i, j]
                end
            end

            # The following matrices are constructed:
            #   M0 = (C - beta UL) . R
            #   A = I - (C - beta UL) . R = I - M0
            #   B = (C - beta UL) . R . C - beta UL = M0 . R . C - beta UL
            # (note that beta = 1/kT = 1 is not written explicitly)

            m0 = (cmat - umat) * rhomat
            a = LinearAlgebra.I - m0
            b = (m0 * cmat) - umat

            # Solve the equation A.X = B so that
            # X = [I - (C - beta U) . R]^(-1) . [(C - beta U) . R . C - beta U]

            local X
            try
                F = lu(a)
                X = F \ b
            catch ex
                if isa(ex, SingularException)
                    w.return_code = AXEQB_ERROR
                    w.error_msg = "axeqb encountered singular problem in oz_solve"
                    if !w.suppress_msgs && !w.silent
                        println(stdout, "** error: ", w.error_msg)
                        println(stdout, "** further messages of this kind will be suppressed")
                        w.suppress_msgs = true
                    end
                    return nothing
                else
                    rethrow(ex)
                end
            end

            # Now X(I, :) = B(PERM(I), :) is the new estimate for the
            # reciprocal space functions ek.  They are built from the
            # upper triangle of the matrix.

            @views for (i, j, ij) in upperTriangularIndices(w)
                w.ek[ik, ij] = X[i, j]
            end

        end # loop over k vectors

    end # select single component or multicomponent case

    # Do the Fourier back transforms

    fftwx = similar(w.fftw)
    fftwy = similar(w.fftw)

    for i = 1:w.nfnc
        fftwx .= w.k .* w.ek[:, i]
        fftwy .= w.planFFT * fftwx
        w.e[:, i, i1] .= w.deltak / ((2π)^2) * fftwy ./ w.r
    end

end

"""
    oz_solve2!(w::Wizard)

Solve an alternate version of the Ornstein-Zernicke equation 
to determine `c` and `e` from the reference `h0`.  
In practice as always we actually 
calculate `c' = c + Ulong` and `e' = e - Ulong`.
Note `h = e + c = e' + c'`.  
The result is saved to position 1 in the history trajectory.  
Note that the offset `ulongk` in.
"""
function oz_solve2!(w::Wizard)

    hk = zeros(w.ng - 1, w.nfnc)
    fftwx = similar(w.fftw)
    fftwy = similar(w.fftw)
    for i = 1:w.nfnc
        fftwx .= w.r .* w.h0[:, i]
        fftwy .= w.planFFT * fftwx
        hk[:, i] .= ((2π) * w.deltar) * fftwy ./ w.k
    end

    if w.ncomp == 1

        # In the one component case the OZ solution is trivial.
        w.ck[:, 1] .= hk[:, 1] ./ (dp(1.0) .+ w.rho[1] * hk[:, 1])

    else # Multicomponent OZ solution
        # irc = 0
        # perm = zeros(Int64, ncomp)
        X = zeros(w.ncomp, w.ncomp)
        # As above set up a unit matrix, and the diagonal R matrix
        rhomat = Diagonal(w.rho)

        # Do the matrix calculations for each wavevector k.
        # hmat_parent = zeros(w.ncomp, w.ncomp)
        hmat_parent = similar(X)
        hmat = Symmetric(hmat_parent)

        for ik = 1:w.ng-1

            # Convert the reciprocal space functions into matrices.

            @views for (i, j, ij) in upperTriangularIndices(w)
                hmat_parent[i, j] = hk[ik, ij]
                #if (i < j)
                #    hmat[j, i] = hmat[i, j]
                #end
            end

            # Construct A = I + H . R, and B = H
            a = LinearAlgebra.I + (hmat * rhomat)
            b = hmat

            # Solve A.X = B so that X = (I + H.R)^(-1) . H.
            # axeqb_reduce(a, ncomp, b, ncomp, perm, irc)
            try
                F = lu(a)
                X = F \ b
            catch ex
                if isa(ex, SingularException)
                    w.return_code = AXEQB_ERROR
                    w.error_msg = "axeqb encountered singular problem in in oz_solve2"
                    if !w.suppress_msgs && !w.silent
                        println(stdout, "** error: ", w.error_msg)
                        println(stdout, "** further messages of this kind will be suppressed")
                        w.suppress_msgs = true
                    end
                    return nothing
                elseif isa(ex, DimensionMismatch)
                    exit()
                else
                end
            end

            # Now compute C = (I + H.R)^(-1) . H + beta UL
            # (map back to functions, and unravel the pivoting)
            # The + beta UL is done afterwards.

            @views for (i, j, ij) in upperTriangularIndices(w)
                w.ck[ik, ij] = X[i, j]
            end

        end

    end

    # Do the Fourier back transforms.

    fftwx = similar(w.fftw)
    fftwy = similar(w.fftw)
    for i = 1:w.nfnc
        fftwx .= w.k .* w.ck[:, i]
        fftwy .= w.planFFT * fftwx
        w.c[:, i, 1] .= w.deltak / ((2π)^2) * fftwy ./ w.r
    end

    # Add + beta UL (since we have the exact real and reciprocal space
    # functions we don't need to pass this through the DDFT).
    # Strictly speaking this is a flourish that isn't necessary for
    # the structural features (pair functions and structure factors)
    # since it cancels out again.  However several of the
    # thermodynamic calculations assume this offset is present.

    w.ck .+= w.ulongk
    w.c[:, :, 1] .+= w.ulong

    # Recover the indirect correlation functions.  This means the
    # results can be used in the structure and thermodynamics routines
    # as though they had come from the HNC/MSA/RPA solution.

    w.ek .= hk - w.ck
    w.e[:, :, 1] .= w.h0[:, :] - w.c[:, :, 1]

end

"""
    hnc_picard!(w::Wizard)

This method implements the HNC closure expressed 
as `c = exp( -beta v + e) - e - 1` 
where `e = h - c` is the indirect correlation function, 
`c` is the direct correlation function from the Ornstein-Zernicke
relation, `h = g - 1`, and `g` is the pair distribution function.  
One can show this is equivalent to `g = exp(-v + h - c)` in Hansen + McDonald.  
As above, the routine actually works 
with `c' = c + Ulong` and `e' = e - Ulong` 
where `Ulong` is the long-range part of the
potential for which the Fourier transform is simple.  This means
that `v` in the above expression is the short-range part of the
potential only.  As we iterate we move forward in the history
trajectory, cyclically.
"""
function hnc_picard!(w::Wizard)
    w.istep += 1
    i1 = mod(w.istep - 1, w.nps) + 1
    # i0 = i1 - 1
    # i0 = (i0 == 0) ? nps : i0
    i0 = (i1 == 1) ? w.nps : i1 - 1
    @views for i = 1:w.nfnc
        w.c[:, i, i1] .=
            w.alpha *
            (w.expnegus[:, i] .* exp.(w.e[:, i, i0]) .- w.e[:, i, i0] .- dp(1.0)) .+
            (dp(1.0) - w.alpha) * w.c[:, i, i0]
    end
end

"""
    hnc_ng!(w::Wizard)

This method implements the Ng method K-C Ng, J. Chem. Phys. v61, 2680 (1974), doi:[10.1063/1.1682399](https://doi.org/10.1063/1.1682399)
as an accelerated solver for the HNC closure. 
As we iterate we move forward in the history trajectory, cyclically.
"""
function hnc_ng!(w::Wizard)
    w.istep += 1
    i1 = mod(w.istep - 1, w.nps) + 1
    # i0 = i1 - 1
    # i0 = (i0 == 0) ? nps : i0
    i0 = (i1 == 1) ? w.nps : i1 - 1
    nd = (w.istep <= w.nps) ? w.istep - 2 : w.nps - 1
    dc = zeros(w.ng - 1, w.nfnc, w.nps - 1)
    de = zeros(w.ng - 1, w.nfnc, w.nps - 1)

    @views for p = 1:nd
        j1 = i0 - p
        j1 = (j1 <= 0) ? (w.nps + j1) : j1
        dc[:, :, p] .= w.c[:, :, i0] .- w.c[:, :, j1]
        de[:, :, p] .= w.e[:, :, i0] .- w.e[:, :, j1]
    end
    a = zeros(w.nps - 1, w.nps - 1)
    x = zeros(w.nps - 1)
    y = zeros(w.nps - 1)

    @views for icp = 1:w.nfnc
        for j = 1:w.ng-1
            aux = w.expnegus[j, icp] * exp(w.e[j, icp, i0]) - dp(1.0)
            for j1 = 1:nd
                y[j1] = aux * de[j, icp, j1] - dc[j, icp, j1]
            end
            yy = aux - w.e[j, icp, i0] - w.c[j, icp, i0]
            for j1 = 1:nd
                for j2 = j1:nd
                    a[j1, j2] += y[j1] * y[j2]
                end
                x[j1] += y[j1] * yy
            end
        end
    end
    # ipiv = zeros(Int64, nps - 1)
    # info = 0
    # DSYSV("U", nd, 1, a, nps - 1, ipiv, x, nps - 1, work, lwork, info)
    # https://support.nag.com/numeric/fl/nagdoc_fl26/pdf/f07/f07maf.pdf
    try
        LinearAlgebra.LAPACK.sysv!('U', a, x)
    catch ex
        if isa(ex, SingularException)
            w.return_code = DSYSV_ERROR
            w.error_msg = "DSYSV encountered singular problem in hnc_ng"
            if !w.suppress_msgs && !w.silent
                println(stdout, "** error: ", w.error_msg)
                println(stdout, "** further messages of this kind will be suppressed")
                w.suppress_msgs = true
            end
            return nothing
        elseif isa(ex, DimensionMismatch)
            exit()
        else
        end
    end

    @views for icp = 1:w.nfnc
        for j = 1:w.ng-1
            aux = w.e[j, icp, i0]
            for j1 = 1:nd
                aux -= de[j, icp, j1] * x[j1]
            end
            w.c[j, icp, i1] = w.expnegus[j, icp] * exp(aux) - aux - dp(1.0)
        end
    end
end

"""
    hnc_solve!(w::Wizard)

Basic driver routine for solving HNC: take a number of Picard
iterations to pump-prime the Ng method.  Stop when error is less
than tolerance, or when exceed maximum number of iterations.  The
flag `cold_start` indicates whether the direct correlation function
should be re-initialised.  

The initial guess to the direct correlation function is either 
* `0` (zero) (`start_type = 1`), or 
* `c = - Ushort` (`start_type = 2`), or 
* `c = e^(-Ushort)-1` (`start_type = 3`).   

Any of these should do in principle, 
but the initial convergence may be different.
Note from above that `c` is actually defined `c' = c + Ulong`, 
ie with the long-range part of the potential added.  

Note that we always start from position `1` in the history trajectory, and
the history trajectory is pump-primed by Picard steps before
attempting the Ng accelerator.  
At the end, the final solution is copied back to position `1` in the history trajectory.
"""
function hnc_solve!(w::Wizard)
    w.return_code = NO_ERROR
    w.istep = 1
    if w.cold_start
        if w.start_type == 1
            w.c[:, :, 1] .= dp(0.0)
            w.verbose && println(stdout, "HNC cold start c' = 0")
        elseif w.start_type == 2
            w.c[:, :, 1] .= -w.ushort[:, :]
            w.verbose && println(stdout, "HNC cold start c\" = -v\"")
        elseif w.start_type == 3
            w.c[:, :, 1] .= w.expnegus[:, :] .- dp(1.0)
            w.verbose &&
                println(stdout, "HNC cold start c\" = e^(-v\")-1")
        else
            throw(ArgumentError("start_type should be either 1, 2 or 3."))
        end
        w.cold_start = false
    else
        (w.verbose) && println(stdout, "HNC warm start c\" = previous c\"")
    end

    oz_solve!(w)
    if (w.return_code > NO_ERROR)
        return
    end
    i = 1
    for outer i = 1:w.maxsteps
        if i <= w.npic
            hnc_picard!(w)
        else
            hnc_ng!(w)
        end
        oz_solve!(w)
        if (w.return_code > NO_ERROR)
            return
        end
        conv_test!(w)
        if w.verbose
            if i <= w.npic
                println(stdout, i, " HNC Picard, error = ", w.ERROR)
            else
                println(stdout, i, " HNC     Ng, error = ", w.ERROR)
            end
        end
        if (w.ERROR < w.tol)
            break
        end
    end
    w.nsteps = i
    if w.ERROR > w.tol
        w.return_code = CONVERGENCE_ERROR
        w.error_msg = "error > tol in hnc_solve"
        if (!w.silent)
            println(stdout, "** warning: ", w.error_msg)
        end
    end
    i1 = mod(w.istep - 1, w.nps) + 1
    if i1 != 1 # copy solution to position 1
        w.c[:, :, 1] .= w.c[:, :, i1]
        w.e[:, :, 1] .= w.e[:, :, i1]
    end
    w.closure_name = "HNC"
    w.closure_type = HNC_CLOSURE
    if w.auto_fns
        make_pair_functions!(w)
        make_structure_factors!(w)
        make_thermodynamics!(w)
    end
end

"""
    conv_test!(w::Wizard)

Calculate the difference between the direct correlation functions
for the current and previous iteration, used as a convergence test;
return answer in variable 'error'.
"""
function conv_test!(w::Wizard)
    i1 = mod(w.istep - 1, w.nps) + 1
    i0 = i1 - 1
    i0 = (i0 == 0) ? w.nps : i0
    w.ERROR = sqrt(w.deltar * sum((w.c[:, :, i1] .- w.c[:, :, i0]) .^ 2))
    #    norm = sqrt(deltar * sum( c(:, :, i1)**2 ))
    #    print *, "conv_test: norm = ", norm
end


"""
    msa_picard!(w::Wizard)

This method implements the MSA closure expressed as `c' = - e' - 1`
within the hard core, in similar terms to the HNC closure.  
Outwith the hard core, `c' = - beta  v'` is left untouched, 
presuming it is set
correctly in the MSA initialisation step.  As we iterate we move
forward in the history trajectory, cyclically.
"""
function msa_picard!(w::Wizard)
    w.istep = w.istep + 1
    i1 = mod(w.istep - 1, w.nps) + 1
    i0 = i1 - 1
    i0 = (i0 == 0) ? w.nps : i0
    @views for i = 1:w.nfnc
        irc = NINT(w.dd[i] / w.deltar) # Only work inside the hard core
        w.c[1:irc, i, i1] .=
            w.alpha * (-w.e[1:irc, i, i0] .- dp(1.0)) .+
            (dp(1.0) - w.alpha) * w.c[1:irc, i, i0]
    end
end

"""
    msa_ng!(w::Wizard)

This method implements the Ng method K-C Ng, J. Chem. Phys. v61, 2680 (1974), doi: [10.1063/1.1682399](https://doi.org/10.1063/1.1682399) as an accelerated solver for the MSA closure (cf HNC above).  As we iterate we move forward in the history
trajectory, cyclically.
"""
function msa_ng!(w::Wizard)
    w.istep += 1
    i1 = mod(w.istep - 1, w.nps) + 1
    i0 = i1 - 1
    i0 = (i0 == 0) ? w.nps : i0
    nd = (w.istep <= w.nps) ? w.istep - 2 : w.nps - 1
    dc = similar(w.c)
    de = similar(w.e)
    for p = 1:nd
        j1 = i0 - p
        j1 = (j1 <= 0) ? w.nps + j1 : j1
        dc[:, :, p] .= w.c[:, :, i0] .- w.c[:, :, j1]
        de[:, :, p] .= w.e[:, :, i0] .- w.e[:, :, j1]
    end
    a = zeros(w.nps - 1, w.nps - 1)
    x = zeros(w.nps - 1)
    y = zeros(w.nps - 1)
    # , y[nps-1], yy, aux

    for icp = 1:w.nfnc
        irc = NINT(w.dd[icp] / w.deltar) # Only work inside the hard core
        for j = 1:irc
            aux = -dp(1.0)
            for j1 = 1:nd
                y[j1] = aux * de[j, icp, j1] - dc[j, icp, j1]
            end
            yy = aux - w.e[j, icp, i0] - w.c[j, icp, i0]
            for j1 = 1:nd
                for j2 = j1:nd
                    a[j1, j2] = a[j1, j2] + y[j1] * y[j2]
                end
                x[j1] = x[j1] + y[j1] * yy
            end
        end
    end

    try
        LinearAlgebra.LAPACK.sysv!('U', a, x)
    catch ex
        if isa(ex, SingularException)
            w.return_code = DSYSV_ERROR
            w.error_msg = "DSYSV encountered singular problem in msa_ng"
            if !w.suppress_msgs && !w.silent
                println(stdout, "** error: ", w.error_msg)
                println(stdout, "** further messages of this kind will be suppressed")
                w.suppress_msgs = true
            end
            return nothing
        else
            rethrow(ex)
        end
    end

    #if info > 0
    #    w.return_code = DSYSV_ERROR
    #    error_msg = "DSYSV encountered singular problem in msa_ng"
    #    if !w.suppress_msgs && !w.silent
    #        println(stdout, "** error: ", error_msg)
    #        println(stdout, "** further messages of this kind will be suppressed")
    #        w.suppress_msgs = true
    #    end
    #    return nothing
    # end
    for icp = 1:w.nfnc
        irc = NINT(w.dd[icp] / w.deltar) # Only work inside the hard core
        for j = 1:irc
            aux = w.e[j, icp, i0]
            for j1 = 1:nd
                aux = aux - de[j, icp, j1] * x[j1]
            end
            w.c[j, icp, i1] = -aux - dp(1.0)
        end
    end
end

"""
    msa_solve!(w::Wizard)

Basic driver routine for solving MSA: take a number of Picard
iterations to pump-prime the Ng method.  Stop when error is less
than tolerance, or when exceed maximum number of iterations.  The
flag `cold_start` indicates whether the direct correlation function
should be re-initialised.  

For a cold start, the initial guess to
the direct correlation function inside the hard core is `c' = -1`.

Irrespective of this, 
we need to make sure `c'` is correctly initialised outwith the hard core.  

Note that we always start from
position `1` in the history trajectory, and the history trajectory is
pump-primed by Picard steps before attempting the Ng accelerator.
At the end, the final solution is copied back to position `1` in the history trajectory.
"""
function msa_solve!(w::Wizard)
    w.return_code = NO_ERROR
    @views for i = 1:w.nfnc
        irc = NINT(w.dd[i] / w.deltar) # Reset everywhere outwith hard core
        for p = 1:w.nps
            w.c[irc+1:w.ng-1, i, p] .= -w.ushort[irc+1:w.ng-1, i]
        end
    end
    w.istep = 1
    @views if w.cold_start
        for i = 1:w.nfnc
            irc = NINT(w.dd[i] / w.deltar) # Only initialise inside the hard core
            w.c[1:irc, i, 1] .= -dp(1.0)
        end
        w.cold_start = false
        w.verbose && println(stdout, "MSA cold start c' = -1")
    else
        w.verbose && println(stdout, "MSA warm start c\" = previous c\"")
    end
    oz_solve!(w)
    if (w.return_code > NO_ERROR)
        return
    end
    i = 1
    for outer i = 1:w.maxsteps
        if i <= w.npic
            msa_picard!(w)
        else
            msa_ng!(w)
        end
        oz_solve!(w)
        if (w.return_code > NO_ERROR)
            return
        end
        conv_test!(w)
        if w.verbose
            if i <= w.npic
                println(stdout, i, "MSA Picard, error = ", w.ERROR)
            else
                println(stdout, i, "MSA     Ng, error = ", w.ERROR)
            end
        end
        if (w.ERROR < w.tol)
            break
        end
    end
    w.nsteps = i
    if w.ERROR > w.tol
        w.return_code = CONVERGENCE_ERROR
        w.error_msg = "error > tol in msa_solve"
        (!w.silent) && println(stdout, "** warning: ", w.error_msg)
    end
    i1 = mod(w.istep - 1, w.nps) + 1
    if i1 != 1 # copy solution to position 1
        @views begin
            w.c[:, :, 1] .= w.c[:, :, i1]
            w.e[:, :, 1] .= w.e[:, :, i1]
        end
    end
    w.closure_name = "MSA"
    w.closure_type = MSA_CLOSURE
    if w.auto_fns
        make_pair_functions!(w)
        make_structure_factors!(w)
        make_thermodynamics!(w)
    end
end

"""
    rpa_solve!(w::Wizard)

Given the HNC machinery, the implementation of the RPA is almost
completely trivial and corresponds to one iteration through the
Ornstein-Zernike solver given the choice c = - Ushort.  We save the
result to position 1 in the history trajectory.
"""
function rpa_solve!(w::Wizard)

    w.return_code = NO_ERROR
    w.istep = 1
    @views w.c[:, :, 1] .= -w.ushort[:, :]
    oz_solve!(w)
    if (w.return_code > NO_ERROR)
        return
    end
    w.ERROR = dp(0.0)
    w.closure_name = "RPA"
    w.closure_type = RPA_CLOSURE
    if w.auto_fns
        make_pair_functions!(w)
        make_structure_factors!(w)
        make_thermodynamics!(w)
    end
end

"""
    save_reference(w::Wizard)

Save the reference state, assuming the `c` and `e` functions 
are those in the position 1 in the history trajectory.  
The corresponding `c` and `e` functions can be restored by a call to `oz_solve2`.
"""
function save_reference(w::Wizard)
    @views w.h0[:, :] .= w.c[:, :, 1] + w.e[:, :, 1]
end

"""
    exp_refine!(w::Wizard)

The EXP approximation refines the current RPA/MSA solution by using
the current solution (`h = c + e`) and a reference solution (`h0`) to
send `h0 --> (1 + h0) exp(h - h0) - 1`.  A round trip through the
alternate version of the Ornstein-Zernike relation re-calculates the
direct and indirect correlation functions.  We assume the current
solution is in position 1 in the history trajectory.  The reference
state is reset to `h0 = 0` after a call to this function, for safety!
"""
function exp_refine!(w::Wizard)

    hsave = zeros(w.ng - 1, w.nfnc)
    hsave .= w.h0
    w.h0 .=
        (dp(1.0) .+ w.h0) .*
        exp.(w.c[:, :, 1] + w.e[:, :, 1] - w.h0) .- dp(1.0)
    oz_solve2!(w)
    if (w.return_code > NO_ERROR)
        return
    end
    w.h0 = hsave
    w.closure_name = "EXP"
    w.closure_type = EXP_CLOSURE
    if w.auto_fns
        make_pair_functions!(w)
        make_structure_factors!(w)
        make_thermodynamics!(w)
    end
end

"""
    make_structure_factors!(w::Wizard)

Construct the structure factors out of the transform of the total correlation function.  
Note that `ck` and `ek` are available after a call to the OZ solver.

# Fields calculated by this method
* `hk` : fourier transformed pair functions
* `sk` : fourier transformed partial structure factors
"""
function make_structure_factors!(w::Wizard)

    # hk = zeros(ng - 1, w.nfnc)
    hk = similar(w.ck)
    hk .= w.ck + w.ek
    @views for (i, j, ij) in upperTriangularIndices(w)
        if i == j
            w.sk[:, i, i] .= w.rho[i] .* (dp(1.0) .+ w.rho[i] * hk[:, ij])
        else
            w.sk[:, i, j] .= w.rho[i] * w.rho[j] * hk[:, ij]
            w.sk[:, j, i] .= w.sk[:, i, j]
        end
    end
end

"""
    make_pair_functions!(w::Wizard)

Construct the total correlation functions out of the direct
correlation functions.
Note that the above routines actually works
with `c' = c + Ulong` and `e' = e - Ulong` 
where `Ulong` is the long-range part of the potential, 
but `h = g - 1 = e + c = e' + c'`.  
The pair correlation functions are `g = 1 + h` - the addition of `1` is left
for the user to implement.  
We assume the `c` and `e` functions are those in the position `1` in the history trajectory.

# Fields calculated by this method
* `hc` : current contact values
* `hr` : current total correlation functions
"""
function make_pair_functions!(w::Wizard)

    @views for (i, j, ij) in upperTriangularIndices(w)
        w.hr[:, i, j] .= w.c[:, ij, 1] + w.e[:, ij, 1]
        w.hr[:, j, i] .= w.hr[:, i, j]
    end

    # The contact value in the case of hard cores.  We extrapolate
    # from the two nearest points.  The same extrapolation is done in
    # make_thermodynamics but we want to keep these subroutines
    # independent.

    @views for (i, j, ij) in upperTriangularIndices(w)
        w.hc[i, j] = if ispositive(w.dd[ij])
            let irc = NINT(w.dd[ij] / w.deltar),
                r1 = w.r[irc+1],
                r2 = w.r[irc+2],
                h1 = w.c[irc+1, ij, 1] + w.e[irc+1, ij, 1],
                h2 = w.c[irc+2, ij, 1] + w.e[irc+2, ij, 1]
                # interpolated value
                ((h1 - h2) * w.dd[ij] + h2 * r1 - h1 * r2) / (r1 - r2)
            end
        else
            dp(0.0)
        end
        if i < j
            w.hc[j, i] = w.hc[i, j]
        end
    end
end

"""
    make_thermodynamics!(w::Wizard)

Calculate various thermodynamics properties by spatial integration
(as contrasted to thermodynamic integration).  We use the trapezium
rule, taking account where necessary the end-point values (see
intro), `h_0 = 2 h_1 - h_2` at `r = 0` (`i = 1`), 
and `h_n = 0` at `r = L` (`i = ng`).  
Also we have `r = i * deltar` for `i = 1` to `ng-1`.  
Note that the above routines actually 
work with `c' = c + Ulong` and  `e' = e - Ulong`
where `Ulong` is the long-range part of the potential, 
so we have `h = g - 1 = e + c = e' + c'`.  
See also Vrbka et al, J. Chem. Phys. 131, 154109 (2009), doi:[10.1063/1.3248218](https://doi.org/10.1063/1.3248218).
We assume the `c` and `e` functions are those in the position `1` in the history trajectory.

The mean-field thermodynamic expressions can often be obtained
analytically from the potential. In this routine they are calculated
from species pair contributions, which are themselves calculated in
the potential routines (which do not have access to the densities).

# Fields calculated by this method
* `cf_mf` : Compressibility factor: mean field contribution
* `cf_xc` : Compressibility factor: contact contribution
* `cf_gc` : Compressibility factor: correlation contribution
* `pex` :  Excess pressure, virial route
* `press` : Pressure, virial route
* `comp` : Compressibility: total
* `comp_xc` : Compressibility: correlation contribution
* `uex` : Internal energy density
* `uex_mf` : Internal energy density: mean field contribution
* `uex_xc` : Internal energy density: correlation contribution
* `aex_rl` : HNC excess free energy density, real space (long)
* `aex_rs` : HNC excess free energy density, real space (short)
* `aex_kd` : HNC excess free energy density, k-space log(det)
* `aex_kl` : HNC excess free energy density, k-space tr (long)
* `aex_ks` : HNC excess free energy density, k-space tr (short)
* `aex` : HNC total free energy density
* `deficit` : HNC excess free energy deficit (see docs)
"""
function make_thermodynamics!(w::Wizard)

    rhotot = sum(w.rho)

    # rhoxx is rho x_i x_j, doubled up for the off-diagonal components
    rhoxx = zeros(w.nfnc)
    @views for (i, j, ij) in upperTriangularIndices(w)
        rhoxx[ij] = if i == j
            w.rho[i]^2 / rhotot
        else
            dp(2.0) * w.rho[i] * w.rho[j] / rhotot
        end
    end

    # Calculate the various contributions to the virial-route
    # pressure.  This is the mean field contribution.

    w.cf_mf = sum(rhoxx .* w.tp)

    # Evaluate t_ij = - (2pi/3) int_d^inf d(U_ij)/dr h_ij r^3 dr. The
    # contribution from both end-points r = 0 (i = 0) and r = L (i =
    # ng) vanishes, hence the sum just consists of the middle part of
    # the trapezium rule.  With a hard core the first point
    # contributes half.

    t = zeros(w.nfnc)
    @views for i = 1:w.nfnc
        let integrand = J ->
                (w.dushort[J, i] + w.dulong[J, i]) *
                (w.c[J, i, 1] + w.e[J, i, 1]) * w.r[J]^3
            t[i] = if ispositive(w.dd[i])
                let irc = NINT(w.dd[i] / w.deltar)
                    # trapezium rule
                    integrand(irc + 1) / 2 +
                    sum(integrand(J) for J in irc+2:length(w.r))
                end
            else
                sum(integrand(J) for J in eachindex(w.r))
            end
        end
    end
    t .*= w.deltar * (-2π / 3)

    # The correlation contribution is sum_ij rho x_i x_j t_ij.
    w.cf_xc = sum(rhoxx .* t)

    # The contact contribution in the case of hard cores.  We
    # extrapolate the contact value of the pair distribution function
    # from the two nearest outside points.  The same extrapolation is
    # done in make_pair_functions but we want to keep these
    # subroutines independent.

    @views for i = 1:w.nfnc
        t[i] = if ispositive(w.dd[i])
            let irc = NINT(w.dd[i] / w.deltar),
                r1 = w.r[irc+1],
                r2 = w.r[irc+2],
                h1 = w.c[irc+1, i, 1] + w.e[irc+1, i, 1],
                h2 = w.c[irc+2, i, 1] + w.e[irc+2, i, 1],
                hc = ((h1 - h2) * w.dd[i] + h2 * r1 - h1 * r2) / (r1 - r2)
                # value
                w.dd[i]^3 * (dp(1.0) + hc) * (2π / 3)
            end
        else
            dp(0.0)
        end
    end
    w.cf_gc = sum(rhoxx .* t)

    # This is the final pressure.

    w.pex = rhotot * (w.cf_gc + w.cf_mf + w.cf_xc)
    w.press = rhotot + w.pex

    # Now we do the compressibility (not to be confused with the above
    # compressibility factor).  The long range part is accounted for
    # separately.

    # Evaluate t_ij = 4 pi int_0^inf c_ij r^2 dr.  Again the
    # contribution from both endpoints vanishes, hence the sum just
    # consists of the middle part of the trapezium rule.

    @views for i = 1:w.nfnc
        t[i] = sum(w.c[:, i, 1] .* w.r .^ 2)
    end
    t .*= (4π) * w.deltar

    # The compressibility is 1 - sum_ij rho x_i x_j t_ij

    w.comp_xc = -sum(rhoxx .* (t .- w.tl))
    w.comp = dp(1.0) + w.comp_xc

    # Now we do the energy density.  First the mean-field
    # contribution.

    w.uex_mf = rhotot * sum(rhoxx .* w.tu)

    # Evaluate t_ij = 2 pi int_0^inf U_ij h_ij r^2 dr.  Note that with
    # a hard core the contribution from the inner end-point is halved
    # (trapezium rule) and the outer vanishes, otherwise the
    # contribution from both ends vanishes.

    @views for i = 1:w.nfnc
        let integrand = J ->
                (w.ushort[J, i] + w.ulong[J, i]) *
                (w.c[J, i, 1] + w.e[J, i, 1]) * w.r[J]^2
            t[i] = if ispositive(w.dd[i])
                let irc = NINT(w.dd[i] / w.deltar)
                    integrand(irc + 1) / 2 +
                    sum(integrand(J) for J in irc+2:length(w.r))
                end
            else
                sum(integrand(J) for J in eachindex(w.r))
            end
        end
    end
    t .*= w.deltar * (2π)

    # The extra contribution is sum_ij rho x_i x_j t_ij
    w.uex_xc = rhotot * sum(rhoxx .* t)

    # This is the final energy density.
    w.uex = w.uex_mf + w.uex_xc

    # Finally do the free energy and chemical potentials for HNC

    if w.closure_type == HNC_CLOSURE

        # Evaluate `t_ij = 2 pi int_0^inf ( h_ij^2 / 2 - c_ij) r^2 dr` where
        # we use c = c' - `Ulong`  in the second term and pull this out as a
        # long range contribution.

        @views for i = 1:w.nfnc
            t[i] = sum(
                (dp(0.5) * (w.c[:, i, 1] .+ w.e[:, i, 1]) .^ 2 .-
                 w.c[:, i, 1]) .* w.r .^ 2)
        end
        t .*= (4π) * w.deltar
        w.aex_rs = sum(rhoxx .* t) * rhotot / 2
        w.aex_rl = sum(rhoxx .* w.tl) * rhotot / 2

        # Call out to compute the reciprocal space contributions

        hnc_kspace_integrals!(w)

        # Final total excess free energy density
        w.aex = w.aex_rs + w.aex_rl + w.aex_kd + w.aex_ks + w.aex_kl

        # Evaluate `t_ij = 4 pi int_0^inf (h_ij e_ij / 2 - c_ij) r^2 dr`.
        # Note that the contribution from both end-points again vanishes.
        # We can use `c'` for `c` in the second term because the long range
        # part is accounted for analytically (it is the same integral as
        # appears in the compressibility), 
        # but we must use `e = e' + Ulong  for the first term.

        @views for i = 1:w.nfnc
            let integrand = J ->
                    (dp(0.5) * (w.c[J, i, 1] + w.e[J, i, 1]) *
                     (w.e[J, i, 1] + w.ulong[J, i]) -
                     w.c[J, i, 1]) * w.r[J]^2
                t[i] = sum(integrand(J) for J in eachindex(w.r))
            end
        end
        t .*= (4π) * w.deltar

        # The excess chemical potential of the `i`-th component is then `sum_j rho_j t_ij`

        w.muex .= dp(0.0)
        @views for i = 1:w.ncomp
            for j = 1:w.ncomp
                if i <= j
                    ij = i + j * (j - 1) ÷ 2
                else
                    ij = j + i * (i - 1) ÷ 2
                end
                w.muex[i] += w.rho[j] * (t[ij] + w.tl[ij])
            end
        end

        # The free energy density should satisfy
        # `f - sum_mu rho_mu mu_mu + p = 0`.
        # The deficit/excess here is a test of the numerics

        w.deficit = w.aex - sum(w.rho .* w.muex) + w.pex
    end
end

"""
    hnc_kspace_integrals!(w::Wizard)

Calculate the reciprocal space contributions to the HNC free energy.  

Note that `c = c' - Ulong`, and `ck` is available after a call to OZ solver.
"""
function hnc_kspace_integrals!(w::Wizard)

    # Calculate the short range trace fn in k-space

    TR = zeros(size(w.ck, 1))
    @views for (i, ij) in diagonalUpperTriangularIndices(w)
        TR .+= w.rho[i] * w.ck[:, ij]
    end

    # Calculate log(det) fn in k-space
    lndet = zeros(w.ng - 1)

    if w.ncomp == 1 # One-component case

        @views lndet[:] .= log.(dp(1.0) .- w.rho[1] * (w.ck[:, 1] .- w.ulongk[:, 1]))

    else # Multicomponent case

        # First set up a unit matrix, and the diagonal R matrix

        rhomat = Diagonal(w.rho)

        # Do the matrix calculations for each wavevector k.
        csubumat_parent = zeros(w.ncomp, w.ncomp)
        csubumat = Symmetric(csubumat_parent) # construct U from LU part
        for ik = 1:w.ng-1

            # Unpack the reciprocal space functions into matrices.
            # Calculate C - beta UL first.

            @views for (i, j, ij) in upperTriangularIndices(w)
                csubumat_parent[i, j] = w.ck[ik, ij] - w.ulongk[ik, ij]
            end

            # Construct M = I - R . (C - beta UL)
            m = LinearAlgebra.I - (rhomat * csubumat)

            # Find the eigenvalues of M (which is _not_ symmetric)
            # DGEEV("N", "N", ncomp, m, ncomp, wr, wi, m, ncomp, m, ncomp, work, lwork, info)
            # 
            local ww
            try
                ww = eigvals(m)
            catch ex
                # info = 0
                # if info > 0
                if isa(ex, SingularException)
                    w.return_code = DGEEV_ERROR
                    w.error_msg = "DGEEV failed to converge in hnc_kspace_integrals"
                    if !w.silent
                        println(stdout, "** error: ", w.error_msg)
                    end
                    return nothing
                else
                    rethrow(ex)
                end
            end

            # Calculate ln det M = 1/2 sum_i ln |lambda_i|^2
            lndet[ik] = sum(log ∘ abs2, ww) / 2

        end # loop over k vectors

    end # select single component or multicomponent case

    # Finally the contribution is the integral of ln(det) + trace.
    # The long range contribution in the trace is taken care of
    # analytically as it resolves to the r = 0 value of the long range
    # potential (see documentation).

    prefac = dp(1.0) / ((4π) * (π))

    w.aex_kd = prefac * w.deltak * sum(lndet .* w.k .^ 2)
    w.aex_ks = prefac * w.deltak * sum(TR .* w.k .^ 2)
    w.aex_kl = sum(w.rho .* w.u0) / (-2)
end

"""
    write_thermodynamics(w::Wizard)

Write out parameters for thermodynamics.
"""
function write_thermodynamics(w::Wizard)
    #    integer :: i, j, ij

    if w.closure_type == NO_CLOSURE
        println(stdout, "No closure = no thermodynamics")
        return
    end

    if w.model_type == DPD_LINEAR_CHARGES
        println(stdout, "No thermodynamics for ", w.model_name)
        return
    end

    println(stdout, "Thermodynamics for ", w.model_name)
    println(stdout, w.closure_name, " closure, convergence = ", w.ERROR)
    println(stdout, "Compressibility factor: mean field contribution = ", w.cf_mf)
    if any(ispositive.(w.dd))
        println(stdout, "Compressibility factor: contact contribution = ", w.cf_gc)
    end
    println(stdout, "Compressibility factor: correlation contribution = ", w.cf_xc)
    println(stdout, "Compressibility factor: total = ",
        dp(1.0) + w.cf_mf + w.cf_gc + w.cf_xc)
    println(stdout, "Excess pressure (virial route) = ", w.pex)
    println(stdout, "Pressure (virial route) = ", w.press)
    println(stdout, "Compressibility: correlation contribution = ", w.comp_xc)
    println(stdout, "Compressibility: total = ", w.comp)
    println(stdout, "Internal energy density: mean field contribution = ", w.uex_mf)
    println(stdout, "Internal energy density: correlation contribution = ", w.uex_xc)
    println(stdout, "Internal energy density = ", w.uex)
    sum_rho = sum(w.rho)
    println(stdout, "Internal energy per particle = ", w.uex / sum_rho)
    println(stdout, "Internal energy per particle / 3 = ", w.uex / (dp(3.0) * sum_rho))

    @views for (i, j, ij) in upperTriangularIndices(w)
        if ispositive(w.dd[ij])
            println(stdout, "Contact g", i, j, " = 1 + ", w.hc[i, j], " = ", 1 + w.hc[i, j])
        end
    end

    if w.closure_type == HNC_CLOSURE
        for i = 1:w.ncomp
            println(stdout, "HNC: species", i, " muex = ", w.muex[i])
        end
        println(stdout, "HNC: real space (short),  aex_rs = ", w.aex_rs)
        println(stdout, "HNC: real space (long),   aex_rl = ", w.aex_rl)
        println(stdout, "HNC: k-space, log(det),   aex_kd = ", w.aex_kd)
        println(stdout, "HNC: k-space, tr (short), aex_ks = ", w.aex_ks)
        println(stdout, "HNC: k-space, tr (long),  aex_kl = ", w.aex_kl)
        println(stdout, "HNC: total free energy density = ", w.aex)

        println(stdout, "HNC:   --ditto--  per particle = ", w.aex / sum_rho)
        println(stdout, "HNC: sum rho muex - pex        = ", sum(w.rho .* w.muex) - w.pex)
        println(stdout, "HNC: deficit = ", w.deficit)
    end
end



# This is a Julia port of SunlightHNC 
# Copyright (C) 2025 Hiroharu Sugawara, Tokyo Metropolitan University <hsugawa@tmu.ac.jp>

"""
    Series_hr

Plots recipe to generate data series of total correlation function `hr` vs `r`.

# Example
```julia
w = Wizard()
using Plots
plot(Series_hr(), w; rcut=30)
```
# Keyword arguments
* `rcut=Inf` : Restrict `r < rcut`. No restriction for default.
* `nlimit`: Plot only the data of the component <= `nlimit <= w.ncomp`
"""
struct Series_hr end
@recipe function f(::Series_hr, w::Wizard; rcut=Inf, nlimit=100)
    the_nlimit = nlimit
    delete!(plotattributes, :nlimit)
    the_nlimit = min(the_nlimit, w.ncomp)
    IR = w.r .< rcut
    delete!(plotattributes, :rcut)
    title --> "Total correlation function"
    xlabel --> L"$r$"
    for (i, j, ij) in upperTriangularIndices(w, the_nlimit)
        @series begin
            label --> "hr$(i)$(j)"
            w.r[IR], w.hr[IR, i, j]
        end
    end
end


"""
    Series_gr

Plots recipe to generate data series of pair distribution function `gr = 1 + hr` vs `r`.

# Example
```julia
w = Wizard()
using Plots
plot(Series_gr(), w; rcut=30)
```
# Keyword arguments
* `rcut=Inf` : Restrict `r < rcut`. No restriction for default.
* `nlimit`: Plot only the data of the component <= `nlimit <= w.ncomp`
"""
struct Series_gr end
@recipe function f(::Series_gr, w::Wizard; rcut=Inf, nlimit=100)
    the_nlimit = nlimit
    delete!(plotattributes, :nlimit)
    the_nlimit = min(the_nlimit, w.ncomp)
    IR = w.r .< rcut
    delete!(plotattributes, :rcut)
    title --> "Pair distribution function"
    xlabel --> L"$r$"
    for (i, j, ij) in upperTriangularIndices(w, the_nlimit)
        @series begin
            label --> "gr$(i)$(j)"
            w.r[IR], 1e0 .+ w.hr[IR, i, j]
        end
    end
end

"""
    Series_cr

Plots recipe to generate data series of `c` vs `r`.

# Example
```julia
using Plots
w = Wizard()
plot(Series_cr(), w; rcut=30)
```
# Keyword arguments
* `rcut=Inf` : Restrict `r < rcut`. No restriction for default.
* `nlimit`: Plot only the data of the component <= `nlimit <= w.ncomp`
"""
struct Series_cr end
@recipe function f(::Series_cr, w::Wizard; rcut=Inf, nlimit=100)
    the_nlimit = nlimit
    delete!(plotattributes, :nlimit)
    the_nlimit = min(the_nlimit, w.ncomp)
    IR = w.r .< rcut
    delete!(plotattributes, :rcut)
    xlabel --> L"$r$"
    for (i, j, ij) in upperTriangularIndices(w, the_nlimit)
        @series begin
            label --> "c$(i)$(j)"
            w.r[IR], w.c[IR, i, j]
        end
    end
end

"""
    Series_vs_r

Plots recipe to generate data series vs `r`.

# Example
```julia
w = Wizard()
v = similar(w.r)
using Plots
plot(Series_vs_r(), w, v; rcut=30)
```
# Keyword arguments
* `rcut=Inf` : Restrict `r < rcut`. No restriction for default.
* `nlimit`: Plot only the data of the component <= `nlimit <= w.ncomp`
"""
struct Series_vs_r end
@recipe function f(::Series_vs_r, w::Wizard, v::AbstractVector; rcut=Inf, nlimit=100)
    @assert size(w.r) == size(v)
    the_nlimit = nlimit
    delete!(plotattributes, :nlimit)
    the_nlimit = min(the_nlimit, w.ncomp)
    IR = w.r .< rcut
    delete!(plotattributes, :rcut)
    xlabel --> L"$r$"
    @series begin
        w.k[IR], v[IR]
    end
end

"""
   calc_dcr

Calculate `(c11+c12)/2` and `(c11-c12)/2`
"""
function calc_dcr(w::Wizard)
    c12ave = (w.c[:, 1, 1] .- w.ulong[:, 1] .+ w.c[:, 2, 1] .- w.ulong[:, 2]) / 2
    c12diff = (w.c[:, 1, 1] .- w.ulong[:, 1] .- w.c[:, 2, 1] .+ w.ulong[:, 2]) / 2
    return c12ave, c12diff
end


"""
    Series_rhr_tail

Plots recipe to generate data series of `log_{10}(r*hr)` vs `r`.

# Example
```julia
w = Wizard()
using Plots
plot(Series_hr(), w; rcut=30)
```
# Keyword arguments
* `eps=1e-20`
* `rcut=Inf` : Restrict `r < rcut`. No restriction for default.
* `nlimit`: Plot only the data of the component <= `nlimit <= w.ncomp`
"""
struct Series_rhr_tail end
@recipe function f(::Series_rhr_tail, w::Wizard; rcut=Inf, nlimit=100)
    the_nlimit = nlimit
    delete!(plotattributes, :nlimit)
    the_nlimit = min(the_nlimit, w.ncomp)
    IR = w.r .< rcut
    delete!(plotattributes, :rcut)
    the_eps = eps
    delete!(plotattributes, :eps)
    xlabel --> L"$r$"
    for (i, j, ij) in upperTriangularIndices(w, the_nlimit)
        @series begin
            label --> "r*hr$(i)$(j)"
            w.r[IR], log10.(the_eps .+ abs.(w.r[IR] .* w.hr[IR, i, j]))
        end
    end
end


"""
    Series_sk

Plots recipe to generate data series of structure factor `sk(k)` vs `k`.

# Example
```julia
w = Wizard()
using Plots
plot(Series_sk(), w; kcut=30)
```
# Keyword arguments
* `kcut=Inf` : Restrict `k < kcut`. No restriction for default.
* `nlimit`: Plot only the data of the component <= `nlimit <= w.ncomp`
"""
struct Series_sk end
@recipe function f(::Series_sk, w::Wizard; kcut=Inf, nlimit=100)
    the_nlimit = nlimit
    delete!(plotattributes, :nlimit)
    the_nlimit = min(the_nlimit, w.ncomp)
    IK = w.k .< kcut
    delete!(plotattributes, :kcut)
    title --> "Structure factor"
    xlabel --> L"$k$"
    for (i, j, ij) in upperTriangularIndices(w, the_nlimit)
        @series begin
            label --> "sk$(i)$(j)"
            w.k[IK], w.sk[IK, i, j]
        end
    end
end

"""
    Series_vs_k

Plots recipe to generate data series vs `k`.

# Example
```julia
w = Wizard()
v = similar(w.k)
using Plots
plot(Series_vs_k(), w, v; kcut=30)
```
# Keyword arguments
* `kcut=Inf` : Restrict `k < kcut`. No restriction for default.
* `nlimit`: Plot only the data of the component <= `nlimit <= w.ncomp`
"""
struct Series_vs_k end
@recipe function f(::Series_vs_k, w::Wizard, v::AbstractVector; kcut=Inf, nlimit=100)
    @assert size(w.k) == size(v)
    the_nlimit = nlimit
    delete!(plotattributes, :nlimit)
    the_nlimit = min(the_nlimit, w.ncomp)
    IK = w.k .< kcut
    delete!(plotattributes, :kcut)
    xlabel --> L"$k$"
    @series begin
        w.k[IK], v[IK]
    end
end


"""
    Series_vs_ksq

Plots recipe to generate data series vs `k^2`.

# Example
```julia
w = Wizard()
v = similar(w.k)
using Plots
plot(Series_vs_ksq(), w, v; ksqcut=30)
```
# Keyword arguments
* `ksqcut=Inf` : Restrict `k^2 < ksqcut`. No restriction for default.
* `nlimit`: Plot only the data of the component <= `nlimit <= w.ncomp`
"""
struct Series_vs_ksq end
@recipe function f(::Series_vs_ksq, w::Wizard, v::AbstractVector; kscut=Inf, nlimit=100)
    @assert size(w.k) == size(v)
    the_nlimit = nlimit
    delete!(plotattributes, :nlimit)
    the_nlimit = min(the_nlimit, w.ncomp)
    IK = w.k .^ 2 .< ksqcut
    delete!(plotattributes, :ksqcut)
    xlabel --> L"$k^2$"
    @series begin
        w.k[IK] .^ 2, v[IK]
    end
end


"""
    Series_hk

Plots recipe to generate data series of Fourier transformed total correlation function `hk = ek + ck` vs `k`.

# Example
```julia
w = Wizard()
using Plots
plot(Series_hk(), w; kcut=30)
```
# Keyword arguments
* `kcut=Inf` : Restrict `k < kcut`. No restriction for default.
* `nlimit`: Plot only the data of the component <= `nlimit <= w.ncomp`
"""
struct Series_hk end
@recipe function f(::Series_hk, w::Wizard; kcut=Inf, nlimit=100)
    the_nlimit = nlimit
    delete!(plotattributes, :nlimit)
    the_nlimit = min(the_nlimit, w.ncomp)
    IK = w.k .< kcut
    delete!(plotattributes, :kcut)
    title --> "Fourier transformed total correlation function"
    xlabel --> L"$k$"
    for (i, j, ij) in upperTriangularIndices(w, the_nlimit)
        @series begin
            label --> "hk$(i)$(j)"
            w.k[IK], w.ek[IK, ij] .+ w.ck[IK, ij]
        end
    end
end


"""
    calc_snn

Alias to [`calc_ddsf`](@ref)
"""
calc_snn(W::Wizard) = calc_ddsf(w)

"""
    calc_ddsf(w::Wizard)

Calculate density-density structure factor (ddsf)
"""
function calc_ddsf(w::Wizard)
    # snn = sum(sum(w.sk, dims=3), dims=2) / sum(w.rho)
    ddsf = sum(w.sk, dims=[2, 3]) / sum(w.rho)
    @assert size(ddsf) == (w.ng - 1, 1, 1)
    return view(ddsf, :, 1, 1)
end

"""
    calc_ccsf(w::Wizard)

Calculate charge-charge structure factor (ccsf)
"""
function calc_ccsf(w::Wizard)
    @assert w.comp >= 2
    ccsf = zeros(w.ng - 1)
    ccsf .= w.sk[:, 1, 1] .- w.sk[:, 1, 2] .- w.sk[:, 2, 1] .+ w.sk[:, 2, 2]
    return ccsf
end


"""
    calc_szz(w::Wizard)

Calculate `szz`
"""
function calc_szz(w::Wizard)
    # szz = dot(dot(w.z, w.sk), w.z) / sum(w.rho .* w.z .^ 2)
    szz = similar(w.k)
    for i in eachindex(w.k)
        szz[i] = dot(w.z, w.sk[i, :, :], w.z) / sum(w.rho .* w.z .^ 2)
    end
    # @show size(w.z), size(w.sk), size(snn), size(szz)
    return szz
end


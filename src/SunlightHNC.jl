
# This is a Julia port of SunlightHNC 
# Copyright (C) 2025 Hiroharu Sugawara, Tokyo Metropolitan University <hsugawa@tmu.ac.jp>

module SunlightHNC

export dp, NINT, ispositive
export NO_CLOSURE,
    HNC_CLOSURE,
    RPA_CLOSURE,
    MSA_CLOSURE,
    EXP_CLOSURE
export DPD_GAUSSIAN_CHARGES,
    DPD_BESSEL_CHARGES,
    DPD_LINEAR_CHARGES,
    DPD_SLATER_APPROX_CHARGES,
    DPD_SLATER_EXACT_CHARGES,
    URPM_WITHOUT_USHORT,
    URPM_WITH_USHORT,
    RPM_WITHOUT_USHORT,
    RPM_WITH_USHORT,
    HARD_SPHERES
export NO_ERROR,
    CONVERGENCE_ERROR,
    AXEQB_ERROR,
    DSYSV_ERROR,
    DGEEV_ERROR,
    MISMATCH_ERROR
export NO_CLOSURE,
    HNC_CLOSURE,
    RPA_CLOSURE,
    MSA_CLOSURE,
    EXP_CLOSURE

export Wizard,
    upperTriangularIndices,
    diagonalUpperTriangularIndices,
    write_params,
    dpd_potential!,
    urpm_potential!,
    rpm_potential!,
    hs_potential!,
    hnc_solve!,
    msa_solve!,
    save_reference,
    exp_refine!,
    write_thermodynamics

export Model,
    additive_primitive_model,
    restricted_primitive_model,
    solve!

export
    Series_hr,
    Series_gr,
    Series_sk, 
    Series_vs_k, 
    Series_vs_ksq,
    calc_ccsf,
    calc_ddsf,
    calc_szz

using FFTW
using LinearAlgebra
using SpecialFunctions
using StaticArrays
include("./wizard.jl")

using StructArrays
using StaticArrays
include("./model.jl")

using Plots
using RecipesBase
using LaTeXStrings
include("./plots.jl")

end # SunlightHNC


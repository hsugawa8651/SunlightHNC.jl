```@meta
CurrentModule = SunlightHNC
```

## Core

```@docs
SunlightHNC.Wizard
SunlightHNC.Wizard()
```


## Potentials

```@docs
SunlightHNC.dpd_potential!
SunlightHNC.rpm_potential!
SunlightHNC.urpm_potential!
SunlightHNC.hs_potential!
```

## [Constants for `charge_type`](@id charge_types)

```@docs
SunlightHNC.NO_MODEL_TYPE
SunlightHNC.DPD_GAUSSIAN_CHARGES
SunlightHNC.DPD_BESSEL_CHARGES
SunlightHNC.DPD_LINEAR_CHARGES
SunlightHNC.DPD_SLATER_APPROX_CHARGES
SunlightHNC.DPD_SLATER_EXACT_CHARGES
SunlightHNC.URPM_WITHOUT_USHORT
SunlightHNC.URPM_WITH_USHORT
SunlightHNC.RPM_WITHOUT_USHORT
SunlightHNC.RPM_WITH_USHORT
SunlightHNC.HARD_SPHERES
```

## Solvers

```@docs
SunlightHNC.oz_solve!
SunlightHNC.oz_solve2!
SunlightHNC.hnc_ng!
SunlightHNC.hnc_picard!
SunlightHNC.hnc_solve!
SunlightHNC.msa_ng!
SunlightHNC.msa_picard!
SunlightHNC.msa_solve!
SunlightHNC.rpa_solve!
```


## [Constants for `return_code`](@id return_codes)

```@docs
SunlightHNC.NO_ERROR
SunlightHNC.CONVERGENCE_ERROR
SunlightHNC.AXEQB_ERROR
SunlightHNC.DSYSV_ERROR
SunlightHNC.DGEEV_ERROR
SunlightHNC.MISMATCH_ERROR
```

## [Constants for `closure_type`](@id closure_types)

```@docs
SunlightHNC.NO_CLOSURE
SunlightHNC.HNC_CLOSURE
SunlightHNC.RPA_CLOSURE
SunlightHNC.MSA_CLOSURE
SunlightHNC.EXP_CLOSURE
```

## Solver utilities
```@docs
SunlightHNC.conv_test!
SunlightHNC.hnc_kspace_integrals!
SunlightHNC.exp_refine!
```

## Statistical Mechanics

```@docs
SunlightHNC.make_pair_functions!
SunlightHNC.make_structure_factors!
SunlightHNC.make_thermodynamics!
```

## Reports

```@docs
SunlightHNC.save_reference
SunlightHNC.write_params
SunlightHNC.write_thermodynamics
```

## Utilities

```@docs
SunlightHNC.upperTriangularIndices
SunlightHNC.diagonalUpperTriangularIndices
```

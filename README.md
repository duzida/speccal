# speccal

**Calibrated inference for specification curve (multiverse) analysis.**

Existing tools (`specr`, `multiverse`) enumerate specifications, fit them and draw the curve. `speccal`
adds what comes after the curve and what complex-survey data need before it:

- **Survey-design grid fitting** (`spec_fit`): `survey::svyglm` with subsets as domains, marginal risk
  differences and risk ratios by standardisation, and convergence / quasi-separation diagnostics for every fit.
- **Permutation nulls that preserve confounding** (`spec_null`): within-stratum Freedman–Lane permutation,
  saved specification by specification, alongside the simple permutation for comparison.
- **A null-based criterion for defensible specifications** (`spec_defensible`): the false-positive rate of
  each level of each analytic choice under the null, with a prespecified threshold and sensitivity values.
- **Calibrated summaries** (`spec_summary`): share of significant specifications with a permutation P value,
  the effective number of independent specifications, and approximate intervals.
- **Calibrated variance decomposition** (`spec_decompose`): which analytic choice moves the estimate, with
  pairwise interactions, equal-cardinality sensitivity, and the same decomposition under the null.
- **Constraint test for composite indices** (`composite_constraint`): whether the fixed ratio a composite
  imposes on its constituents is rejected by the data.
- **Breakpoint boundary diagnostics** (`breakpoint_boot`): naive and Rao–Wu bootstrap, share of replicates on
  the search bounds, and sensitivity to the search range.

Status: **v0.0.0.9000, all planned functions implemented** (`spec_axes`, `spec_fit`, `as_spec_grid`,
`spec_null`, `spec_defensible`, `spec_keep`, `spec_summary`, `spec_decompose`, `composite_constraint`,
`composite_coefficients`, `breakpoint_boot`, and `plot()` methods). Every function reproduces the
corresponding result of the accompanying paper from the archived ledger to machine precision
(`inst/reproduce/01`–`05`; 15,840 specifications, 300-replicate null, 1,056 constraint tests, 200
Rao–Wu replicates). 109 unit tests.

## Installation

```r
# install.packages("remotes")
remotes::install_github("duzida/speccal")
```

## Minimal example

```r
library(speccal); library(survey)
des <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC6YR, nest = TRUE, data = nh)
ax  <- spec_axes(
  outcome    = c(any = "y_cdc_any", severe = "y_cdc_severe"),
  covariates = list(crude = character(0), demo = c("Age", "Gender", "Race"),
                    full  = c("Age", "Gender", "Race", "Smoking", "BMI")),
  coding     = c("per_sd", "log2", "q4_vs_q1"),
  sample     = list(all = TRUE, age40 = quote(Age >= 40)))
g   <- spec_fit(ax, nh, design = des, exposure = c("WBC", "NLR"), cores = 8)
nul <- spec_null(g, scheme = "fl", B = 300, stratum = "SDMVSTRA", covariates = c("Age", "Gender", "Race", "Smoking", "BMI"),
                 restrict = list(weight = "design"), cores = 8)
def <- spec_defensible(nul, threshold = 0.10, within = list(covset = "full"))
keep <- spec_keep(def, g)
sm  <- spec_summary(g, nul, keep = keep)       # share significant, permutation P, effective number of specifications
vd  <- spec_decompose(g, nul, interactions = TRUE, equal_cardinality = list(axis = "outcome", k = 1))
ct  <- composite_constraint(nh, des, composite = "NLR", constituents = c(NEU = 1, LYM = -1), axes = ax)
bp  <- breakpoint_boot(nh, des, "NLR", "y_cdc_any", covariates = c("Age", "Gender", "Race", "Smoking", "BMI"))
plot(g, keep = keep); plot(def); plot(sm, nul, keep, g); plot(vd); plot(bp)
```

## Reproducing the paper

`inst/reproduce/` contains scripts that rebuild every number in the accompanying manuscript from the
public NHANES 2009–2014 files and compare them with the archived results ledger.

## Citation

Huang P, Zhang R. speccal: calibrated inference for specification curve analysis. R package. [manuscript in preparation]

## License

MIT © Huang Peng, Zhang Rui

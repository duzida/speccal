# speccal 0.1.1

Data corrections to `nhanes_perio`; no change to the package functions.

* `MON` is now the absolute monocyte count (monocyte percentage times the white
  blood cell count / 100). In 0.1.0 it held the monocyte percentage, so `LMR`,
  `SIRI` and `PIV` were wrong; these three columns are recomputed from the
  corrected count. All other columns, including the white blood cell and
  neutrophil counts, NLR, dNLR, PLR, SII and AGR, are unchanged.
* `Arthritis`: the labels of the self-reported arthritis item (MCQ160A) were
  reversed. The values are swapped and the factor levels keep their order, so
  every exposure estimate is identical and only the sign of the arthritis
  coefficient changes.
* `Educational_level`: the two lower categories are relabelled to match their
  NHANES coding ("Less than 9th grade"; "9th-12th grade or high school
  graduate"). Values and level order are unchanged.
* The covariates were imputed once, before the monocyte correction; the
  imputation was not repeated (the exposures were complete and not imputed).

# speccal 0.1.0

First release, accompanying the paper on blood-derived inflammatory markers and
periodontitis in NHANES 2009-2014.

* `spec_null()`: the confounding-preserving null is now called exposure-residual
  permutation, `scheme = "resid"` (default). `"fl"` is kept as a deprecated alias.
  The null is conditional on the residualisation model; non-linear terms such as
  `splines::ns()` may be given in `covariates` to rebuild it as a sensitivity analysis.
* `inst/reproduce/revision/`: scripts for the residualisation sensitivity analyses, the
  simulation, the independent permutations with adjustment for multiplicity, the Rao-Wu
  bootstrap intervals for shares and the sup-likelihood-ratio breakpoint test.

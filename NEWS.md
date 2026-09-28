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

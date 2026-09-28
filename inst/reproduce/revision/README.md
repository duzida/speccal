# Analyses added in revision

Scripts for the analyses added after methodological review. They run on the analysis
scripts of the paper (`spec_core.R`, `spec_core2.R`; set `NH_BASE` to the project root)
rather than on the package API, and were run with `OMP_NUM_THREADS=1` and
`mclapply` workers. Seeds are fixed inside each script.

| Script | Analysis | Seeds |
|---|---|---|
| `resid_sens.R` | Null rebuilt with three alternative residualisation models (`RMODEL=socio`, `g2core`, `spline`), 300 permutations each | 1e6/2e6/3e6 + b |
| `var_diag.R` | Source of the 6.7% false-positive rate of the full model: linearisation, jackknife and unweighted variance | 500000 + b |
| `sim_study.R` | Simulation: measured, non-linear, unmeasured and collider confounding; power | fixed per condition and replicate |
| `summarise_null.R` | False-positive rates, defensible set and permutation P for any null | - |
| `perm_ext.R` | 2,000 permutations independent of the selection of defensible levels (b = 301..2300) | 500000 + b |
| `share_boot.R` | Rao-Wu bootstrap intervals for the share of significant specifications | 20260928 |
| `suplr.R` | sup-likelihood-ratio test for a breakpoint, parametric bootstrap | 700000 + b |
| `summarise_step3.R` | Holm and minP adjustment; bootstrap and sup-LR summaries | - |

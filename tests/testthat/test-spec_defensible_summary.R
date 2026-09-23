test_that("spec_defensible flags under-adjusted levels under a confounded fl null and spec_keep selects rows", {
  d <- make_toy(n = 1500); d$x1 <- d$x1 * exp(0.06 * d$Age); d$x2 <- d$x2 * exp(0.06 * d$Age)
  d$y_any <- stats::rbinom(nrow(d), 1, stats::plogis(-2 + 0.03 * d$Age))
  ax <- spec_axes(outcome = c(any = "y_any"), covariates = list(crude = character(0), demo = c("Age", "Gender")),
                  coding = c("per_sd", "log2"), sample = list(all = TRUE, age40 = quote(Age >= 40)), weighting = "unweighted")
  g <- spec_fit(ax, d, exposure = c("x1", "x2"), marginal = FALSE)
  nul <- spec_null(g, scheme = "fl", reps = 1:8, stratum = "strat", covariates = c("Age", "Gender"), seed = 11)
  def <- spec_defensible(nul, threshold = 0.10, within = list(covset = "demo"))
  expect_s3_class(def, "spec_defensible")
  d0 <- as.data.frame(def)
  expect_true(all(c("axis", "level", "fpr", "lo", "hi", "defensible", "defensible_0.075", "defensible_0.15") %in% names(d0)))
  expect_false(d0$defensible[d0$axis == "covset" & d0$level == "crude"])   # confounded crude model: many false positives
  expect_true(d0$defensible[d0$axis == "covset" & d0$level == "demo"])
  expect_true(all(d0$lo <= d0$fpr & d0$fpr <= d0$hi))
  expect_true(all(d0$n_specs[d0$axis == "covset"] == nrow(g) / 2))
  expect_warning(spec_defensible(nul, within = list(covset = c("crude", "demo"))), "themselves not defensible")
  keep <- spec_keep(def, g)
  expect_equal(length(keep), nrow(g))
  expect_error(spec_summary(g, keep = rep(FALSE, nrow(g))), "no specifications selected")
  expect_true(all(g$covset[keep] == "demo")); expect_false(any(keep & g$covset == "crude"))
})

test_that("spec_summary computes shares, direction, Janus and permutation calibration", {
  d <- make_toy(n = 1500)
  ax <- spec_axes(outcome = c(any = "y_any", sev = "y_sev"), covariates = list(crude = character(0), demo = c("Age", "Gender")),
                  coding = c("per_sd", "log2"), weighting = "unweighted")
  g <- spec_fit(ax, d, exposure = c("x1", "x2"), marginal = FALSE)
  s0 <- spec_summary(g)
  expect_s3_class(s0, "spec_summary"); expect_equal(nrow(s0), 2)
  expect_equal(s0$n, c(8, 8)); expect_true(all(s0$sig >= s0$sig_modal))
  expect_equal(s0$direction[s0$exposure == "x1"], 1)          # x1 raises the outcome in the toy
  expect_true(s0$sig_modal[s0$exposure == "x1"] > s0$sig_modal[s0$exposure == "x2"])
  nul <- spec_null(g, scheme = "simple", reps = 1:10, stratum = "strat", seed = 5)
  s1 <- spec_summary(g, nul)
  expect_true(all(c("p_perm", "n_eff_var", "n_eff_liji", "ci_lo", "ci_hi", "null_median_modal") %in% names(s1)))
  expect_true(s1$p_perm[s1$exposure == "x1"] <= s1$p_perm[s1$exposure == "x2"])
  expect_true(all(s1$p_perm >= 1 / 11 & s1$p_perm <= 1))
  ok <- is.na(s1$ci_lo) | (s1$ci_lo <= s1$sig_modal & s1$sig_modal <= s1$ci_hi)   # NA when the null share has zero variance
  expect_true(all(ok))
  expect_true(all(s1$n_eff_liji <= 8 + 1e-9))
  ## keep restricts both grid and null
  s2 <- spec_summary(g, nul, keep = g$covset == "demo")
  expect_equal(s2$n, c(4, 4)); expect_equal(s2$n_null_specs, c(4, 4))
})

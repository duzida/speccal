test_that("spec_null runs both schemes, is reproducible by seed, and preserves the exposure-covariate association under fl", {
  d <- make_toy(n = 1500); d$x1 <- d$x1 * exp(0.06 * d$Age)   # exposure strongly confounded by Age
  d$y_any <- stats::rbinom(nrow(d), 1, stats::plogis(-2 + 0.03 * d$Age))   # outcome depends on Age only
  ax <- spec_axes(outcome = c(any = "y_any"), covariates = list(crude = character(0), demo = c("Age", "Gender")),
                  coding = c("per_sd", "log2"), sample = list(all = TRUE, age40 = quote(Age >= 40)), weighting = "unweighted")
  g <- spec_fit(ax, d, exposure = c("x1", "x2"))
  nf <- spec_null(g, scheme = "fl", reps = 1:3, stratum = "strat", covariates = c("Age", "Gender"), seed = 7)
  ns <- spec_null(g, scheme = "simple", reps = 1:3, stratum = "strat", seed = 7)
  expect_s3_class(nf, "spec_null"); expect_s3_class(ns, "spec_null")
  expect_equal(nrow(nf$perm), 3 * nrow(g)); expect_equal(nrow(ns$perm), 3 * nrow(g))
  expect_true(all(c("rep", "exposure", "outcome", "covset", "coding", "sample", "weight", "estimate", "p") %in% names(nf$perm)))
  nf2 <- spec_null(g, scheme = "fl", reps = 1:3, stratum = "strat", covariates = c("Age", "Gender"), seed = 7)
  expect_equal(nf2$perm$estimate, nf$perm$estimate)
  nf3 <- spec_null(g, scheme = "fl", reps = 1:3, stratum = "strat", covariates = c("Age", "Gender"), seed = 8)
  expect_false(isTRUE(all.equal(nf3$perm$estimate, nf$perm$estimate)))
  ## fl keeps the exposure-Age association: the crude (confounded) estimate stays positive under the fl null,
  ## whereas simple permutation removes it; the adjusted estimate is near 0 under both
  sel <- function(nn, cv) nn$perm$estimate[nn$perm$covset == cv & nn$perm$exposure == "x1" & nn$perm$coding == "log2" & nn$perm$sample == "all"]
  expect_gt(mean(sel(nf, "crude")), 0.08)
  expect_lt(abs(mean(sel(ns, "crude"))), 0.08)
  expect_lt(abs(mean(sel(nf, "demo"))), 0.08)
})

test_that("restrict_axes keeps only the requested levels and spec_null honours restrict", {
  ax <- toy_axes()
  r <- restrict_axes(ax, list(sample = "all", coding = c("per_sd", "log2")))
  expect_equal(names(r$subsets$sample), "all"); expect_equal(names(r$coding), c("per_sd", "log2"))
  expect_equal(n_specs(r), 2 * 2 * 2 * 1 * 2)
  expect_error(restrict_axes(ax, list(nope = "a")), "unknown axis")
  d <- make_toy(n = 600)
  des <- survey::svydesign(id = ~psu, strata = ~strat, weights = ~w, data = d, nest = TRUE)
  g <- spec_fit(ax, d, des, exposure = "x1", marginal = FALSE)
  nn <- spec_null(g, scheme = "simple", reps = 1:2, stratum = "strat", restrict = list(sample = "all", weight = "design"), seed = 1)
  expect_equal(nrow(nn$perm), 2 * 2 * 2 * 3)
  expect_true(all(nn$perm$weight == "design") && all(nn$perm$sample == "all"))
})

test_that("spec_null resumes from a directory", {
  d <- make_toy(n = 500)
  ax <- spec_axes(outcome = c(any = "y_any"), covariates = list(crude = character(0)), coding = "log2", weighting = "unweighted")
  g <- spec_fit(ax, d, exposure = "x1", marginal = FALSE)
  td <- tempfile(); dir.create(td)
  n1 <- spec_null(g, scheme = "simple", reps = 1:2, seed = 3, dir = td)
  expect_equal(sort(list.files(td)), c("rep_0001.csv", "rep_0002.csv"))
  n2 <- spec_null(g, scheme = "simple", reps = 1:3, seed = 3, dir = td)
  expect_equal(n2$perm$estimate[n2$perm$rep <= 2], n1$perm$estimate)
  expect_equal(length(unique(n2$perm$rep)), 3)
})

test_that("composite_constraint rejects a wrong ratio and accepts a right one; coefficients line up", {
  set.seed(3); n <- 2000
  d <- data.frame(strat = rep(1:10, each = n / 10), psu = rep(1:2, n / 2), w = stats::runif(n, .5, 2),
                  Age = round(stats::runif(n, 30, 80)), N = stats::rlnorm(n, 1.5, .4), L = stats::rlnorm(n, 0.5, .4))
  d$NLR <- d$N / d$L
  d$y_ok  <- stats::rbinom(n, 1, stats::plogis(-1 + 0.8 * (log(d$N) - log(d$L))))   # outcome follows the ratio exactly
  d$y_bad <- stats::rbinom(n, 1, stats::plogis(-1 + 0.8 * log(d$N)))                 # outcome ignores L
  des <- survey::svydesign(id = ~psu, strata = ~strat, weights = ~w, data = d, nest = TRUE)
  ax <- spec_axes(outcome = c(ok = "y_ok", bad = "y_bad"), covariates = list(crude = character(0), age = "Age"), coding = "log2")
  ct <- composite_constraint(d, des, composite = "NLR", constituents = c(N = 1, L = -1), axes = ax)
  expect_s3_class(ct, "composite_constraint"); expect_equal(nrow(ct), 4)
  expect_true(all(ct$p_constraint[ct$outcome == "bad"] < 0.01))
  expect_true(all(ct$p_constraint[ct$outcome == "ok"] > 0.05))
  expect_equal(ct$k_free, rep(1, 4))
  s <- summary(ct, exclude = list(outcome = "bad")); expect_equal(s$n_excluding, 2); expect_equal(s$rejected_excluding, 0)
  ## omitting the other constituent gives the same statistic
  ct2 <- composite_constraint(d, des, "NLR", c(N = 1, L = -1), ax, free = "L")
  expect_equal(ct2$F, ct$F, tolerance = 1e-8)
  expect_error(composite_constraint(d, des, "NLR", c(N = 1, L = 1), ax), "not the stated combination")
  cc <- composite_coefficients(d, des, "NLR", c(N = 1, L = -1), outcome = "y_ok", covariates = "Age")
  expect_equal(nrow(cc), 2); expect_equal(cc$implied[1], -cc$implied[2])
  expect_true(abs(cc$free[1] - 0.8) < 0.15 && abs(cc$free[2] + 0.8) < 0.15)
})

test_that("breakpoint_boot recovers a planted breakpoint and reports bound shares", {
  set.seed(5); n <- 3000
  d <- data.frame(strat = rep(1:15, each = n / 15), psu = rep(1:2, n / 2), w = stats::runif(n, .5, 2), x = stats::rlnorm(n, 0, .6))
  lx <- log2(d$x); d$y <- stats::rbinom(n, 1, stats::plogis(-1.5 + 0.1 * lx + 1.6 * pmax(lx - 0.5, 0)))   # break at log2(x) = 0.5
  des <- survey::svydesign(id = ~psu, strata = ~strat, weights = ~w, data = d, nest = TRUE)
  bp <- breakpoint_boot(d, des, "x", "y", grid = c(0.05, 0.95, 21), B = 30, resampling = "rao-wu", seed = 4)
  expect_s3_class(bp, "breakpoint_boot")
  expect_equal(nrow(bp$boot), 30); expect_equal(length(bp$grid), 21)
  expect_true(abs(log2(bp$observed$bp) - 0.5) < 0.4)
  expect_true(bp$summary$lo <= bp$observed$bp && bp$observed$bp <= bp$summary$hi)
  expect_true(bp$summary$at_lower_bound + bp$summary$at_upper_bound < 0.5)
  expect_true(bp$observed$p_slope_change < 0.01)
  bn <- breakpoint_boot(d, des, "x", "y", grid = c(0.05, 0.95, 11), B = 4, resampling = "naive", seed = 4, sensitivity_grid = c(0.2, 0.8, 7))
  expect_equal(nrow(bn$boot), 4); expect_equal(length(bn$sensitivity$grid), 7)
  bn2 <- breakpoint_boot(d, des, "x", "y", grid = c(0.05, 0.95, 11), B = 4, resampling = "naive", seed = 4)
  expect_equal(bn2$boot$bp, bn$boot$bp)
})

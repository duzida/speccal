make_toy <- function(n = 1500, seed = 1) {
  set.seed(seed)
  d <- data.frame(strat = rep(1:10, each = n / 10), psu = rep(1:2, n / 2), w = stats::runif(n, .5, 2),
                  Age = round(stats::runif(n, 30, 80)), Gender = factor(sample(c("F", "M"), n, TRUE)),
                  x1 = stats::rlnorm(n), x2 = stats::rlnorm(n))
  d$y_any <- stats::rbinom(n, 1, stats::plogis(-1 + .4 * log(d$x1) + .01 * d$Age))
  d$y_sev <- stats::rbinom(n, 1, stats::plogis(-2 + .3 * log(d$x1)))
  d
}
toy_axes <- function() spec_axes(
  outcome = c(any = "y_any", severe = "y_sev"),
  covariates = list(crude = character(0), demo = c("Age", "Gender")),
  coding = c("per_sd", "log2", "q4_vs_q1"),
  sample = list(all = TRUE, age40 = quote(Age >= 40)))

test_that("spec_axes counts specifications and rejects reserved names", {
  ax <- toy_axes()
  expect_s3_class(ax, "spec_axes")
  expect_equal(n_specs(ax), 2 * 2 * 3 * 2 * 2)
  expect_error(spec_axes(outcome = c(a = "y"), covariates = list(c = NULL), coding = "per_sd", outcome2 = list(a = TRUE)), NA)
  expect_error(spec_axes(outcome = c(a = "y"), covariates = list(c = NULL), coding = "nope"), "unknown built-in coding")
  expect_error(spec_axes(outcome = c(a = "y"), covariates = list(c = NULL), coding = "per_sd", covset = list(a = TRUE)), "reserved")
})

test_that("spec_fit returns one row per specification with diagnostics, and refit reproduces", {
  d <- make_toy(); ax <- toy_axes()
  des <- survey::svydesign(id = ~psu, strata = ~strat, weights = ~w, data = d, nest = TRUE)
  g <- spec_fit(ax, d, des, exposure = c("x1", "x2"))
  expect_s3_class(g, "spec_grid")
  expect_equal(nrow(g), 2 * n_specs(ax))
  expect_setequal(axis_cols(g), c("outcome", "covset", "coding", "sample", "weight"))
  expect_true(all(c("estimate", "se", "or", "p", "rr", "p_rr", "rd", "p_rd", "converged", "sep_flag") %in% names(g)))
  expect_true(all(g$converged))
  expect_equal(g$or, exp(g$estimate))
  ## domain estimation: restricted sample keeps the full design's degrees of freedom
  expect_true(all(g$design_df[g$weight == "design"] == survey::degf(des)))
  expect_true(all(is.infinite(g$design_df[g$weight == "unweighted"])))
  expect_true(all(g$n[g$sample == "age40"] < g$n[g$sample == "all"]))
  g2 <- attr(g, "refit")(d)
  expect_equal(g2$estimate, g$estimate)
  expect_equal(g2$p_rr, g$p_rr)
})

test_that("weighted point estimate equals svyglm on the domain, unweighted equals glm", {
  d <- make_toy(); ax <- toy_axes()
  des <- survey::svydesign(id = ~psu, strata = ~strat, weights = ~w, data = d, nest = TRUE)
  g <- spec_fit(ax, d, des, exposure = "x1", marginal = FALSE)
  r <- g[g$outcome == "any" & g$covset == "demo" & g$coding == "log2" & g$sample == "age40", ]
  de <- subset(stats::update(des, ..x = log2(d$x1)), d$Age >= 40)
  m <- survey::svyglm(y_any ~ ..x + Age + Gender, design = de, family = stats::quasibinomial())
  expect_equal(r$estimate[r$weight == "design"], unname(stats::coef(m)["..x"]))
  expect_equal(r$p[r$weight == "design"], unname(summary(m)$coefficients["..x", 4]))
  dd <- d[d$Age >= 40, ]; dd$..x <- log2(dd$x1)
  m2 <- stats::glm(y_any ~ ..x + Age + Gender, data = dd, family = stats::binomial())
  expect_equal(r$estimate[r$weight == "unweighted"], unname(stats::coef(m2)["..x"]))
})

test_that("marginal effects: binary contrast standardises over the domain", {
  d <- make_toy(); ax <- spec_axes(outcome = c(any = "y_any"), covariates = list(crude = character(0)),
                                  coding = "q4_vs_q1", weighting = "unweighted")
  g <- spec_fit(ax, d, exposure = "x1")
  expect_equal(nrow(g), 1)
  expect_true(is.finite(g$rr) && g$rr > 0)
  expect_equal(sign(g$rd), sign(g$estimate))
})

test_that("as_spec_grid wraps external results", {
  df <- expand.grid(a = c("a1", "a2"), b = c("b1", "b2", "b3"), stringsAsFactors = FALSE)
  df$est <- stats::rnorm(6); df$se <- 0.1; df$pv <- stats::runif(6)
  g <- as_spec_grid(df, axes = c("a", "b"), estimate = "est", se = "se", p = "pv")
  expect_s3_class(g, "spec_grid"); expect_equal(nrow(g), 6); expect_equal(g$or, exp(df$est))
  expect_equal(axis_cols(g), c("a", "b"))
})

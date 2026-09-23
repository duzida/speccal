#' Test the fixed-ratio constraint implied by a composite index
#'
#' On the log scale a composite index such as a ratio or product of measured
#' quantities is a linear combination of its constituents with coefficients
#' fixed at +1 or -1 (for example log NLR = log N - log L). Entering the
#' composite in a regression therefore equals entering the constituents with
#' their coefficients constrained to that ratio. The constraint is testable:
#' the free constituents (all but one) are added to a model that already
#' contains the composite and tested jointly with a design-based Wald test.
#' Rejection means the data prefer a different weighting of the constituents
#' than the composite imposes. The test is a within-sample comparison of
#' nested models; a composite can fail it and still predict another outcome
#' well.
#'
#' @param data Data frame.
#' @param design A `survey.design` on `data` (row order must match).
#' @param composite Name of the composite column.
#' @param constituents Named numeric vector of constituent column names and
#'   their weights on the transformed scale, e.g. `c(NEU = 1, LYM = -1)`.
#'   The identity `transform(composite) == sum(w * transform(constituent))` is
#'   checked.
#' @param axes A [spec_axes()] object; its `outcome`, `covariates` and subset
#'   axes define the specifications in which the constraint is tested (coding
#'   and weighting axes are ignored: the composite enters on the transformed
#'   scale, design-based).
#' @param transform Transformation under which the composite is linear in its
#'   constituents (default `log`).
#' @param free Names of the constituents entered freely (default: all but the
#'   last). The Wald statistic does not depend on which one is omitted.
#' @param alpha Significance level for the summary.
#' @param cores Cores for `parallel::mclapply` over specifications.
#' @return A data frame of class `composite_constraint` with one row per
#'   specification: axis columns, `k_free`, `F`, `p_constraint`,
#'   `p_composite` (P value of the composite in the model without the free
#'   constituents), `sep_flag`; attribute `summary` (rejection rate, rejection
#'   rate excluding levels named in `exclude` of `summary()`, composite
#'   significance rate).
#' @export
composite_constraint <- function(data, design, composite, constituents, axes, transform = log,
                                 free = NULL, alpha = 0.05, cores = 1L) {
  stopifnot(is.data.frame(data), inherits(design, "survey.design"), inherits(axes, "spec_axes"),
            is.numeric(constituents), !is.null(names(constituents)), length(constituents) >= 2)
  cn <- names(constituents)
  for (v in c(composite, cn)) if (is.null(data[[v]])) stop("column not found: ", v)
  tv <- paste0(".t_", c(composite, cn))
  for (i in seq_along(tv)) data[[tv[i]]] <- transform(data[[c(composite, cn)[i]]])
  lc <- as.matrix(data[, tv[-1], drop = FALSE]) %*% constituents
  dev <- max(abs(data[[tv[1]]] - lc), na.rm = TRUE)
  if (dev > 1e-6) stop(sprintf("composite is not the stated combination of its constituents on the transformed scale (max deviation %.3g)", dev))
  if (is.null(free)) free <- cn[-length(cn)]
  stopifnot(all(free %in% cn), length(free) == length(cn) - 1)
  tfree <- paste0(".t_", free); tcomp <- tv[1]
  design$variables <- data
  subs <- .subset_levels(axes, data)
  sub_grid <- if (length(subs)) expand.grid(lapply(subs, names), stringsAsFactors = FALSE) else data.frame(row.names = 1L)[1L, , drop = FALSE]
  specs <- list()
  for (si in seq_len(max(1L, nrow(sub_grid)))) for (ou in names(axes$outcome)) for (cv in names(axes$covariates))
    specs[[length(specs) + 1]] <- list(si = si, ou = ou, cv = cv)
  one <- function(sp) {
    keep <- rep(TRUE, nrow(data)); lab <- list()
    for (nm in names(subs)) { lv <- sub_grid[[nm]][sp$si]; keep <- keep & subs[[nm]][[lv]]; lab[[nm]] <- lv }
    y <- axes$outcome[[sp$ou]]; cvs <- axes$covariates[[sp$cv]]
    keep <- keep & !is.na(data[[y]])
    de <- subset(design, keep)
    cvs_txt <- if (length(cvs)) paste("+", paste(cvs, collapse = "+")) else ""
    m  <- survey::svyglm(stats::as.formula(paste(y, "~", tcomp, "+", paste(tfree, collapse = "+"), cvs_txt)), design = de, family = stats::quasibinomial())
    tt <- survey::regTermTest(m, stats::as.formula(paste("~", paste(tfree, collapse = "+"))), method = "Wald")
    m1 <- survey::svyglm(stats::as.formula(paste(y, "~", tcomp, cvs_txt)), design = de, family = stats::quasibinomial())
    r <- data.frame(composite = composite, outcome = sp$ou, covset = sp$cv, stringsAsFactors = FALSE)
    for (nm in names(lab)) r[[nm]] <- lab[[nm]]
    r$k_free <- length(free); r$F <- as.numeric(tt$Ftest); r$p_constraint <- as.numeric(tt$p)
    r$p_composite <- summary(m1)$coefficients[tcomp, 4]; r$n <- sum(keep)
    r$sep_flag <- isTRUE(max(abs(stats::coef(m))) > 10)
    r
  }
  res <- if (cores > 1) parallel::mclapply(specs, function(s) try(one(s), silent = TRUE), mc.cores = cores) else lapply(specs, function(s) try(one(s), silent = TRUE))
  bad <- vapply(res, function(r) inherits(r, "try-error"), TRUE)
  out <- do.call(rbind, res[!bad]); rownames(out) <- NULL
  if (any(bad)) warning(sum(bad), " specification(s) failed to fit and were dropped")
  structure(out, class = c("composite_constraint", "data.frame"), alpha = alpha,
            axis_cols = c("outcome", "covset", names(axes$subsets)), constituents = constituents, free = free)
}

#' @export
summary.composite_constraint <- function(object, exclude = NULL, ...) {
  a <- attr(object, "alpha"); d <- as.data.frame(object)
  keep <- rep(TRUE, nrow(d))
  if (!is.null(exclude)) for (nm in names(exclude)) keep <- keep & !d[[nm]] %in% exclude[[nm]]
  data.frame(composite = d$composite[1], n_tests = nrow(d), rejected = mean(d$p_constraint < a),
             rejected_excluding = if (is.null(exclude)) NA else mean(d$p_constraint[keep] < a),
             n_excluding = sum(keep), median_F = stats::median(d$F), composite_significant = mean(d$p_composite < a),
             stringsAsFactors = FALSE)
}

#' @export
print.composite_constraint <- function(x, ...) {
  s <- summary(x)
  cat(sprintf("Constraint test for %s = %s (transformed scale), %d specifications\n", s$composite,
              paste(sprintf("%+d*%s", as.integer(attr(x, "constituents")), names(attr(x, "constituents"))), collapse = " "), s$n_tests))
  cat(sprintf("  constraint rejected in %.1f%% (median F %.2f); composite itself significant in %.1f%%\n",
              100 * s$rejected, s$median_F, 100 * s$composite_significant))
  invisible(x)
}

#' Constituent coefficients: imposed by the composite versus estimated freely
#'
#' Fits, for one specification, the model with the composite alone and the
#' model with all constituents entered freely, and returns for each constituent
#' the coefficient the composite implies (its coefficient times the fixed
#' weight) and the freely estimated coefficient, with standard errors.
#' @inheritParams composite_constraint
#' @param outcome Outcome column name.
#' @param covariates Character vector of covariate names.
#' @param subset Optional logical vector selecting the analysis domain.
#' @return A data frame with one row per constituent.
#' @export
composite_coefficients <- function(data, design, composite, constituents, outcome, covariates = character(0),
                                   transform = log, subset = NULL) {
  cn <- names(constituents); tv <- paste0(".t_", c(composite, cn))
  for (i in seq_along(tv)) data[[tv[i]]] <- transform(data[[c(composite, cn)[i]]])
  design$variables <- data
  de <- if (is.null(subset)) design else subset(design, subset)
  cvs <- if (length(covariates)) paste("+", paste(covariates, collapse = "+")) else ""
  m1 <- survey::svyglm(stats::as.formula(paste(outcome, "~", tv[1], cvs)), design = de, family = stats::quasibinomial())
  m2 <- survey::svyglm(stats::as.formula(paste(outcome, "~", paste(tv[-1], collapse = "+"), cvs)), design = de, family = stats::quasibinomial())
  b1 <- stats::coef(m1)[tv[1]]; se1 <- sqrt(stats::vcov(m1)[tv[1], tv[1]])
  cf <- summary(m2)$coefficients[tv[-1], , drop = FALSE]
  data.frame(composite = composite, constituent = cn, weight = unname(constituents),
             implied = unname(b1 * constituents), implied_se = unname(se1 * abs(constituents)),
             free = unname(cf[, 1]), free_se = unname(cf[, 2]), free_p = unname(cf[, 4]), df_resid = m2$df.residual,
             stringsAsFactors = FALSE)
}

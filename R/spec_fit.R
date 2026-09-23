#' Fit every specification of a grid
#'
#' Fits one logistic model per specification and exposure. Survey-weighted
#' specifications use `survey::svyglm` on the full design with subsets applied
#' as domains (`subset(design, keep)`), so that design-based variances are
#' correct for restricted samples. Every fit records convergence, iteration
#' count, the largest absolute coefficient and standard error, and a
#' quasi-separation flag (either exceeding 10). Marginal risk differences and
#' risk ratios are obtained by standardisation over the domain with delta-method
#' variances.
#'
#' @param axes A [spec_axes()] object.
#' @param data Data frame containing outcomes, covariates, exposures and design
#'   variables.
#' @param design A `survey.design` object built on the full `data` (row order
#'   must match), or `NULL` when only unweighted specifications are requested.
#' @param exposure Character vector of exposure column names.
#' @param marginal Logical; compute marginal RD and RR (default `TRUE`).
#' @param cores Number of cores for `parallel::mclapply` over exposure x coding
#'   (set BLAS threads to 1 when `cores > 1`).
#' @param quiet Suppress progress messages.
#' @return A data frame of class `spec_grid` with one row per specification and
#'   exposure, columns `exposure`, the axis columns, `estimate` (log odds ratio),
#'   `se`, `or`, `lo`, `hi`, `p`, `rd`, `rd_se`, `p_rd`, `rr`, `lrr_se`, `p_rr`,
#'   `base_risk`, `prevalence`, `n`, `design_df`, `converged`, `iter`,
#'   `max_abs_b`, `max_abs_se`, `sep_flag`. Attributes `axes`, `exposure`,
#'   `refit` (a function `refit(newdata, axes_new, cores)` that repeats the fit
#'   on new data with the same design structure, optionally on a restricted set
#'   of axes; used by [spec_null()]) and `context` (the data and design).
#' @export
spec_fit <- function(axes, data, design = NULL, exposure, marginal = TRUE, cores = 1L, quiet = FALSE) {
  stopifnot(inherits(axes, "spec_axes"), is.data.frame(data), is.character(exposure), length(exposure) >= 1)
  if ("design" %in% axes$weighting) {
    if (is.null(design)) stop("design is required for weighted specifications")
    stopifnot(inherits(design, "survey.design"), nrow(design$variables) == nrow(data))
  }
  if (cores > 1 && .blas_threads_unset() && !quiet)
    message("cores > 1: set OMP_NUM_THREADS=1 (and OPENBLAS/MKL) to avoid oversubscription")
  jobs <- expand.grid(exposure = exposure, coding = names(axes$coding), stringsAsFactors = FALSE)
  one <- function(j) .fit_exposure_coding(axes, data, design, jobs$exposure[j], jobs$coding[j], marginal)
  res <- if (cores > 1) parallel::mclapply(seq_len(nrow(jobs)), one, mc.cores = cores, mc.preschedule = FALSE)
         else lapply(seq_len(nrow(jobs)), one)
  bad <- vapply(res, function(r) inherits(r, "try-error") || is.null(r), TRUE)
  if (any(bad)) stop("fitting failed for ", sum(bad), " exposure x coding job(s): ",
                     paste(jobs$exposure[bad], jobs$coding[bad], sep = "/", collapse = ", "))
  out <- do.call(rbind, res); rownames(out) <- NULL
  out <- .order_grid(out, axes, exposure)
  ax_cols <- c("outcome", "covset", "coding", names(axes$subsets), "weight")
  refit <- function(newdata, axes_new = axes, cores = 1L) {
    stopifnot(nrow(newdata) == nrow(data))
    des <- design
    if (!is.null(des)) des$variables <- newdata
    spec_fit(axes_new, newdata, des, exposure, marginal = marginal, cores = cores, quiet = TRUE)
  }
  structure(out, class = c("spec_grid", "data.frame"), axes = axes, exposure = exposure,
            axis_cols = ax_cols, refit = refit, context = list(data = data, design = design))
}

.blas_threads_unset <- function() {
  all(!nzchar(Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"))))
}

.order_grid <- function(out, axes, exposure) {
  key <- list(match(out$exposure, exposure), match(out$coding, names(axes$coding)))
  for (nm in names(axes$subsets)) key[[length(key) + 1]] <- match(out[[nm]], names(axes$subsets[[nm]]))
  key <- c(key, list(match(out$outcome, names(axes$outcome)), match(out$covset, names(axes$covariates)),
                     match(out$weight, axes$weighting)))
  out[do.call(order, key), , drop = FALSE]
}

## one exposure x one coding: the coded variable is computed once on the full sample
.fit_exposure_coding <- function(axes, data, design, ex, cd, marginal) {
  v <- data[[ex]]; if (is.null(v)) stop("exposure '", ex, "' not in data")
  x <- axes$coding[[cd]]$fun(v); ctr <- axes$coding[[cd]]$contrast
  des_x <- if (!is.null(design)) stats::update(design, ..x = x) else NULL
  subs <- .subset_levels(axes, data)
  sub_grid <- if (length(subs)) expand.grid(lapply(subs, names), stringsAsFactors = FALSE)
              else data.frame(row.names = 1L)[1L, , drop = FALSE]
  out <- list(); k <- 0L
  for (si in seq_len(max(1L, nrow(sub_grid)))) {
    keep_sub <- rep(TRUE, nrow(data)); sub_lab <- list()
    for (nm in names(subs)) { lv <- sub_grid[[nm]][si]; keep_sub <- keep_sub & subs[[nm]][[lv]]; sub_lab[[nm]] <- lv }
    for (ou in names(axes$outcome)) {
      y <- axes$outcome[[ou]]
      keep <- keep_sub & !is.na(x) & !is.na(data[[y]])
      for (cv in names(axes$covariates)) {
        f <- stats::as.formula(paste(y, "~ ..x", if (length(axes$covariates[[cv]]))
          paste("+", paste(axes$covariates[[cv]], collapse = "+")) else ""))
        for (w in axes$weighting) {
          r <- try(.fit_one(f, y, x, keep, data, des_x, w, ctr, marginal), silent = TRUE)
          if (inherits(r, "try-error")) next
          k <- k + 1L
          lab <- data.frame(exposure = ex, outcome = ou, covset = cv, coding = cd, stringsAsFactors = FALSE)
          for (nm in names(sub_lab)) lab[[nm]] <- sub_lab[[nm]]
          lab$weight <- w
          out[[k]] <- cbind(lab, r)
        }
      }
    }
  }
  do.call(rbind, out)
}

.subset_levels <- function(axes, data) {
  lapply(axes$subsets, function(ax) lapply(ax, function(e) {
    v <- if (is.logical(e) && length(e) == 1L) rep(e, nrow(data))
         else if (is.logical(e)) e
         else eval(e, data, parent.frame())
    v <- as.logical(v); v[is.na(v)] <- FALSE
    if (length(v) != nrow(data)) stop("subset expression does not evaluate to one value per row"); v }))
}

.fit_one <- function(f, y, x, keep, data, des_x, w, ctr, marginal) {
  if (w == "design") {
    de <- subset(des_x, keep)
    m  <- survey::svyglm(f, design = de, family = stats::quasibinomial())
    ww <- stats::weights(de); df <- survey::degf(de)
    prev <- as.numeric(survey::svymean(stats::as.formula(paste0("~", y)), de))
  } else {
    dd <- data[keep, , drop = FALSE]; dd$..x <- x[keep]
    m  <- stats::glm(f, data = dd, family = stats::binomial())
    ww <- rep(1, nrow(stats::model.matrix(m))); df <- Inf
    prev <- mean(dd[[y]])
  }
  cf <- summary(m)$coefficients
  if (!"..x" %in% rownames(cf)) stop("exposure dropped from model")
  s <- cf["..x", ]
  mf <- if (marginal) marginal_effects(m, ww, stats::vcov(m), contrast = ctr)
        else c(rd = NA, rd_se = NA, lrr = NA, lrr_se = NA, mu0 = NA)
  tstat <- function(est, se) if (is.finite(df)) 2 * stats::pt(-abs(est / se), df) else 2 * stats::pnorm(-abs(est / se))
  data.frame(estimate = unname(s[1]), se = unname(s[2]),
             or = exp(unname(s[1])), lo = exp(unname(s[1] - 1.96 * s[2])), hi = exp(unname(s[1] + 1.96 * s[2])), p = unname(s[4]),
             rd = unname(mf["rd"]), rd_se = unname(mf["rd_se"]), p_rd = unname(tstat(mf["rd"], mf["rd_se"])),
             rr = exp(unname(mf["lrr"])), lrr_se = unname(mf["lrr_se"]), p_rr = unname(tstat(mf["lrr"], mf["lrr_se"])),
             base_risk = unname(mf["mu0"]), prevalence = prev, n = sum(keep), design_df = df,
             converged = isTRUE(m$converged), iter = as.integer(m$iter),
             max_abs_b = max(abs(cf[, 1])), max_abs_se = max(abs(cf[, 2])),
             sep_flag = isTRUE(max(abs(cf[, 1])) > 10 | max(abs(cf[, 2])) > 10))
}

#' Marginal risk difference and risk ratio by standardisation
#'
#' Standardises fitted probabilities over the (weighted) analysis population
#' under a one-unit shift or a 0/1 contrast of the exposure, with delta-method
#' variances from the model's covariance matrix.
#' @param m A fitted `glm` or `svyglm`.
#' @param w Weights of the rows in `model.matrix(m)`.
#' @param V Covariance matrix of the coefficients.
#' @param xname Name of the exposure column in the model matrix.
#' @param contrast `"shift"` (x + 1) or `"binary"` (x = 1 vs x = 0).
#' @return Named vector `rd`, `rd_se`, `lrr`, `lrr_se`, `mu0`.
#' @export
marginal_effects <- function(m, w, V, xname = "..x", contrast = c("shift", "binary")) {
  contrast <- match.arg(contrast)
  X0 <- stats::model.matrix(m); b <- stats::coef(m)
  if (!xname %in% colnames(X0) || nrow(X0) != length(w)) return(c(rd = NA, rd_se = NA, lrr = NA, lrr_se = NA, mu0 = NA))
  X1 <- X0
  if (contrast == "shift") X1[, xname] <- X0[, xname] + 1 else { X0[, xname] <- 0; X1[, xname] <- 1 }
  ok <- is.finite(X0 %*% b) & is.finite(X1 %*% b)
  if (!all(ok)) { X0 <- X0[ok, , drop = FALSE]; X1 <- X1[ok, , drop = FALSE]; w <- w[ok] }
  p0 <- stats::plogis(drop(X0 %*% b)); p1 <- stats::plogis(drop(X1 %*% b)); sw <- sum(w)
  mu0 <- sum(w * p0) / sw; mu1 <- sum(w * p1) / sw
  g0 <- colSums(w * p0 * (1 - p0) * X0) / sw; g1 <- colSums(w * p1 * (1 - p1) * X1) / sw
  g_rd <- g1 - g0; g_lrr <- g1 / mu1 - g0 / mu0
  c(rd = mu1 - mu0, rd_se = sqrt(max(0, drop(t(g_rd) %*% V %*% g_rd))),
    lrr = log(mu1) - log(mu0), lrr_se = sqrt(max(0, drop(t(g_lrr) %*% V %*% g_lrr))), mu0 = mu0)
}

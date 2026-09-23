#' Breakpoint of a two-segment model with bootstrap boundary diagnostics
#'
#' Fits a two-segment (broken-stick) logistic model on the transformed
#' exposure, choosing the breakpoint by minimum deviance over a fixed grid of
#' candidates at equally spaced percentiles of the exposure in the analysis
#' sample. The bootstrap distribution of the selected breakpoint is reported
#' together with the share of replicates that selected the first or last
#' candidate: a breakpoint whose interval is set by the search range rather
#' than by the data piles up at the bounds. Two resampling schemes are
#' available: `"rao-wu"` (Rao-Wu rescaled bootstrap of primary sampling units
#' within strata via `survey::as.svrepdesign(type = "subbootstrap")`, fitted by
#' weighted maximum likelihood) and `"naive"` (primary sampling units resampled
#' with replacement within strata without rescaling, refitted with
#' `survey::svyglm`; understates the variance when strata contain few units,
#' provided for comparison).
#'
#' @param data Data frame.
#' @param design A `survey.design` on `data`.
#' @param exposure,outcome Column names.
#' @param covariates Character vector of covariate names.
#' @param transform Transformation of the exposure (default `log2`).
#' @param grid `c(lower percentile, upper percentile, number of candidates)`.
#' @param B Number of bootstrap replicates.
#' @param resampling `"rao-wu"` or `"naive"`.
#' @param seed Integer seed (`set.seed(seed)` before the Rao-Wu replicate
#'   weights are drawn; `set.seed(seed + b)` before naive replicate `b`).
#' @param cores Cores for `parallel::mclapply` over replicates.
#' @param sensitivity_grid Optional second grid, e.g. `c(0.15, 0.85, 25)`, on
#'   which the analysis is repeated to show how the interval moves with the
#'   search range.
#' @return A list of class `breakpoint_boot`: `observed` (breakpoint on the
#'   original scale, whether it lies on a bound, deviance improvement over the
#'   linear model and the P value of the slope change), `boot` (one row per
#'   replicate: `bp`, `at_lo`, `at_hi`, `ddev`), `summary` (percentile
#'   interval, shares at the bounds, share of replicates in which the
#'   two-segment model improved on the linear one), `grid` (candidate values on
#'   the original scale), and the same for `sensitivity` when requested.
#' @export
breakpoint_boot <- function(data, design, exposure, outcome, covariates = character(0), transform = log2,
                            grid = c(0.05, 0.95, 41), B = 200L, resampling = c("rao-wu", "naive"), seed = 1L,
                            cores = 1L, sensitivity_grid = NULL) {
  resampling <- match.arg(resampling)
  stopifnot(inherits(design, "survey.design"), nrow(design$variables) == nrow(data))
  run <- function(gr) {
    v <- transform(data[[exposure]]); y <- data[[outcome]]
    cand <- as.numeric(stats::quantile(v, seq(gr[1], gr[2], length.out = gr[3]), na.rm = TRUE))
    inv <- function(c) if (identical(transform, log2)) 2^c else if (identical(transform, log)) exp(c) else if (identical(transform, log10)) 10^c else c
    ## observed: design-based fits with the slope-change test
    d <- data; d$..x <- v; des <- design; des$variables <- d
    f0 <- stats::as.formula(paste(outcome, "~ ..x", if (length(covariates)) paste("+", paste(covariates, collapse = "+")) else ""))
    m0 <- survey::svyglm(f0, design = des, family = stats::quasibinomial()); dev0 <- m0$deviance
    best <- NULL
    for (i in seq_along(cand)) {
      des2 <- stats::update(des, ..seg = pmax(v - cand[i], 0))
      m <- try(survey::svyglm(stats::update(f0, . ~ . + ..seg), design = des2, family = stats::quasibinomial()), silent = TRUE)
      if (inherits(m, "try-error")) next
      if (is.null(best) || m$deviance < best$dev) best <- list(dev = m$deviance, i = i, p = summary(m)$coefficients["..seg", 4])
    }
    observed <- data.frame(bp = inv(cand[best$i]), at_lo = best$i == 1, at_hi = best$i == length(cand),
                           ddev = dev0 - best$dev, p_slope_change = best$p)
    ## bootstrap
    Xc <- stats::model.matrix(stats::as.formula(paste("~", if (length(covariates)) paste(covariates, collapse = "+") else "1")), data = data)
    fit_w <- function(w) {                        # weighted ML on the fixed grid (Rao-Wu)
      w <- w / mean(w); X0 <- cbind(Xc, v)
      m0 <- stats::glm.fit(X0, y, weights = w, family = stats::quasibinomial(), control = stats::glm.control(maxit = 100)); dv0 <- m0$deviance
      bs <- NULL
      for (i in seq_along(cand)) {
        m <- try(stats::glm.fit(cbind(X0, pmax(v - cand[i], 0)), y, weights = w, family = stats::quasibinomial(), control = stats::glm.control(maxit = 100)), silent = TRUE)
        if (inherits(m, "try-error") || !m$converged) next
        if (is.null(bs) || m$deviance < bs$dev) bs <- list(dev = m$deviance, i = i)
      }
      if (is.null(bs)) return(c(bp = NA, at_lo = NA, at_hi = NA, ddev = NA))
      c(bp = inv(cand[bs$i]), at_lo = as.numeric(bs$i == 1), at_hi = as.numeric(bs$i == length(cand)), ddev = dv0 - bs$dev)
    }
    if (resampling == "rao-wu") {
      set.seed(seed)
      rep <- survey::as.svrepdesign(design, type = "subbootstrap", replicates = B)
      W <- stats::weights(rep, type = "analysis")
      res <- if (cores > 1) parallel::mclapply(seq_len(B), function(b) fit_w(W[, b]), mc.cores = cores) else lapply(seq_len(B), function(b) fit_w(W[, b]))
    } else {
      st <- design$strata[, 1]; psu <- design$cluster[, 1]
      fit_naive <- function(b) {
        set.seed(seed + b)
        parts <- lapply(split(seq_len(nrow(data)), st), function(s) {        # one draw of units per stratum, relabelled 1..n_h
          ps <- unique(psu[s]); pk <- sample(ps, length(ps), replace = TRUE)
          list(idx = unlist(lapply(seq_along(pk), function(i) s[psu[s] == pk[i]])),
               id  = unlist(lapply(seq_along(pk), function(i) rep(i, sum(psu[s] == pk[i]))))) })
        idx <- unlist(lapply(parts, `[[`, "idx")); newid <- unlist(lapply(parts, `[[`, "id"))
        db <- data[idx, ]; db$..x <- v[idx]
        db$..psu <- newid; db$..st <- st[idx]; db$..w <- stats::weights(design)[idx]
        deb <- survey::svydesign(id = ~..psu, strata = ~..st, weights = ~..w, data = db, nest = TRUE)
        mb0 <- try(survey::svyglm(f0, design = deb, family = stats::quasibinomial()), silent = TRUE)
        dv0 <- if (inherits(mb0, "try-error")) NA else mb0$deviance
        bs <- NULL
        for (i in seq_along(cand)) {
          deb2 <- stats::update(deb, ..seg = pmax(db$..x - cand[i], 0))
          m <- try(survey::svyglm(stats::update(f0, . ~ . + ..seg), design = deb2, family = stats::quasibinomial()), silent = TRUE)
          if (inherits(m, "try-error")) next
          if (is.null(bs) || m$deviance < bs$dev) bs <- list(dev = m$deviance, i = i)
        }
        if (is.null(bs)) return(c(bp = NA, at_lo = NA, at_hi = NA, ddev = NA))
        c(bp = inv(cand[bs$i]), at_lo = as.numeric(bs$i == 1), at_hi = as.numeric(bs$i == length(cand)), ddev = dv0 - bs$dev)
      }
      res <- if (cores > 1) parallel::mclapply(seq_len(B), fit_naive, mc.cores = cores, mc.preschedule = FALSE) else lapply(seq_len(B), fit_naive)
    }
    boot <- as.data.frame(do.call(rbind, res)); boot$b <- seq_len(B)
    ok <- is.finite(boot$bp)
    summ <- data.frame(estimate = observed$bp, lo = unname(stats::quantile(boot$bp[ok], .025)), hi = unname(stats::quantile(boot$bp[ok], .975)),
                       at_lower_bound = mean(boot$at_lo[ok]), at_upper_bound = mean(boot$at_hi[ok]),
                       segment_improves = mean(boot$ddev[ok] > 0), n_failed = sum(!ok), B = B)
    list(observed = observed, boot = boot, summary = summ, grid = inv(cand), grid_spec = gr)
  }
  out <- run(grid)
  if (!is.null(sensitivity_grid)) out$sensitivity <- run(sensitivity_grid)
  structure(c(out, list(exposure = exposure, outcome = outcome, resampling = resampling, transform = transform)), class = "breakpoint_boot")
}

#' @export
print.breakpoint_boot <- function(x, ...) {
  s <- x$summary; o <- x$observed; g <- x$grid_spec
  cat(sprintf("Breakpoint for %s -> %s on a %d-point grid (percentiles %g-%g), %s bootstrap, B = %d\n",
              x$exposure, x$outcome, g[3], 100 * g[1], 100 * g[2], x$resampling, s$B))
  cat(sprintf("  estimate %.4g (95%% interval %.4g-%.4g); slope-change P = %.3g; deviance improvement %.2f\n",
              s$estimate, s$lo, s$hi, o$p_slope_change, o$ddev))
  cat(sprintf("  replicates on lower bound %.1f%%, on upper bound %.1f%%; two-segment model improved on linear in %.1f%%\n",
              100 * s$at_lower_bound, 100 * s$at_upper_bound, 100 * s$segment_improves))
  if (!is.null(x$sensitivity)) { t <- x$sensitivity$summary; gs <- x$sensitivity$grid_spec
    cat(sprintf("  sensitivity grid (percentiles %g-%g, %d points): interval %.4g-%.4g, lower bound %.1f%%, upper bound %.1f%%\n",
                100 * gs[1], 100 * gs[2], gs[3], t$lo, t$hi, 100 * t$at_lower_bound, 100 * t$at_upper_bound)) }
  invisible(x)
}

#' Calibrated summary of a specification curve
#'
#' For each exposure: the share of specifications significant on the odds-ratio
#' scale (and on the marginal risk-ratio and risk-difference scales when
#' available), the modal direction and the share of estimates on that side,
#' the share significant in the modal direction, and the vibration-of-effects
#' summaries (1st and 99th percentile odds ratios, relative odds ratio,
#' relative P value, Janus effect). With a [spec_null()], the share significant
#' in the modal direction is compared with its null distribution (permutation
#' P value), and the effective number of independent specifications is
#' estimated in two ways: from the variance of the per-replicate share of
#' significant specifications under the null (`n_eff_var = p(1 - p) / Var`)
#' and from the eigenvalues of the correlation matrix of the estimates across
#' replicates (Li and Ji, 2005). An approximate interval for the observed share
#' uses `n_eff_var` in place of the nominal denominator; the permutation P value
#' remains the formal test.
#'
#' @param grid A `spec_grid`.
#' @param null Optional [spec_null()] computed on the same grid (or a subgrid).
#' @param keep Optional logical vector over the rows of `grid` selecting the
#'   specifications to summarise (for example from [spec_keep()]). Rows not
#'   selected are ignored; with a null, the null is restricted to the same
#'   specifications.
#' @param alpha Significance level.
#' @param direction `"modal"` (majority sign of the estimates among the
#'   summarised specifications), `"positive"` or `"negative"`.
#' @return A data frame of class `spec_summary`, one row per exposure.
#' @references Li J, Ji L. Adjusting multiple testing in multilocus analyses
#'   using the eigenvalues of a correlation matrix. Heredity 2005;95:221-227.
#' @export
spec_summary <- function(grid, null = NULL, keep = NULL, alpha = 0.05, direction = c("modal", "positive", "negative")) {
  stopifnot(inherits(grid, "spec_grid"))
  direction <- match.arg(direction)
  G <- as.data.frame(grid)
  if (!is.null(keep)) { stopifnot(length(keep) == nrow(G)); G <- G[keep, ] }
  ax <- attr(grid, "axis_cols")
  G$.key <- do.call(paste, c(G[ax], sep = "\r"))
  out <- lapply(unique(G$exposure), function(e) {
    s <- G[G$exposure == e, ]
    dir <- switch(direction, modal = if (mean(s$estimate > 0) >= 0.5) 1 else -1, positive = 1, negative = -1)
    hit <- s$p < alpha & sign(s$estimate) == dir
    q <- stats::quantile(s$or, c(.01, .99)); lp <- -log10(s$p); qp <- stats::quantile(lp, c(.01, .99))
    r <- data.frame(exposure = e, n = nrow(s), direction = dir, modal_share = mean(sign(s$estimate) == dir),
                    sig = mean(s$p < alpha), sig_modal = mean(hit),
                    sig_rr = if ("p_rr" %in% names(s)) mean(s$p_rr < alpha, na.rm = TRUE) else NA,
                    sig_rd = if ("p_rd" %in% names(s)) mean(s$p_rd < alpha, na.rm = TRUE) else NA,
                    or_p01 = unname(q[1]), or_p99 = unname(q[2]), relative_or = unname(q[2] / q[1]),
                    relative_p = unname(qp[2] - qp[1]), janus = unname(q[1] < 1 & q[2] > 1), stringsAsFactors = FALSE)
    if (!is.null(null)) {
      P <- null$perm[null$perm$exposure == e, ]
      P$.key <- do.call(paste, c(P[ax], sep = "\r"))
      P <- P[P$.key %in% s$.key, ]
      if (!nrow(P)) stop("null contains none of the summarised specifications for ", e)
      nb <- length(unique(P$rep))
      byrep <- split(P, P$rep)
      h <- vapply(byrep, function(b) mean(b$p < alpha & sign(b$estimate) == dir), 1)   # null share in modal direction
      sg <- vapply(byrep, function(b) mean(b$p < alpha), 1)                            # null share significant
      pbar <- mean(sg); v <- stats::var(sg)
      n_eff_var <- if (v > 0) pbar * (1 - pbar) / v else NA
      M <- do.call(cbind, lapply(byrep, function(b) b$estimate[order(b$.key)]))         # specs x reps
      n_eff_lj <- .li_ji(t(M))
      se <- sqrt(r$sig_modal * (1 - r$sig_modal) / n_eff_var)
      r$n_null_specs <- nrow(P) / nb; r$B <- nb
      r$null_mean_sig <- pbar; r$null_sd_sig <- sqrt(v)
      r$null_median_modal <- stats::median(h); r$null_q95_modal <- unname(stats::quantile(h, .95))
      r$p_perm <- (sum(h >= r$sig_modal - 1e-9) + 1) / (nb + 1)
      r$n_eff_var <- n_eff_var; r$n_eff_liji <- n_eff_lj
      r$ci_lo <- max(0, r$sig_modal - 1.96 * se); r$ci_hi <- min(1, r$sig_modal + 1.96 * se)
    }
    r
  })
  out <- do.call(rbind, out); rownames(out) <- NULL
  structure(out, class = c("spec_summary", "data.frame"), alpha = alpha, direction = direction,
            has_null = !is.null(null), n_total = nrow(G))
}

.li_ji <- function(M) {   # M: replicates x specifications
  if (ncol(M) < 2 || nrow(M) < 3) return(NA_real_)
  C <- stats::cor(M); C[is.na(C)] <- 0
  ev <- pmax(eigen(C, symmetric = TRUE, only.values = TRUE)$values, 0)
  sum(as.integer(ev >= 1) + (ev - floor(ev)))
}

#' @export
print.spec_summary <- function(x, ...) {
  d <- as.data.frame(x); a <- attr(x, "alpha")
  cat(sprintf("Specification curve summary (P < %g), %d specifications\n", a, attr(x, "n_total")))
  fmt <- function(v) sprintf("%5.1f", 100 * v)
  hdr <- sprintf("  %-10s %5s %8s %8s %6s %6s %5s", "exposure", "n", "sig%", "modal%", "ROR", "Janus", "dir")
  if (attr(x, "has_null")) hdr <- paste(hdr, sprintf("%8s %7s %10s %12s", "null-med", "P", "n_eff", "approx CI"))
  cat(hdr, "\n")
  for (i in seq_len(nrow(d))) {
    line <- sprintf("  %-10s %5d %8s %8s %6.1f %6s %5s", d$exposure[i], d$n[i], fmt(d$sig[i]), fmt(d$sig_modal[i]),
                    d$relative_or[i], ifelse(d$janus[i], "yes", "no"), ifelse(d$direction[i] > 0, "+", "-"))
    if (attr(x, "has_null")) line <- paste(line, sprintf("%8s %7.3f %4.0f/%-4.0f %5s-%-5s", fmt(d$null_median_modal[i]), d$p_perm[i],
                                                        d$n_eff_var[i], d$n_eff_liji[i], fmt(d$ci_lo[i]), fmt(d$ci_hi[i])))
    cat(line, "\n")
  }
  invisible(x)
}

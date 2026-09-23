#' Data-driven definition of defensible specifications
#'
#' Measures, for every level of every analytic choice, the false-positive rate
#' under a permutation null (the share of null specifications with
#' `p < alpha`). A level whose rate exceeds a prespecified threshold is classed
#' as not defensible. Because specifications within one permutation are
#' correlated, the confidence interval of each rate is taken from the variation
#' of the per-replicate rate across permutations, not from a binomial formula.
#'
#' The criterion detects levels that fail to control confounding (with a
#' Freedman-Lane null) and outcome definitions that saturate; it cannot detect
#' over-adjustment, because adjusting for a collider or mediator does not
#' create false positives when there is no association.
#'
#' @param null A [spec_null()] object.
#' @param threshold False-positive rate above which a level is not defensible
#'   (set before looking at the observed results; default 0.10).
#' @param sensitivity Additional thresholds reported for sensitivity.
#' @param alpha Significance level used to count a false positive.
#' @param within Optional named list, e.g. `list(covset = c("socio", "full"))`:
#'   all axes other than the named one are assessed only within these levels
#'   (typically the adequately adjusted covariate sets), so that the behaviour
#'   of a definition or coding is separated from residual confounding. The named
#'   axis itself is assessed across all specifications.
#' @param axes Axis columns to assess (default: all axis columns in the null).
#' @return An object of class `spec_defensible`: a data frame with `axis`,
#'   `level`, `n_specs` (per replicate), `fpr`, `se`, `lo`, `hi`, `defensible`
#'   and one logical column per sensitivity threshold; attributes `threshold`,
#'   `alpha`, `within`, `restrict`.
#' @seealso [spec_keep()] to translate the classification into a logical
#'   vector over a `spec_grid`.
#' @export
spec_defensible <- function(null, threshold = 0.10, sensitivity = c(0.075, 0.15), alpha = 0.05,
                            within = NULL, axes = NULL) {
  stopifnot(inherits(null, "spec_null"))
  P <- null$perm; P$.sig <- P$p < alpha
  if (is.null(axes)) axes <- intersect(null$axis_cols, names(P))
  axes <- axes[vapply(axes, function(a) length(unique(P[[a]])) > 1, TRUE)]
  if (!is.null(within)) {
    stopifnot(is.list(within), length(within) == 1, !is.null(names(within)), names(within) %in% names(P))
    wax <- names(within); wlv <- within[[1]]
    if (!all(wlv %in% P[[wax]])) stop("within levels not found in axis ", wax)
  }
  rows <- list()
  for (a in axes) {
    Pa <- if (!is.null(within) && a != wax) P[P[[wax]] %in% wlv, ] else P
    for (lv in unique(Pa[[a]])) {
      s <- Pa[Pa[[a]] == lv, ]
      r <- tapply(s$.sig, s$rep, mean)                  # per-replicate rate
      fpr <- mean(r); se <- stats::sd(r) / sqrt(length(r))
      row <- data.frame(axis = a, level = lv, n_specs = nrow(s) / length(r), fpr = fpr, se = se,
                        lo = fpr - 1.96 * se, hi = fpr + 1.96 * se, defensible = fpr <= threshold, stringsAsFactors = FALSE)
      for (t in sensitivity) row[[sprintf("defensible_%g", t)]] <- fpr <= t
      rows[[length(rows) + 1]] <- row
    }
  }
  out <- do.call(rbind, rows); rownames(out) <- NULL
  structure(out, class = c("spec_defensible", "data.frame"), threshold = threshold, sensitivity = sensitivity,
            alpha = alpha, within = within, restrict = null$restrict, B = length(null$reps))
}

#' @export
print.spec_defensible <- function(x, digits = 1, ...) {
  thr <- attr(x, "threshold")
  cat(sprintf("Defensible specifications: false-positive rate under the null (%d permutations), threshold %.1f%%\n",
              attr(x, "B"), 100 * thr))
  if (!is.null(attr(x, "within"))) cat("  levels of other axes assessed within", names(attr(x, "within")), "=",
                                       paste(attr(x, "within")[[1]], collapse = ", "), "\n")
  d <- as.data.frame(x)
  d$rate <- sprintf("%5.1f%% (%.1f-%.1f)", 100 * d$fpr, 100 * d$lo, 100 * d$hi)
  d$flag <- ifelse(d$defensible, "", "  <- not defensible")
  for (a in unique(d$axis)) {
    cat("  ", a, "\n", sep = "")
    s <- d[d$axis == a, ]
    for (i in seq_len(nrow(s))) cat(sprintf("     %-24s %s%s\n", s$level[i], s$rate[i], s$flag[i]))
  }
  bad <- d[!d$defensible, ]
  if (nrow(bad)) cat(sprintf("Excluded: %s\n", paste(paste(bad$axis, bad$level, sep = " = "), collapse = "; ")))
  invisible(x)
}

#' Logical vector of defensible specifications in a grid
#'
#' A specification is defensible when every assessed level it uses is
#' defensible (at the chosen threshold) and it lies on the subgrid on which the
#' null was computed.
#' @param def A [spec_defensible()] object.
#' @param grid A `spec_grid`.
#' @param threshold Use the primary threshold (`NULL`, default) or one of the
#'   sensitivity thresholds.
#' @return Logical vector with one element per row of `grid`.
#' @export
spec_keep <- function(def, grid, threshold = NULL) {
  stopifnot(inherits(def, "spec_defensible"), inherits(grid, "spec_grid"))
  col <- if (is.null(threshold)) "defensible" else sprintf("defensible_%g", threshold)
  if (!col %in% names(def)) stop("threshold not available: ", threshold)
  keep <- rep(TRUE, nrow(grid))
  d <- as.data.frame(def)
  for (a in unique(d$axis)) {
    ok <- d$level[d$axis == a & d[[col]]]
    keep <- keep & grid[[a]] %in% ok
  }
  r <- attr(def, "restrict")
  if (!is.null(r)) for (nm in names(r)) keep <- keep & grid[[nm]] %in% r[[nm]]
  keep
}

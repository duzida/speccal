#' Which analytic choice moves the estimate: variance decomposition of a grid
#'
#' Partitions the variance of the specification-level estimates by analysis of
#' variance on the axes of the grid, separately for each exposure. The grid is
#' fully crossed, so sequential sums of squares do not depend on the order of
#' entry. The estimates are deterministic functions of one dataset, so this is
#' a descriptive partition and no F tests are attached.
#'
#' Three refinements guard the interpretation:
#' * `interactions = TRUE` adds all pairwise interactions, so the share of main
#'   effects, pairwise interactions and higher-order terms can be reported.
#' * `equal_cardinality` re-runs the decomposition on random draws of `k`
#'   levels of one axis (all combinations when feasible), to check whether an
#'   axis dominates merely because it has more levels.
#' * `null` computes the same decomposition in every permutation replicate, so
#'   the observed share of an axis can be compared with its distribution when
#'   no association exists. A large share under the null means the axis
#'   transmits confounding and sampling variability, not that it carries
#'   evidence about the association.
#'
#' @param grid A `spec_grid`.
#' @param null Optional [spec_null()] (its `perm` must carry `estimate`).
#' @param scale Column to decompose: `"estimate"` (log odds ratio), `"lrr"`
#'   (log marginal risk ratio, derived from `rr`) or `"rd"`. For the null only
#'   `"estimate"` is available.
#' @param keep Optional logical vector selecting rows of `grid` (and the
#'   matching null specifications).
#' @param interactions Logical; add pairwise interactions.
#' @param equal_cardinality `NULL` or `list(axis = "<axis>", k = <levels>,
#'   max_combos = 500)`.
#' @param axes Axis columns to use (default: all axis columns with more than
#'   one level among the selected rows).
#' @return A list of class `spec_decompose` with `shares` (data frame: exposure,
#'   term, share, df, mean_sq), `main` (wide table of axis shares per exposure),
#'   `null` (per-exposure quantiles of each axis share under the null and the
#'   observed share), `equal_cardinality` (per-exposure median and 5th-95th
#'   percentiles of the shares of the reduced axis and of every other axis, and
#'   the proportion of draws in which the reduced axis is the largest).
#' @export
spec_decompose <- function(grid, null = NULL, scale = "estimate", keep = NULL, interactions = FALSE,
                           equal_cardinality = NULL, axes = NULL) {
  stopifnot(inherits(grid, "spec_grid"))
  G <- as.data.frame(grid); if (!is.null(keep)) G <- G[keep, ]
  if (scale == "lrr" && !"lrr" %in% names(G) && "rr" %in% names(G)) G$lrr <- log(G$rr)
  if (!scale %in% names(G)) stop("column not found: ", scale)
  ax_all <- attr(grid, "axis_cols")
  if (is.null(axes)) axes <- ax_all[vapply(ax_all, function(a) length(unique(G[[a]])) > 1, TRUE)]
  if (length(axes) < 1) stop("no axis varies among the selected specifications")
  f_main <- stats::as.formula(paste("y ~", paste(axes, collapse = " + ")))
  f_int  <- stats::as.formula(paste("y ~ (", paste(axes, collapse = " + "), ")^2"))
  dec <- function(s, f) {
    s$y <- s[[scale]]; s <- s[is.finite(s$y), ]
    a <- stats::anova(stats::lm(f, data = s)); ss <- a[["Sum Sq"]]
    data.frame(term = rownames(a), share = ss / sum(ss), df = a[["Df"]], mean_sq = a[["Mean Sq"]], stringsAsFactors = FALSE)
  }
  exps <- unique(G$exposure)
  shares <- do.call(rbind, lapply(exps, function(e) cbind(exposure = e, dec(G[G$exposure == e, ], if (interactions) f_int else f_main))))
  rownames(shares) <- NULL
  main <- do.call(rbind, lapply(exps, function(e) {
    s <- shares[shares$exposure == e & shares$term %in% axes, ]
    r <- as.data.frame(as.list(stats::setNames(s$share, s$term))); r$residual <- 1 - sum(s$share)
    if (interactions) { i <- shares[shares$exposure == e & grepl(":", shares$term), ]; r$pairwise <- sum(i$share); r$higher_order <- r$residual - r$pairwise }
    cbind(exposure = e, r, stringsAsFactors = FALSE) }))
  out <- list(shares = shares, main = main, axes = axes, scale = scale, interactions = interactions)
  ## null calibration
  if (!is.null(null)) {
    if (scale != "estimate") stop("null calibration is available for scale = 'estimate' only")
    G$.key <- do.call(paste, c(G[ax_all], sep = "\r"))
    P <- null$perm; P$.key <- do.call(paste, c(P[intersect(ax_all, names(P))], sep = "\r"))
    P <- P[P$.key %in% G$.key, ]
    q <- do.call(rbind, lapply(exps, function(e) {
      Pe <- P[P$exposure == e, ]; if (!nrow(Pe)) return(NULL)
      byrep <- split(Pe, Pe$rep)
      sh <- do.call(rbind, lapply(byrep, function(b) { d <- dec(b, f_main); stats::setNames(d$share[match(axes, d$term)], axes) }))
      do.call(rbind, lapply(axes, function(a) data.frame(exposure = e, axis = a, observed = main[[a]][main$exposure == e],
        null_median = stats::median(sh[, a]), null_q05 = unname(stats::quantile(sh[, a], .05)), null_q95 = unname(stats::quantile(sh[, a], .95)),
        p_perm = (sum(sh[, a] >= main[[a]][main$exposure == e] - 1e-12) + 1) / (nrow(sh) + 1), B = nrow(sh), stringsAsFactors = FALSE)))
    }))
    rownames(q) <- NULL; out$null <- q
  }
  ## equal-cardinality sensitivity
  if (!is.null(equal_cardinality)) {
    ec <- equal_cardinality; stopifnot(ec$axis %in% axes, is.numeric(ec$k))
    maxc <- if (is.null(ec$max_combos)) 500 else ec$max_combos
    lv <- unique(G[[ec$axis]]); if (ec$k >= length(lv)) stop("k must be smaller than the number of levels of ", ec$axis)
    ncomb <- choose(length(lv), ec$k)
    combs <- if (ncomb <= maxc) utils::combn(lv, ec$k, simplify = FALSE) else { set.seed(1); replicate(maxc, sample(lv, ec$k), simplify = FALSE) }
    ecr <- do.call(rbind, lapply(exps, function(e) {
      Ge <- G[G$exposure == e, ]
      sh <- do.call(rbind, lapply(combs, function(cc) { d <- dec(Ge[Ge[[ec$axis]] %in% cc, ], f_main); stats::setNames(d$share[match(axes, d$term)], axes) }))
      r <- data.frame(exposure = e, n_draws = length(combs), stringsAsFactors = FALSE)
      for (a in axes) { r[[paste0(a, "_median")]] <- stats::median(sh[, a]); r[[paste0(a, "_q05")]] <- unname(stats::quantile(sh[, a], .05)); r[[paste0(a, "_q95")]] <- unname(stats::quantile(sh[, a], .95)) }
      others <- setdiff(axes, ec$axis)
      r$prop_reduced_axis_largest <- mean(apply(sh, 1, function(v) all(v[ec$axis] > v[others])))
      attr(r, "draws") <- sh; r }))
    out$equal_cardinality <- ecr; out$equal_cardinality_spec <- ec
    out$equal_cardinality_draws <- stats::setNames(lapply(exps, function(e) { Ge <- G[G$exposure == e, ]
      do.call(rbind, lapply(combs, function(cc) { d <- dec(Ge[Ge[[ec$axis]] %in% cc, ], f_main); stats::setNames(d$share[match(axes, d$term)], axes) })) }), exps)
  }
  structure(out, class = "spec_decompose")
}

#' @export
print.spec_decompose <- function(x, ...) {
  cat("Variance decomposition of", x$scale, "by", paste(x$axes, collapse = ", "), if (x$interactions) "(with pairwise interactions)" else "", "\n")
  m <- x$main; num <- setdiff(names(m), "exposure")
  m[num] <- lapply(m[num], function(v) sprintf("%5.1f", 100 * v))
  print(m, row.names = FALSE)
  if (!is.null(x$null)) {
    cat("\nUnder the null (median [5th-95th]) vs observed, share of", x$scale, "variance:\n")
    n <- x$null; n$null <- sprintf("%5.1f [%4.1f-%4.1f]", 100 * n$null_median, 100 * n$null_q05, 100 * n$null_q95); n$observed <- sprintf("%5.1f", 100 * n$observed)
    print(n[, c("exposure", "axis", "observed", "null", "p_perm")], row.names = FALSE)
  }
  if (!is.null(x$equal_cardinality)) {
    ec <- x$equal_cardinality_spec
    cat(sprintf("\nEqual-cardinality sensitivity: %s reduced to %d levels (%d draws); proportion of draws in which it has the largest share:\n", ec$axis, ec$k, x$equal_cardinality$n_draws[1]))
    print(data.frame(exposure = x$equal_cardinality$exposure, prop = sprintf("%.2f", x$equal_cardinality$prop_reduced_axis_largest)), row.names = FALSE)
  }
  invisible(x)
}

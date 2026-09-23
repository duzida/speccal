#' @importFrom ggplot2 .data
NULL

## ggplot2 methods. Colours: one signal colour for significance, one for the Freedman-Lane null, neutral greys otherwise.
.plain <- function(x) { x <- as.data.frame(x); class(x) <- "data.frame"; attributes(x)[setdiff(names(attributes(x)), c("names", "row.names", "class"))] <- NULL; x }
.pal <- list(sig = "#b8322a", fl = "#1f6f8b", simple = "grey62", off = "grey82", ns = "grey30", def = "#b8322a", cov = "#d4912a")

.theme_speccal <- function(base_size = 9) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_text(face = "bold"),
                   legend.title = ggplot2::element_blank(), legend.position = "top", legend.justification = "left",
                   plot.title = ggplot2::element_text(face = "bold", size = base_size + 1))
}

#' Plot a specification curve
#'
#' Odds ratios ranked within exposure, with confidence intervals, coloured by
#' significance and, when `keep` is supplied, by whether the specification is
#' defensible.
#' @param x A `spec_grid`.
#' @param keep Optional logical vector (e.g. from [spec_keep()]).
#' @param alpha Significance level.
#' @param ylim Limits of the odds-ratio axis (intervals are truncated there).
#' @param nrow Rows of facets.
#' @param ... Ignored.
#' @return A ggplot object.
#' @export
plot.spec_grid <- function(x, keep = NULL, alpha = 0.05, ylim = c(0.2, 10), nrow = 2, ...) {
  d <- .plain(x)
  d$keep <- if (is.null(keep)) TRUE else keep
  d$cls <- ifelse(!d$keep, "Not defensible", ifelse(d$p < alpha, sprintf("Defensible, P < %g", alpha), "Defensible, not significant"))
  if (is.null(keep)) d$cls <- ifelse(d$p < alpha, sprintf("P < %g", alpha), "Not significant")
  d <- d[order(d$exposure, d$or), ]; d$i <- stats::ave(seq_len(nrow(d)), d$exposure, FUN = seq_along)
  d$lo <- pmax(d$lo, ylim[1]); d$hi <- pmin(d$hi, ylim[2])
  cols <- if (is.null(keep)) stats::setNames(c(.pal$sig, .pal$ns), c(sprintf("P < %g", alpha), "Not significant")) else
    stats::setNames(c(.pal$sig, .pal$ns, .pal$off), c(sprintf("Defensible, P < %g", alpha), "Defensible, not significant", "Not defensible"))
  d$exposure <- factor(d$exposure, levels = unique(x$exposure))
  ggplot2::ggplot(d, ggplot2::aes(.data$i, .data$or)) +
    ggplot2::geom_hline(yintercept = 1, linetype = 2, linewidth = .3, colour = "grey45") +
    ggplot2::geom_linerange(ggplot2::aes(ymin = .data$lo, ymax = .data$hi, colour = .data$cls), alpha = .25, linewidth = .3) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$cls), size = .5) +
    ggplot2::facet_wrap(~exposure, nrow = nrow) +
    ggplot2::scale_colour_manual(values = cols, breaks = names(cols)) +
    ggplot2::scale_y_log10() + ggplot2::coord_cartesian(ylim = ylim) +
    ggplot2::scale_x_continuous(breaks = NULL) +
    ggplot2::labs(x = "Specification (ranked by odds ratio)", y = "Odds ratio (95% CI)") +
    .theme_speccal() + ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2, alpha = 1)))
}

#' Plot false-positive rates by level
#' @param x A `spec_defensible` object.
#' @param ... Ignored.
#' @return A ggplot object.
#' @export
plot.spec_defensible <- function(x, ...) {
  d <- .plain(x); thr <- attr(x, "threshold"); sens <- attr(x, "sensitivity")
  d$status <- ifelse(d$defensible, sprintf("At or below %g%%", 100 * thr), sprintf("Above %g%% (not defensible)", 100 * thr))
  d$level <- factor(d$level, levels = rev(unique(d$level))); d$axis <- factor(d$axis, levels = unique(d$axis))
  ggplot2::ggplot(d, ggplot2::aes(100 * .data$fpr, .data$level)) +
    ggplot2::geom_vline(xintercept = 100 * attr(x, "alpha"), linetype = 3, colour = "grey55") +
    ggplot2::geom_vline(xintercept = 100 * sens, linetype = 2, colour = "grey80") +
    ggplot2::geom_vline(xintercept = 100 * thr, colour = "grey35") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = 100 * .data$lo, xmax = 100 * .data$hi), width = 0, orientation = "y", colour = "grey30") +
    ggplot2::geom_point(ggplot2::aes(colour = .data$status), size = 2) +
    ggplot2::scale_colour_manual(values = stats::setNames(c(.pal$sig, "grey20"), c(sprintf("Above %g%% (not defensible)", 100 * thr), sprintf("At or below %g%%", 100 * thr)))) +
    ggplot2::facet_grid(axis ~ ., scales = "free_y", space = "free_y") +
    ggplot2::labs(x = "False-positive rate under the null (%)", y = NULL,
                  title = sprintf("Dotted: nominal %g%%; solid: threshold %g%%; dashed: sensitivity thresholds", 100 * attr(x, "alpha"), 100 * thr)) +
    .theme_speccal() + ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0))
}

#' Plot observed shares of significant specifications against the null
#' @param x A `spec_summary` computed with a null.
#' @param null The [spec_null()] used, to draw the null distribution.
#' @param keep The same `keep` used in [spec_summary()].
#' @param grid The `spec_grid` (needed with `keep`).
#' @param ... Ignored.
#' @return A ggplot object.
#' @export
plot.spec_summary <- function(x, null = NULL, keep = NULL, grid = NULL, ...) {
  d <- .plain(x); d$exposure <- factor(d$exposure, levels = rev(d$exposure))
  p <- ggplot2::ggplot() + ggplot2::labs(x = sprintf("Specifications significant in the modal direction (%%)"), y = NULL) + .theme_speccal()
  if (!is.null(null) && "p_perm" %in% names(d)) {
    P <- null$perm; a <- attr(x, "alpha")
    if (!is.null(keep)) { stopifnot(!is.null(grid)); ax <- attr(grid, "axis_cols"); G <- as.data.frame(grid)[keep, ]
      key <- do.call(paste, c(G[ax], sep = "\r")); P <- P[do.call(paste, c(P[ax], sep = "\r")) %in% key, ] }
    dir <- stats::setNames(d$direction, as.character(d$exposure))
    nb <- stats::aggregate(list(h = P$p < a & sign(P$estimate) == dir[P$exposure]), list(rep = P$rep, exposure = P$exposure), mean)
    nb$exposure <- factor(nb$exposure, levels = levels(d$exposure))
    p <- p + ggplot2::geom_boxplot(data = nb, ggplot2::aes(100 * .data$h, .data$exposure, colour = "Null"), width = .5, outlier.size = .4, fill = "grey92") +
      ggplot2::geom_text(data = d, ggplot2::aes(x = 104, y = .data$exposure, label = sprintf("P = %.3f", .data$p_perm)), hjust = 0, size = 2.8, colour = "grey20") +
      ggplot2::scale_x_continuous(limits = c(0, 125), breaks = seq(0, 100, 25))
  }
  p + ggplot2::geom_point(data = d, ggplot2::aes(100 * .data$sig_modal, .data$exposure, colour = "Observed"), size = 2.2) +
    ggplot2::scale_colour_manual(values = c(Observed = .pal$sig, Null = .pal$fl))
}

#' Plot a variance decomposition
#' @param x A `spec_decompose` object.
#' @param ... Ignored.
#' @return A ggplot object (patchwork when a null is present).
#' @export
plot.spec_decompose <- function(x, ...) {
  m <- x$main; ax <- x$axes
  long <- do.call(rbind, lapply(c(ax, "residual"), function(a) data.frame(exposure = m$exposure, term = a, share = m[[a]])))
  long$term <- factor(long$term, levels = rev(c(ax, "residual"))); long$exposure <- factor(long$exposure, levels = m$exposure)
  cols <- stats::setNames(c("grey88", grDevices::grey.colors(length(ax) - 2, 0.35, 0.7), .pal$cov, .pal$def)[seq_len(length(ax) + 1)], rev(c(ax, "residual")))
  pA <- ggplot2::ggplot(long, ggplot2::aes(.data$exposure, 100 * .data$share, fill = .data$term)) +
    ggplot2::geom_col(width = .7, colour = "white", linewidth = .2) +
    ggplot2::scale_fill_manual(values = cols, breaks = c(ax, "residual")) +
    ggplot2::labs(x = NULL, y = sprintf("Share of variance in %s (%%)", x$scale)) + .theme_speccal()
  if (is.null(x$null)) return(pA)
  n <- x$null; n$exposure <- factor(n$exposure, levels = m$exposure)
  pB <- ggplot2::ggplot(n, ggplot2::aes(y = .data$exposure)) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = 100 * .data$null_q05, xmax = 100 * .data$null_q95, colour = "Null, 5th-95th percentile"), width = .3, orientation = "y") +
    ggplot2::geom_point(ggplot2::aes(x = 100 * .data$observed, colour = "Observed"), size = 2) +
    ggplot2::facet_wrap(~axis, nrow = 1) +
    ggplot2::scale_colour_manual(values = c(Observed = .pal$sig, "Null, 5th-95th percentile" = .pal$fl)) +
    ggplot2::labs(x = "Share of variance (%)", y = NULL) + .theme_speccal()
  patchwork::wrap_plots(pA, pB, ncol = 1) + patchwork::plot_annotation(tag_levels = "a")
}

#' Plot the bootstrap distribution of a breakpoint
#' @param x A `breakpoint_boot` object.
#' @param ... Ignored.
#' @return A ggplot object.
#' @export
plot.breakpoint_boot <- function(x, ...) {
  g <- x$grid; cnt <- data.frame(bp = g, n = vapply(g, function(v) sum(abs(x$boot$bp - v) < 1e-9 * v, na.rm = TRUE), 1L))
  cnt$edge <- ifelse(seq_along(g) %in% c(1, length(g)), "Bound of the search grid", "Interior candidate")
  s <- x$summary
  ggplot2::ggplot(cnt, ggplot2::aes(.data$bp, .data$n)) +
    ggplot2::geom_segment(ggplot2::aes(xend = .data$bp, y = 0, yend = .data$n, colour = .data$edge), linewidth = 1) +
    ggplot2::annotate("segment", x = s$lo, xend = s$hi, y = max(cnt$n) * 1.08, yend = max(cnt$n) * 1.08, colour = "grey15") +
    ggplot2::annotate("point", x = s$estimate, y = max(cnt$n) * 1.08, size = 2, colour = "grey15") +
    ggplot2::scale_colour_manual(values = c("Interior candidate" = "grey45", "Bound of the search grid" = .pal$sig)) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(x = sprintf("Breakpoint in %s (log scale)", x$exposure), y = sprintf("Bootstrap replicates (of %d)", s$B),
                  title = sprintf("%s bootstrap: %.1f%% at lower bound, %.1f%% at upper bound; estimate and 95%% interval above",
                                  x$resampling, 100 * s$at_lower_bound, 100 * s$at_upper_bound)) +
    .theme_speccal()
}

#' Coerce a table of specification-level results into a `spec_grid`
#'
#' Use this for results produced outside speccal (for example by `specr` or a
#' hand-written loop). The table needs one row per specification, columns
#' identifying the axis levels, and an estimate on the log scale with its
#' standard error and P value.
#'
#' @param df A data frame.
#' @param axes Character vector of column names that identify the axes.
#' @param exposure Column name identifying the exposure (or `NULL` for a single
#'   exposure).
#' @param estimate,se,p Column names of the log-scale estimate, its standard
#'   error and its P value.
#' @param refit Optional function `refit(newdata)` returning a `spec_grid` with
#'   the same rows for new data; required by [spec_null()].
#' @return A `spec_grid`.
#' @export
as_spec_grid <- function(df, axes, exposure = NULL, estimate = "estimate", se = "se", p = "p", refit = NULL) {
  stopifnot(is.data.frame(df), all(c(axes, exposure, estimate, se, p) %in% names(df)))
  out <- df
  if (is.null(exposure)) out$exposure <- "exposure" else if (exposure != "exposure") out$exposure <- df[[exposure]]
  if (estimate != "estimate") out$estimate <- df[[estimate]]
  if (se != "se") out$se <- df[[se]]
  if (p != "p") out$p <- df[[p]]
  if (!"or" %in% names(out)) out$or <- exp(out$estimate)
  for (nm in axes) out[[nm]] <- as.character(out[[nm]])
  structure(out, class = c("spec_grid", "data.frame"), axes = NULL, exposure = unique(out$exposure),
            axis_cols = axes, refit = refit)
}

#' @export
print.spec_grid <- function(x, ...) {
  ax <- attr(x, "axis_cols"); ex <- unique(x$exposure)
  cat("spec_grid:", nrow(x), "specifications,", length(ex), "exposure(s),", length(ax), "axes:", paste(ax, collapse = ", "), "\n")
  if ("sep_flag" %in% names(x)) cat("  converged:", sum(x$converged, na.rm = TRUE), "/", nrow(x),
                                   "| quasi-separation flags:", sum(x$sep_flag, na.rm = TRUE), "\n")
  if (!is.null(attr(x, "refit"))) cat("  refit closure available (spec_null can be run)\n")
  invisible(x)
}

#' Axis columns of a `spec_grid`
#' @param grid A `spec_grid`.
#' @export
axis_cols <- function(grid) attr(grid, "axis_cols")

## keep attributes through base subsetting
#' @export
`[.spec_grid` <- function(x, i, j, ..., drop = FALSE) {
  at <- attributes(x)[c("axes", "exposure", "axis_cols", "refit")]
  out <- NextMethod()
  if (is.data.frame(out)) { for (nm in names(at)) attr(out, nm) <- at[[nm]]; class(out) <- c("spec_grid", "data.frame") }
  out
}

#' Define the axes of a specification grid
#'
#' An axis is one analytic choice with two or more levels. The grid is the full
#' cross of all axes. Every level must be something a reasonable analyst could
#' have chosen; the grid does not endorse levels, it measures their consequence.
#'
#' @param outcome Named character vector or list: level name -> outcome variable
#'   name in `data` (binary 0/1).
#' @param covariates Named list: level name -> character vector of covariate
#'   names (`character(0)` for a crude model).
#' @param coding Exposure parameterisation. Either a character vector of built-in
#'   codings (`"per_sd"`, `"log2"`, `"q4_vs_q1"`, `"raw"`) or a named list of
#'   `list(fun = function(v) ..., contrast = "shift" | "binary")`. `fun` maps the
#'   raw exposure to the analysed variable (computed on the full sample before
#'   any subsetting); `contrast` tells the marginal-effect routine whether the
#'   effect is a one-unit shift or a 0/1 contrast.
#' @param weighting Character vector with any of `"design"` (survey-weighted,
#'   design-based variance) and `"unweighted"` (ordinary logistic regression).
#' @param ... Additional subset axes, each a named list of quoted expressions
#'   (or logical vectors) evaluated in `data`, e.g.
#'   `sample = list(all = TRUE, age40 = quote(Age >= 40))`,
#'   `missing = list(imputed = TRUE, complete = quote(complete_case == 1))`.
#'   Subsets are applied as survey domains (`subset()` on the design), never by
#'   filtering the data before the design is built.
#' @return An object of class `spec_axes`.
#' @examples
#' ax <- spec_axes(
#'   outcome    = c(any = "y_any", severe = "y_severe"),
#'   covariates = list(crude = character(0), demo = c("Age", "Gender")),
#'   coding     = c("per_sd", "log2"),
#'   sample     = list(all = TRUE, age40 = quote(Age >= 40)))
#' ax
#' @export
spec_axes <- function(outcome, covariates, coding = c("per_sd", "log2", "q4_vs_q1"),
                      weighting = c("design", "unweighted"), ...) {
  outcome <- unlist(outcome)
  stopifnot(is.character(outcome), length(outcome) >= 1, !is.null(names(outcome)), all(nzchar(names(outcome))))
  stopifnot(is.list(covariates), length(covariates) >= 1, !is.null(names(covariates)))
  covariates <- lapply(covariates, function(v) if (is.null(v)) character(0) else as.character(v))
  coding <- .resolve_coding(coding)
  weighting <- match.arg(weighting, c("design", "unweighted"), several.ok = TRUE)
  extra <- list(...)
  if (length(extra)) {
    stopifnot(!is.null(names(extra)), all(nzchar(names(extra))))
    for (nm in names(extra)) {
      ax <- extra[[nm]]
      stopifnot(is.list(ax), !is.null(names(ax)), all(nzchar(names(ax))))
      if (nm %in% c("outcome", "covset", "coding", "weight", "exposure"))
        stop("axis name '", nm, "' is reserved")
    }
  }
  structure(list(outcome = outcome, covariates = covariates, coding = coding,
                 subsets = extra, weighting = weighting), class = "spec_axes")
}

.builtin_codings <- list(
  per_sd   = list(fun = function(v) as.numeric(scale(v)), contrast = "shift"),
  log2     = list(fun = function(v) log2(v), contrast = "shift"),
  raw      = list(fun = function(v) as.numeric(v), contrast = "shift"),
  q4_vs_q1 = list(fun = function(v) {
    q <- cut(v, stats::quantile(v, 0:4 / 4, na.rm = TRUE), include.lowest = TRUE, labels = paste0("Q", 1:4))
    ifelse(q == "Q4", 1, ifelse(q == "Q1", 0, NA_real_))
  }, contrast = "binary"))

.resolve_coding <- function(coding) {
  if (is.character(coding)) {
    bad <- setdiff(coding, names(.builtin_codings))
    if (length(bad)) stop("unknown built-in coding: ", paste(bad, collapse = ", "))
    return(.builtin_codings[coding])
  }
  stopifnot(is.list(coding), !is.null(names(coding)))
  lapply(coding, function(cd) {
    stopifnot(is.function(cd$fun)); cd$contrast <- match.arg(cd$contrast, c("shift", "binary")); cd })
}

#' @export
print.spec_axes <- function(x, ...) {
  n <- c(outcome = length(x$outcome), covset = length(x$covariates), coding = length(x$coding),
         vapply(x$subsets, length, 1L), weight = length(x$weighting))
  cat("Specification axes:", length(n), "axes,", prod(n), "specifications per exposure\n")
  for (i in seq_along(n)) cat(sprintf("  %-12s %d level%s\n", names(n)[i], n[i], if (n[i] > 1) "s" else ""))
  invisible(x)
}

#' Number of specifications per exposure
#' @param axes A `spec_axes` object.
#' @export
n_specs <- function(axes) {
  prod(c(length(axes$outcome), length(axes$covariates), length(axes$coding),
         vapply(axes$subsets, length, 1L), length(axes$weighting)))
}

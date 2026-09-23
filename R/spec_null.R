#' Permutation null distribution for a specification grid
#'
#' Refits the grid on permuted data to obtain, for every specification, the
#' distribution of estimates and P values under a null of no exposure-outcome
#' association. Two schemes are available:
#'
#' * `"fl"` (default): within-stratum Freedman-Lane permutation. Each
#'   (transformed) exposure is regressed on the covariates in `covariates`; the
#'   fitted part is kept and only the residuals are permuted within strata. This
#'   removes the exposure-outcome association while preserving the association
#'   of the exposure with the covariates, so the null is "no association given
#'   the covariates", i.e. confounding is intact.
#' * `"simple"`: within-stratum permutation of the exposure itself, which also
#'   breaks the exposure-covariate association and therefore tests a null in
#'   which confounding does not exist. It is anticonservative when the exposure
#'   is confounded; it is provided for comparison.
#'
#' All exposures in the grid are permuted jointly with the same permutation, so
#' their correlation structure is preserved. Every replicate is saved
#' specification by specification, which is what the calibration functions
#' need.
#'
#' @param grid A `spec_grid` produced by [spec_fit()] (or by [as_spec_grid()]
#'   with a `refit` function and `data`).
#' @param scheme `"fl"` or `"simple"`.
#' @param B Number of permutations; ignored when `reps` is given.
#' @param reps Integer vector of replicate indices to run (for resuming or
#'   for verification against archived replicates). Defaults to `1:B`.
#' @param stratum Name of the stratum variable in the data; permutation is
#'   within its levels. `NULL` permutes across the whole sample.
#' @param covariates Character vector of covariate names for the Freedman-Lane
#'   regression (typically the fullest covariate set). Required for `"fl"`.
#' @param transform,inverse Functions applied to the exposure before the
#'   Freedman-Lane regression and to map back afterwards (default `log`/`exp`,
#'   appropriate for positive, right-skewed markers).
#' @param restrict Named list restricting subset and weighting axes for the
#'   null, e.g. `list(missing = "imputed", weight = "design")`; the null is
#'   then computed on that subgrid only. `NULL` uses the full grid.
#' @param seed Integer; replicate `b` uses `set.seed(seed + b)`.
#' @param dir Optional directory; each replicate is written to
#'   `sprintf("rep_%04d.csv", b)` there and skipped if the file exists, so long
#'   jobs can be resumed.
#' @param cores Cores for `parallel::mclapply` over replicates.
#' @param data Data frame; required only for grids from [as_spec_grid()].
#' @param keep Columns of the refitted grid to store per replicate.
#' @return An object of class `spec_null`: a list with `perm` (data frame:
#'   `rep`, `exposure`, axis columns, kept columns), `scheme`, `reps`, `seed`,
#'   `stratum`, `restrict`, `axis_cols`.
#' @export
spec_null <- function(grid, scheme = c("fl", "simple"), B = 300L, reps = NULL, stratum = NULL,
                      covariates = NULL, transform = log, inverse = exp, restrict = NULL,
                      seed = 1L, dir = NULL, cores = 1L, data = NULL,
                      keep = c("estimate", "p", "rr", "p_rr", "sep_flag")) {
  scheme <- match.arg(scheme)
  refit <- attr(grid, "refit"); if (is.null(refit)) stop("grid has no refit function; use spec_fit() or supply refit to as_spec_grid()")
  ctx <- attr(grid, "context"); if (is.null(data)) data <- ctx$data
  if (is.null(data)) stop("data is required (not stored in this grid)")
  exposures <- unique(grid$exposure)
  for (ex in exposures) if (is.null(data[[ex]])) stop("exposure '", ex, "' not in data")
  if (scheme == "fl" && is.null(covariates)) stop("covariates are required for the Freedman-Lane scheme")
  if (is.null(reps)) reps <- seq_len(B)
  axes <- attr(grid, "axes")
  axes_null <- if (!is.null(restrict)) restrict_axes(axes, restrict) else axes
  n <- nrow(data)
  strata_idx <- if (is.null(stratum)) list(seq_len(n)) else split(seq_len(n), data[[stratum]])
  if (scheme == "fl") {
    Xc <- stats::model.matrix(stats::as.formula(paste("~", paste(covariates, collapse = "+"))), data = data)
    if (nrow(Xc) != n) stop("covariates contain missing values; the Freedman-Lane regression needs complete covariates")
    L <- sapply(exposures, function(ex) transform(data[[ex]]))
    fitL <- Xc %*% qr.coef(qr(Xc), L); resL <- L - fitL
  }
  perm_within <- function() { idx <- seq_len(n); for (s in strata_idx) idx[s] <- if (length(s) > 1) sample(s) else s; idx }
  if (!is.null(dir)) dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  one <- function(b) {
    fp <- if (!is.null(dir)) file.path(dir, sprintf("rep_%04d.csv", b)) else NULL
    if (!is.null(fp) && file.exists(fp)) return(utils::read.csv(fp, stringsAsFactors = FALSE))
    set.seed(seed + b); idx <- perm_within(); d <- data
    if (scheme == "simple") { for (ex in exposures) d[[ex]] <- data[[ex]][idx] }
    else { Lp <- fitL + resL[idx, , drop = FALSE]; for (j in seq_along(exposures)) d[[exposures[j]]] <- inverse(Lp[, j]) }
    g <- refit(d, axes_new = axes_null, cores = 1L)
    out <- as.data.frame(g)[, c("exposure", attr(g, "axis_cols"), intersect(keep, names(g)))]
    out <- cbind(rep = b, out)
    if (!is.null(fp)) utils::write.csv(out, fp, row.names = FALSE)
    out
  }
  res <- if (cores > 1) parallel::mclapply(reps, function(b) try(one(b), silent = TRUE), mc.cores = cores, mc.preschedule = FALSE)
         else lapply(reps, function(b) try(one(b), silent = TRUE))
  bad <- vapply(res, function(r) inherits(r, "try-error"), TRUE)
  if (any(bad)) stop("replicates failed: ", paste(reps[bad], collapse = ", "), "\n", res[[which(bad)[1]]])
  perm <- do.call(rbind, res); rownames(perm) <- NULL
  structure(list(perm = perm, scheme = scheme, reps = reps, seed = seed, stratum = stratum,
                 restrict = restrict, axis_cols = attr(grid, "axis_cols"), exposures = exposures),
            class = "spec_null")
}

#' Restrict subset and weighting axes to chosen levels
#' @param axes A `spec_axes` object.
#' @param restrict Named list: axis name -> level name(s) to keep. Applies to
#'   subset axes and to `weight`; `outcome`, `covset` and `coding` may also be
#'   restricted by level name.
#' @return A `spec_axes` object.
#' @export
restrict_axes <- function(axes, restrict) {
  stopifnot(inherits(axes, "spec_axes"), is.list(restrict), !is.null(names(restrict)))
  for (nm in names(restrict)) {
    lv <- restrict[[nm]]
    if (nm == "weight") { bad <- setdiff(lv, axes$weighting); if (length(bad)) stop("unknown weighting level: ", bad); axes$weighting <- lv }
    else if (nm == "outcome") axes$outcome <- axes$outcome[lv]
    else if (nm == "covset") axes$covariates <- axes$covariates[lv]
    else if (nm == "coding") axes$coding <- axes$coding[lv]
    else if (nm %in% names(axes$subsets)) axes$subsets[[nm]] <- axes$subsets[[nm]][lv]
    else stop("unknown axis: ", nm)
    if (any(is.na(names(unlist(lapply(list(axes$outcome, axes$covariates, axes$coding), names)))))) stop("unknown level in axis ", nm)
  }
  axes
}

#' @export
print.spec_null <- function(x, ...) {
  cat("spec_null:", x$scheme, "permutation,", length(x$reps), "replicates,",
      nrow(x$perm) / length(x$reps), "specifications per replicate,", length(x$exposures), "exposure(s)\n")
  if (!is.null(x$stratum)) cat("  permuted within:", x$stratum, "\n")
  if (!is.null(x$restrict)) cat("  restricted to:", paste(names(x$restrict), unlist(x$restrict), sep = " = ", collapse = "; "), "\n")
  invisible(x)
}

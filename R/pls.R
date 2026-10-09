#' Supervised dimension reduction using censored quantile covariance
#'
#' @param X Numeric predictor matrix (observations by predictors).
#' @param Y Numeric response vector, on the scale used for the quantile fit.
#' @param ncomp Maximum number of components.
#' @param weights Finite, nonnegative inverse censoring probability weights.
#'   These may already include the event indicator. Entries with
#'   `delta == 0` are set to zero before constructing the IPCW quantile
#'   and the supervised covariance.
#' @param delta Event indicators, with 1 denoting an observed event.
#' @param lambda Nonnegative penalty parameter.
#' @param P Symmetric positive semidefinite penalty matrix. Required when
#'   `lambda > 0`; the unused vector default from the original script is replaced
#'   by an explicit matrix requirement.
#' @param norm,maxit Retained for source compatibility. Components are always
#'   normalized; this direct covariance algorithm does not iterate.
#' @param tau Target quantile in (0, 1).
#' @param scale.X Whether to standardize predictors as a coordinate
#'   reparameterization. The constraint is defined in the original predictor
#'   coordinates, so standardization preserves the fitted component space.
#' @param grid Retained for source compatibility. The IPCW marginal quantile
#'   is obtained from an empirical CDF and does not use a regression quantile grid.
#' @details The extraction weights in the original residual-predictor
#'   coordinates maximize squared censored quantile covariance subject to
#'   `w' (I + lambda * P) w = 1`. With `D = diag(Xscale)`, the equivalent
#'   standardized metric is `D^(-1) (I + lambda * P) D^(-1)`.
#'   The marginal threshold is the generalized inverse of the normalized
#'   IPCW empirical CDF, without interpolation. A common rescaling of the
#'   IPCW weights is used internally for numerical stability; it changes
#'   neither that quantile nor the normalized extraction directions.
#'   Extraction stops when the residual covariance is numerically zero relative
#'   to the original predictor and weighted-indicator moments.
#'   Supplied weights must be estimated using the original event indicators;
#'   masking censored weights here cannot undo changes to the fitted
#'   censoring distribution. When using [get_weights()], set
#'   `force_last_event = FALSE`.
#' @return A list containing raw and standardized projection matrices `W` and
#'   `W_scaled`, scores `U`, centering and scaling vectors, the scaled penalty,
#'   residual-coordinate directions, and the number of extracted components.
#'   `directions_raw` contains the residual-coordinate directions on the original
#'   predictor scale, `metric` is their `I + lambda * P` metric,
#'   `metric_scaled` is its standardized-coordinate counterpart,
#'   and `ipcw_weights` records the event weights used by the procedure.
#'   In particular, `sweep(X, 2, Xmean, "-") %*% W` equals `U`.
#' @export
#' @examples
#' d <- simulate_sfcqr_data(n = 100, censoring = 0, seed = 1)
#' x <- d$Z %*% d$Bxy
#' p <- plsfit_cqcov(x, log(d$Y), 2, rep(1, 100), d$delta)
#' stopifnot(isTRUE(all.equal(unname(p$U),
#'   unname(sweep(x, 2, p$Xmean, "-") %*% p$W))))
plsfit_cqcov <- function(X, Y, ncomp = 3, weights, delta,
                        lambda = 0, P = NULL, norm = TRUE, maxit = 100,
                        tau = 0.5, scale.X = TRUE, grid = NULL) {
  X <- .numeric_matrix(X, "X")
  response <- .validate_response(Y, delta, positive = FALSE)
  Y <- response$Y
  delta <- response$delta
  if (nrow(X) != length(Y)) stop("nrow(X) must equal length(Y).")
  tau <- .validate_tau(tau)
  if (nrow(X) < 2L) stop("At least two observations are required.")
  .scalar_number(lambda, "lambda", lower = 0)
  .scalar_number(ncomp, "ncomp", lower = 1, upper = .Machine$integer.max)
  if (ncomp != as.integer(ncomp)) stop("ncomp must be an integer.")
  if (!is.logical(scale.X) || length(scale.X) != 1L || is.na(scale.X)) {
    stop("scale.X must be TRUE or FALSE.")
  }
  weights <- .numeric_vector(weights, "weights")
  if (length(weights) != nrow(X) || any(weights < 0)) {
    stop("weights must be nonnegative and have one entry per observation.")
  }
  # Enforce the original observed event indicators in both IPCW moments.
  weights[delta == 0] <- 0
  if (!any(weights > 0)) {
    stop("At least one observed event must have positive IPCW weight.")
  }
  weights_stable <- weights / max(weights)
  tol <- sqrt(.Machine$double.eps)
  npred <- ncol(X)
  Xmean <- colMeans(X)
  X0_raw <- sweep(X, 2, Xmean, "-")
  Xscale <- if (scale.X) apply(X0_raw, 2, stats::sd) else rep(1, npred)
  Xscale[!is.finite(Xscale) | Xscale < tol] <- 1
  X0 <- sweep(X0_raw, 2, Xscale, "/")
  effective_rank <- qr(X0, tol = tol)$rank
  if (effective_rank == 0L) stop("X has no nonconstant predictor directions.")
  ncomp <- min(as.integer(ncomp), effective_rank, nrow(X) - 1L)

  # P and the Euclidean part of the constraint are defined on the raw scale.
  penalty_raw <- matrix(0, npred, npred)
  if (lambda > 0) penalty_raw <- .validate_penalty(P, npred)
  penaltyMatrix <- penalty_raw / outer(Xscale, Xscale)
  A_lambda <- diag(npred) + lambda * penalty_raw
  if (any(!is.finite(A_lambda))) {
    stop("Penalized metric is nonfinite; reduce lambda or rescale P.")
  }
  A_factor <- tryCatch(chol(A_lambda), error = function(e) {
    stop("The penalized metric is not numerically positive definite; ",
         "check P or reduce lambda.", call. = FALSE)
  })
  metric_scaled <- A_lambda / outer(Xscale, Xscale)

  Y0_tau <- .ipcw_marginal_quantile(Y, weights_stable, tau)
  # Detect cancellation at machine precision relative to the original moments.
  # This tolerance depends on data scale, rather than the penalized score scale.
  marker <- weights_stable * (tau - as.numeric(Y <= Y0_tau))
  marker_rms <- sqrt(mean((marker - mean(marker))^2))
  raw_column_scale <- apply(abs(X0_raw), 2, max)
  raw_column_scale[raw_column_scale == 0] <- 1
  predictor_rms <- raw_column_scale * sqrt(colMeans(
    sweep(X0_raw, 2, raw_column_scale, "/")^2))
  q_tolerance <- 64 * .Machine$double.eps * predictor_rms * marker_rms
  W_scaled <- W_raw <- directions <- matrix(0, npred, ncomp)
  U <- matrix(0, nrow(X), ncomp)
  Xwork <- X0
  R_deflation <- diag(npred)
  used <- 0L
  for (a in seq_len(ncomp)) {
    # cq_cov uses 1/(n-1); rescale to the manuscript's 1/n convention.
    q_scaled <- as.numeric(cq_cov(Xwork, Y, Y0_tau, weights_stable, tau)) *
      ((nrow(X) - 1) / nrow(X))
    q_raw <- q_scaled * Xscale
    if (any(!is.finite(q_raw))) {
      stop("Quantile covariance is nonfinite; rescale the predictor values.")
    }
    q_scale <- max(abs(q_raw))
    if (all(abs(q_raw) <= q_tolerance)) break
    # A common positive scaling of q cancels in the normalized optimizer.
    q_unit <- q_raw / q_scale
    a_raw <- as.numeric(backsolve(A_factor,
      forwardsolve(t(A_factor), q_unit)))
    metric_norm_sq <- sum(q_unit * a_raw)
    if (!is.finite(metric_norm_sq) || metric_norm_sq <= 0) {
      stop("Unable to normalize the direction in the penalized metric.")
    }
    w_raw <- a_raw / sqrt(metric_norm_sq)
    w <- w_raw * Xscale
    score <- Xwork %*% w
    score_scale <- max(abs(score))
    if (!is.finite(score_scale)) stop("Extracted scores are nonfinite.")
    if (score_scale == 0) break
    # Deflate using unit-scaled scores: equivalent, without squaring tiny scores.
    score_unit <- score / score_scale
    den <- sum(score_unit^2)
    projection <- as.numeric(R_deflation %*% w)
    loading_unit <- crossprod(score_unit, Xwork) / den
    Xwork <- Xwork - score_unit %*% loading_unit
    R_deflation <- R_deflation -
      (projection / score_scale) %*% loading_unit
    W_scaled[, a] <- projection
    W_raw[, a] <- projection / Xscale
    directions[, a] <- w
    U[, a] <- score
    used <- a
  }
  if (used == 0L) stop("No nonzero quantile covariance components were found.")
  keep <- seq_len(used)
  W_raw <- W_raw[, keep, drop = FALSE]
  W_scaled <- W_scaled[, keep, drop = FALSE]
  directions <- directions[, keep, drop = FALSE]
  U <- U[, keep, drop = FALSE]
  colnames(W_raw) <- colnames(W_scaled) <- colnames(directions) <-
    colnames(U) <- paste0("Component", keep)
  list(W = W_raw, W_scaled = W_scaled,
       U = U, Xmean = Xmean, Xscale = Xscale,
       penaltyMatrix = penaltyMatrix,
       directions = directions,
       directions_raw = sweep(directions, 1, Xscale, "/"),
       metric = A_lambda, metric_scaled = metric_scaled,
       ipcw_weights = weights, ncomp = used,
       tau = tau, Y_quantile = Y0_tau)
}


# Generalized inverse of the normalized IPCW empirical CDF (no interpolation).
.ipcw_marginal_quantile <- function(Y, weights, tau) {
  positive <- which(weights > 0)
  if (!length(positive)) stop("The IPCW empirical CDF has zero total weight.")
  order_events <- positive[order(Y[positive])]
  ordered_weights <- weights[order_events] / max(weights[positive])
  cumulative_weight <- cumsum(ordered_weights)
  cdf <- cumulative_weight / cumulative_weight[length(cumulative_weight)]
  cdf[length(cdf)] <- 1
  Y[order_events[which(cdf >= tau)[1L]]]
}

.validate_penalty <- function(P, size) {
  if (is.null(P) || !is.matrix(P) || !is.numeric(P) ||
      !identical(dim(P), c(as.integer(size), as.integer(size))) ||
      any(!is.finite(P))) {
    stop("P/penalty must be a finite square matrix matching the basis dimension.")
  }
  tol <- sqrt(.Machine$double.eps) * max(1, max(abs(P)))
  if (max(abs(P - t(P))) > tol) stop("P/penalty must be symmetric.")
  P <- (P + t(P)) / 2
  if (min(eigen(P, symmetric = TRUE, only.values = TRUE)$values) < -tol) {
    stop("P/penalty must be positive semidefinite.")
  }
  P
}

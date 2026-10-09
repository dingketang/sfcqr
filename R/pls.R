#' Supervised dimension reduction using censored quantile covariance
#'
#' @param X Numeric predictor matrix (observations by predictors).
#' @param Y Numeric response vector, on the scale used for the quantile fit.
#' @param ncomp Maximum number of components.
#' @param weights Nonnegative inverse censoring probability weights.
#' @param delta Event indicators, with 1 denoting an observed event.
#' @param lambda Nonnegative penalty parameter.
#' @param P Symmetric positive semidefinite penalty matrix. Required when
#'   `lambda > 0`; the unused vector default from the original script is replaced
#'   by an explicit matrix requirement.
#' @param norm,maxit Retained for source compatibility. Components are always
#'   normalized; this direct covariance algorithm does not iterate.
#' @param tau Target quantile in (0, 1).
#' @param scale.X Whether to standardize predictors before dimension reduction.
#' @param grid Optional quantile grid for the Peng-Huang fit.
#' @return A list containing raw and standardized projection matrices `W` and
#'   `W_scaled`, scores `U`, centering and scaling vectors, the scaled penalty,
#'   residual-coordinate directions, and the number of extracted components.
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
  .validate_tau(tau)
  .scalar_number(lambda, "lambda", lower = 0)
  .scalar_number(ncomp, "ncomp", lower = 1, upper = .Machine$integer.max)
  if (ncomp != as.integer(ncomp)) stop("ncomp must be an integer.")
  if (!is.logical(scale.X) || length(scale.X) != 1L || is.na(scale.X)) {
    stop("scale.X must be TRUE or FALSE.")
  }
  if (!is.numeric(weights) || length(weights) != nrow(X) ||
      any(!is.finite(weights)) || any(weights < 0) || !any(weights > 0)) {
    stop("weights must be finite, nonnegative, and have positive total weight.")
  }
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

  penaltyMatrix <- matrix(0, npred, npred)
  M_root <- NULL
  if (lambda > 0) {
    P <- .validate_penalty(P, npred)
    penaltyMatrix <- P / outer(Xscale, Xscale)
    eig <- eigen(diag(npred) + lambda * penaltyMatrix, symmetric = TRUE)
    if (any(eig$values <= 0)) stop("Penalized metric is not positive definite.")
    M_root <- eig$vectors %*%
      diag(1 / sqrt(eig$values), npred, npred) %*% t(eig$vectors)
  }

  Ymean <- mean(Y)
  qfit <- .fit_crq(matrix(numeric(0), length(Y), 0), Y - Ymean, delta, grid)
  Y0_tau <- .crq_coef(qfit, tau, 1L)
  if (any(!is.finite(Y0_tau))) {
    stop("Target tau is outside the estimable censored quantile range.")
  }
  W_scaled <- W_raw <- directions <- matrix(0, npred, ncomp)
  U <- matrix(0, nrow(X), ncomp)
  Xwork <- X0
  A <- diag(npred)
  used <- 0L
  for (a in seq_len(ncomp)) {
    w <- as.numeric(cq_cov(Xwork, Y - Ymean, Y0_tau, weights, tau))
    if (!is.null(M_root)) w <- as.numeric(M_root %*% w)
    magnitude <- sqrt(sum(w^2))
    if (!is.finite(magnitude) || magnitude < tol) break
    w <- w / magnitude
    score <- Xwork %*% w
    den <- sum(score^2)
    if (!is.finite(den) || den < tol) break
    # Map each residual-coordinate direction back to the original design.
    projection <- as.numeric(A %*% w)
    loading <- crossprod(score, Xwork) / den
    Xwork <- Xwork - score %*% loading
    A <- A - projection %*% loading
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
       directions = directions, ncomp = used,
       tau = tau, Y_quantile = as.numeric(Y0_tau) + Ymean)
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

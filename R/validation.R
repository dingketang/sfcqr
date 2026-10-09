.numeric_matrix <- function(x, name, nrow_expected = NULL, allow_empty = FALSE) {
  if (is.data.frame(x)) {
    if (!all(vapply(x, function(z) is.numeric(z) && !is.complex(z), logical(1)))) {
      stop(name, " must contain only numeric columns.", call. = FALSE)
    }
    x <- as.matrix(x)
  } else if (is.numeric(x) && !is.complex(x) && is.null(dim(x))) {
    x <- matrix(x, ncol = 1L)
  }
  if (!is.matrix(x) || !is.numeric(x) || is.complex(x)) {
    stop(name, " must be a numeric matrix, data frame, or vector.", call. = FALSE)
  }
  if (nrow(x) < 1L || (!allow_empty && ncol(x) < 1L)) {
    stop(name, " must have at least one row",
         if (allow_empty) "." else " and one column.", call. = FALSE)
  }
  if (!is.null(nrow_expected) && nrow(x) != nrow_expected) {
    stop("nrow(", name, ") must equal ", nrow_expected, ".", call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop(name, " must contain only finite values.", call. = FALSE)
  }
  storage.mode(x) <- "double"
  x
}

.numeric_vector <- function(x, name, allow_empty = FALSE) {
  if (is.matrix(x) && ncol(x) == 1L) x <- as.vector(x)
  if (!is.numeric(x) || is.complex(x) || !is.null(dim(x)) ||
      (!allow_empty && length(x) == 0L) || any(!is.finite(x))) {
    stop(name, " must be a ", if (allow_empty) "" else "nonempty ",
         "finite numeric vector.", call. = FALSE)
  }
  as.numeric(x)
}

.scalar_number <- function(x, name, lower = -Inf, upper = Inf) {
  if (!is.numeric(x) || is.complex(x) || length(x) != 1L ||
      !is.finite(x) || x < lower || x > upper) {
    stop(name, " must be a finite numeric scalar between ", lower,
         " and ", upper, ".", call. = FALSE)
  }
  as.numeric(x)
}

.validate_response <- function(Y, delta, positive = TRUE) {
  Y <- .numeric_vector(Y, "Y")
  if (positive && any(Y <= 0)) {
    stop("All Y must be strictly positive.", call. = FALSE)
  }
  if (is.logical(delta)) delta <- as.integer(delta)
  delta <- .numeric_vector(delta, "delta")
  if (length(delta) != length(Y) || any(!delta %in% c(0, 1))) {
    stop("delta must contain 0 or 1 and have the same length as Y.", call. = FALSE)
  }
  list(Y = Y, delta = delta)
}

.validate_tau <- function(tau) {
  tau <- .scalar_number(tau, "tau")
  if (tau <= 0 || tau >= 1) {
    stop("tau must be strictly between 0 and 1.", call. = FALSE)
  }
  tau
}

.logical_scalar <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }
  x
}

.positive_integer <- function(x, name) {
  x <- .scalar_number(x, name, lower = 1, upper = .Machine$integer.max)
  if (x != floor(x)) stop(name, " must be an integer.", call. = FALSE)
  as.integer(x)
}

#' Cumulative quantile hazard transform
#'
#' Evaluate the transform used by the Peng--Huang estimating equation.
#'
#' @param tau Numeric quantile levels in \eqn{[0,1)}. Zero is accepted because
#'   the integration starts at zero.
#' @return A numeric vector, \eqn{-\log(1-\tau)}.
#' @export
#' @examples
#' H_tau(c(0, 0.25, 0.5))
H_tau <- function(tau) {
  tau <- .numeric_vector(tau, "tau", allow_empty = TRUE)
  if (any(tau < 0 | tau >= 1)) {
    stop("tau values must be in [0, 1).", call. = FALSE)
  }
  -log1p(-tau)
}

.H_tau <- H_tau

#' Standard Gaussian kernel
#'
#' @param u A finite numeric vector.
#' @return The standard normal density at each element of \code{u}.
#' @export
#' @examples
#' gaussian_kernel(c(-1, 0, 1))
gaussian_kernel <- function(u) {
  u <- .numeric_vector(u, "u", allow_empty = TRUE)
  exp(-0.5 * u^2) / sqrt(2 * pi)
}

#' Reshape an image vector
#'
#' @param vec A finite numeric vector of exactly \code{nrow * ncol} elements.
#' @param nrow,ncol Positive integer image dimensions.
#' @return A numeric matrix filled by columns, following R's default image
#'   vectorization convention.
#' @export
#' @examples
#' vec_to_image(1:12, nrow = 3, ncol = 4)
vec_to_image <- function(vec, nrow = 100, ncol = 100) {
  vec <- .numeric_vector(vec, "vec")
  nrow <- .positive_integer(nrow, "nrow")
  ncol <- .positive_integer(ncol, "ncol")
  if (length(vec) != as.double(nrow) * ncol) {
    stop("length(vec) must equal nrow * ncol.", call. = FALSE)
  }
  matrix(vec, nrow = nrow, ncol = ncol, byrow = FALSE)
}

#' Quantile covariance
#'
#' @param x A finite numeric vector or matrix of predictors.
#' @param y A finite numeric response vector.
#' @param tau A single quantile level strictly between zero and one.
#' @return The sample covariance between \code{x} and the response indicator
#'   \code{as.numeric(y > quantile(y, tau)) - tau}. The quantile uses R's
#'   default type 7 interpolation.
#' @export
#' @examples
#' q_cov(1:6, c(1, 3, 2, 6, 4, 5), tau = 0.5)
q_cov <- function(x, y, tau) {
  y <- .numeric_vector(y, "y")
  if (length(y) < 2L) stop("At least two observations are required.", call. = FALSE)
  tau <- .validate_tau(tau)
  if (is.null(dim(x))) {
    x <- .numeric_vector(x, "x")
    if (length(x) != length(y)) stop("x and y must have the same length.", call. = FALSE)
  } else {
    x <- .numeric_matrix(x, "x", nrow_expected = length(y))
  }
  indicator <- as.numeric(y > stats::quantile(y, tau, names = FALSE)) - tau
  stats::cov(x, indicator)
}

#' Weighted censored quantile covariance
#'
#' @param x A finite numeric vector or matrix of predictors.
#' @param y A finite numeric response vector.
#' @param y_tau A finite scalar estimated response quantile.
#' @param weights A finite, nonnegative vector of observation weights.
#' @param tau A single quantile level strictly between zero and one.
#' @return The sample covariance between \code{x} and
#'   \code{weights * (as.numeric(y > y_tau) - (1 - tau))}.
#' @export
#' @examples
#' cq_cov(1:6, c(1, 3, 2, 6, 4, 5), y_tau = 3.5,
#'        weights = rep(1, 6), tau = 0.5)
cq_cov <- function(x, y, y_tau, weights, tau) {
  y <- .numeric_vector(y, "y")
  if (length(y) < 2L) stop("At least two observations are required.", call. = FALSE)
  tau <- .validate_tau(tau)
  y_tau <- .scalar_number(y_tau, "y_tau")
  weights <- .numeric_vector(weights, "weights")
  if (length(weights) != length(y) || any(weights < 0)) {
    stop("weights must be nonnegative and have the same length as y.", call. = FALSE)
  }
  if (is.null(dim(x))) {
    x <- .numeric_vector(x, "x")
    if (length(x) != length(y)) stop("x and y must have the same length.", call. = FALSE)
  } else {
    x <- .numeric_matrix(x, "x", nrow_expected = length(y))
  }
  stats::cov(x, (as.numeric(y > y_tau) - (1 - tau)) * weights)
}

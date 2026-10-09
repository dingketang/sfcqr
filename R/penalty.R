.image_gradient <- function(image, axis) {
  nr <- nrow(image)
  nc <- ncol(image)
  row_before <- c(1L, seq_len(nr - 1L))
  row_after <- c(seq.int(2L, nr), nr)
  col_before <- c(1L, seq_len(nc - 1L))
  col_after <- c(seq.int(2L, nc), nc)
  # CImg's float-rounded constants, used by imager::imgradient(scheme = 3).
  a <- 0.14644661545753479
  b <- 0.20710676908493042
  if (axis == "x") {
    difference <- image[row_after, , drop = FALSE] - image[row_before, , drop = FALSE]
    a * (difference[, col_before, drop = FALSE] + difference[, col_after, drop = FALSE]) +
      b * difference
  } else {
    difference <- image[, col_after, drop = FALSE] - image[, col_before, drop = FALSE]
    a * (difference[row_before, , drop = FALSE] + difference[row_after, , drop = FALSE]) +
      b * difference
  }
}

#' Image basis roughness penalty
#'
#' Construct the Gram matrix of numerical second derivatives of image basis
#' functions, preserving the penalty in the supplied research code.
#'
#' @param basis_list A nonempty list of finite numeric matrices. Every matrix
#'   must have the same rectangular dimensions, with at least two rows and
#'   two columns.
#' @return A symmetric positive semidefinite matrix with one row and column
#'   for each basis image.
#' @details Each first derivative uses the rotation invariant central
#'   difference stencil of \code{imager::imgradient(scheme = 3)}, with unit
#'   pixel spacing and replicated values outside the image boundaries.
#'   Applying that operator twice gives the second derivatives. The penalty
#'   entry is the pixel sum of products of the two \eqn{xx} derivatives plus
#'   the corresponding sum for \eqn{yy}; no mixed derivative is included.
#'   This implementation does not require the imager package.
#' @export
#' @examples
#' basis <- list(matrix(1, 3, 4), matrix(seq_len(12), 3, 4))
#' P <- get_P_matrix(basis)
#' stopifnot(isTRUE(all.equal(P, t(P))))
get_P_matrix <- function(basis_list) {
  if (!is.list(basis_list) || length(basis_list) == 0L) {
    stop("basis_list must be a nonempty list of numeric image matrices.", call. = FALSE)
  }
  dimensions <- NULL
  derivatives <- vector("list", length(basis_list))
  for (i in seq_along(basis_list)) {
    image <- .numeric_matrix(basis_list[[i]], paste0("basis_list[[", i, "]]"))
    if (nrow(image) < 2L || ncol(image) < 2L) {
      stop("Each basis image must have at least two rows and two columns.", call. = FALSE)
    }
    if (is.null(dimensions)) dimensions <- dim(image)
    if (!identical(dim(image), dimensions)) {
      stop("All basis images must have the same dimensions.", call. = FALSE)
    }
    dxx <- .image_gradient(.image_gradient(image, "x"), "x")
    dyy <- .image_gradient(.image_gradient(image, "y"), "y")
    derivatives[[i]] <- c(as.vector(dxx), as.vector(dyy))
  }
  derivative_matrix <- do.call(cbind, derivatives)
  if (any(!is.finite(derivative_matrix))) {
    stop("Image derivatives overflow; rescale the basis images.", call. = FALSE)
  }
  result <- crossprod(derivative_matrix)
  if (any(!is.finite(result))) {
    stop("The roughness penalty overflows; rescale the basis images.", call. = FALSE)
  }
  if (!is.null(names(basis_list))) {
    dimnames(result) <- list(names(basis_list), names(basis_list))
  }
  result
}

#' Matrix inverse or Moore--Penrose generalized inverse
#'
#' @param X A nonempty finite numeric or complex matrix. A numeric vector is
#'   interpreted as a one column matrix.
#' @param tol A finite nonnegative relative tolerance for retaining singular
#'   values when a generalized inverse is needed.
#' @return The ordinary inverse when a square matrix can be solved; otherwise
#'   its singular value decomposition generalized inverse, with dimensions
#'   \code{ncol(X)} by \code{nrow(X)}.
#' @details A warning reports the effective rank when an ordinary inverse
#'   fails or a rectangular matrix is rank deficient.
#' @export
#' @examples
#' Minverse(diag(c(2, 4)))
#' X <- matrix(c(1, 0, 0, 1, 1, 1), nrow = 3)
#' stopifnot(isTRUE(all.equal(Minverse(X) %*% X, diag(2))))
Minverse <- function(X, tol = sqrt(.Machine$double.eps)) {
  if (is.data.frame(X)) X <- .numeric_matrix(X, "X")
  if (is.numeric(X) && is.null(dim(X))) X <- matrix(X, ncol = 1L)
  if (!is.matrix(X) || !(is.numeric(X) || is.complex(X)) ||
      any(dim(X) == 0L) || any(!is.finite(X))) {
    stop("X must be a nonempty finite numeric or complex matrix.", call. = FALSE)
  }
  tol <- .scalar_number(tol, "tol", lower = 0)
  inverse_failed <- FALSE
  if (nrow(X) == ncol(X)) {
    inverse <- tryCatch(solve(X), error = function(e) NULL)
    if (!is.null(inverse) && all(is.finite(inverse))) return(inverse)
    inverse_failed <- TRUE
  }
  decomposition <- svd(X)
  retained <- decomposition$d > max(tol * decomposition$d[1L], 0)
  rank <- sum(retained)
  if (inverse_failed || rank < min(dim(X))) {
    warning("X is computationally singular; using an SVD generalized inverse ",
            "with effective rank ", rank, ".", call. = FALSE)
  }
  if (rank == 0L) {
    return(matrix(if (is.complex(X)) 0 + 0i else 0, ncol(X), nrow(X)))
  }
  decomposition$v[, retained, drop = FALSE] %*%
    ((1 / decomposition$d[retained]) * Conj(t(decomposition$u[, retained, drop = FALSE])))
}

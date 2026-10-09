#' Simulate the supplied main manuscript setting
#'
#' Reproduces the \code{scenario = "main"} branch of the supplied
#' \code{gene_data()} research code, including its random draw order and
#' fixed censoring shift. Images have 100 by 100 pixels and use all 36 tensor
#' products of six orthonormalized cubic B-spline functions per dimension.
#' Scalar covariates are independent Uniform(0, 1) draws, and image scores
#' are independent Uniform(0, 4) draws multiplied by the basis index to the
#' power -1/4.
#'
#' The log event time is
#' \deqn{\log T = \langle Z,C_2\rangle + \langle Z,C_1\rangle\epsilon +
#'                 X_1\epsilon + 0.5X_2,}
#' where \eqn{\epsilon} is Uniform(-1, 1),
#' \eqn{C_1 = B_1 + B_8}, and \eqn{C_2 = 0.8B_{22} + 1.2B_{29}}.
#' Inner products are sums over pixels, consistent with the supplied code.
#' Event times are truncated below at \eqn{10^{-10}}. Log censoring times are
#' \eqn{4X_1 + 4X_2 + U + \mathtt{censor\_shift}}, with
#' \eqn{U \sim \operatorname{Logistic}(0, 4)} independently of the event error
#' and image scores. The default shift is -2, exactly as in the supplied
#' main setting. The observed censoring proportion varies randomly.
#'
#' The univariate spline basis is
#' \code{splines::bs(1:100, df = 6, degree = 3, intercept = TRUE)}, followed
#' by \code{qr.Q(qr(...))}. Tensor columns follow nested loops with the first
#' basis index outermost and the second index varying fastest. Column
#' \eqn{6(i-1)+j} is the column-major vectorization of
#' \code{matrix(Q[, i]) \%*\% t(matrix(Q[, j]))}. Each tensor column is
#' normalized to Euclidean norm one, as in the research code.
#'
#' @param n Number of observations, at least 30. The manuscript uses 500,
#'   1000, and 2000; the supplied main code defaults to 1000.
#' @param seed Optional nonnegative integer random seed. When supplied, the
#'   caller's random number state is restored on exit.
#' @param censor_shift Additive shift of the log censoring time. The supplied
#'   main setting uses -2.
#' @return A list with the original \code{Y}, \code{X}, logical \code{delta},
#'   \code{Z}, \code{Bxy}, \code{Ti}, \code{A}, and \code{var_epsilon} fields.
#'   \code{Y}, \code{Ti}, and \code{delta} retain their original single-column
#'   matrix shapes; the model fitter accepts them without conversion.
#'   Additional fields include \code{image_dim}, \code{event_time} (an alias
#'   of \code{Ti}), \code{censoring_time}, latent log times,
#'   \code{latent_scores}, and \code{error}. The \code{truth} list contains
#'   coefficient vectors \code{C1} and \code{C2} in the 36-dimensional basis,
#'   \code{taus}, a 36 by 3 matrix \code{coef_tau} of image basis coefficients,
#'   a 2 by 3 matrix \code{coef_X_tau} of scalar coefficients, and zero
#'   \code{intercept_tau}. Columns correspond to 0.3, 0.5, and 0.7.
#'   \code{setting} records basis and censoring choices.
#' @export
#' @examples
#' d <- simulate_sfcqr_main(n = 40, seed = 2026)
#' dim(d$Z)
#' stopifnot(max(abs(crossprod(d$Bxy) - diag(36))) < 1e-10)
#' d$truth$coef_X_tau
simulate_sfcqr_main <- function(n = 1000, seed = 2026, censor_shift = -2) {
  .scalar_number(n, "n", lower = 30, upper = .Machine$integer.max)
  if (n != as.integer(n)) stop("n must be an integer.", call. = FALSE)
  .scalar_number(censor_shift, "censor_shift")
  if (!is.null(seed)) {
    .scalar_number(seed, "seed", lower = 0, upper = .Machine$integer.max)
    if (seed != as.integer(seed)) stop("seed must be an integer.", call. = FALSE)
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  X1 <- stats::runif(n)
  X2 <- stats::runif(n)
  X <- cbind(X1, X2)
  B <- splines::bs(1:100, df = 6, degree = 3, intercept = TRUE)
  Borth <- qr.Q(qr(B))
  Bxy <- matrix(0, nrow = 10000, ncol = 36)
  for (i in 0:5) {
    for (j in 0:5) {
      Bxy[, i * 6 + j + 1] <- c(matrix(Borth[, i + 1]) %*%
                                 t(matrix(Borth[, j + 1])))
    }
  }
  Bxy <- Bxy %*% diag(1 / sqrt(diag(t(Bxy) %*% Bxy)))
  A <- matrix(stats::runif(n * 36) * 4, ncol = 36, nrow = n)
  coefz <- diag(sqrt(sqrt(1 / (1:36))))
  latent_scores <- A %*% coefz
  Z <- latent_scores %*% t(Bxy)
  image_C1 <- Bxy[, 1] + Bxy[, 8]
  image_C2 <- 0.8 * Bxy[, 22] + 1.2 * Bxy[, 29]
  epsilon <- stats::runif(n) * 2 - 1
  var_epsilon <- stats::var(epsilon)
  logTi <- Z %*% image_C2 + Z %*% image_C1 * epsilon + X1 * epsilon + X2 * 0.5
  Ti <- pmax(exp(logTi), 1e-10)
  u <- stats::rlogis(n, location = 0, scale = 4)
  logC <- 4 * X1 + 4 * X2 + u + censor_shift
  C <- exp(logC)
  delta <- C > Ti
  Y <- Ti
  Y[!delta] <- C[!delta]

  basis_names <- paste0("B", 1:36)
  C1 <- C2 <- stats::setNames(numeric(36), basis_names)
  C1[c(1, 8)] <- 1
  C2[c(22, 29)] <- c(0.8, 1.2)
  taus <- c(0.3, 0.5, 0.7)
  quantile_names <- paste0("tau", taus)
  coef_tau <- vapply(taus, function(tau) (2 * tau - 1) * C1 + C2, numeric(36))
  dimnames(coef_tau) <- list(basis_names, quantile_names)
  coef_X_tau <- rbind(X1 = 2 * taus - 1, X2 = rep(0.5, length(taus)))
  colnames(coef_X_tau) <- quantile_names
  tensor_index <- data.frame(row = rep(1:6, each = 6), col = rep(1:6, 6))

  list(Y = Y, X = X, delta = delta, Z = Z, Bxy = Bxy, Ti = Ti, A = A,
       var_epsilon = var_epsilon, image_dim = c(100L, 100L),
       event_time = Ti, censoring_time = C,
       log_event_time = as.numeric(logTi), log_censoring_time = logC,
       latent_scores = latent_scores, error = epsilon,
       truth = list(C1 = C1, C2 = C2, image_C1 = image_C1, image_C2 = image_C2,
                    taus = taus, coef_tau = coef_tau, coef_X_tau = coef_X_tau,
                    intercept_tau = stats::setNames(rep(0, 3), quantile_names)),
       setting = list(grid = 1:100, knots = attr(B, "knots"),
                      boundary_knots = attr(B, "Boundary.knots"),
                      spline_degree = 3L, basis_per_dimension = 6L,
                      tensor_index = tensor_index, score_decay = -1 / 4,
                      censoring_location = 0, censoring_scale = 4,
                      censor_shift = censor_shift))
}

#' Compatibility generator for the supplied main simulation code
#'
#' Calls [simulate_sfcqr_main()] with the main manuscript setting and a fixed
#' log censoring shift of -2. Only the main scenario is exposed by this wrapper.
#'
#' @param seed Nonnegative integer random seed.
#' @param n Optional sample size. \code{NULL} uses the original main default
#'   of 1000 observations.
#' @return The list returned by [simulate_sfcqr_main()], including all original
#'   main-generator fields.
#' @export
#' @examples
#' d <- gene_data(seed = 2026, n = 40)
#' dim(d$Z)
gene_data <- function(seed, n = NULL) {
  if (is.null(n)) n <- 1000
  simulate_sfcqr_main(n = n, seed = seed, censor_shift = -2)
}

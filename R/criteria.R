#' Selection criteria for censored functional quantile regression
#'
#' Evaluate the AIC, BIC, and generalized approximate cross validation
#' criteria implemented in the supplied research code on a fitted quantile
#' coefficient path.
#'
#' @param Y Positive, finite observed times.
#' @param delta Event indicators: one for an observed response event and zero
#'   for a censored response.
#' @param Xfit Finite numeric design matrix excluding the intercept. A matrix
#'   with zero columns represents an intercept only model.
#' @param tau_vec Strictly increasing quantile levels in \eqn{(0,1)}.
#' @param coef_mat Numeric coefficient matrix with \code{ncol(Xfit) + 1} rows
#'   and \code{length(tau_vec)} columns. Its first row is the intercept.
#' @param omega Finite exponent in \eqn{C_n = \log(n)^\omega}.
#' @param H_tau Function mapping quantile levels to cumulative hazard values.
#'   The default is \code{H_tau(tau) = -log1p(-tau)}. It must accept zero,
#'   return finite values, and be nondecreasing along the evaluation grid.
#' @param df_penalty Finite nonnegative effective model dimension. Defaults
#'   to \code{ncol(Xfit)}.
#' @param coef_bound Positive threshold for the sum of absolute coefficients
#'   in each quantile column when \code{remove_bad_tau = TRUE}.
#' @param remove_bad_tau Whether to remove columns with nonfinite coefficients
#'   or an absolute coefficient sum greater than \code{coef_bound}.
#' @param min_tau Positive integer minimum number of retained quantile levels.
#'   The default, 11, preserves the supplied code's \code{L <= 10} rule.
#' @return A list including pointwise criteria, their arithmetic means
#'   \code{IAIC}, \code{IBIC}, and \code{IGACV}, retained quantile levels,
#'   intermediate matrices, and filtering diagnostics. With fewer than
#'   \code{min_tau} retained levels, the three mean criteria are \code{Inf}
#'   and \code{message} describes the reason. A nonpositive or nonfinite
#'   denominator produces \code{Inf} for the affected criterion and a warning.
#' @details The three reported mean criteria average over retained quantile
#'   levels with equal weights. They retain the research code's meaning of
#'   \code{IAIC}, \code{IBIC}, and \code{IGACV}, rather than applying numerical
#'   quadrature over an unequally spaced grid.
#' @export
#' @examples
#' Y <- seq(1, 4, length.out = 20)
#' tau <- seq(0.1, 0.8, length.out = 15)
#' coefficients <- matrix(log(stats::quantile(Y, tau)), nrow = 1)
#' result <- gacv_fcrq(Y, rep(1, 20), matrix(numeric(0), 20, 0),
#'                    tau, coefficients)
#' stopifnot(is.finite(result$IBIC))
gacv_fcrq <- function(Y, delta, Xfit, tau_vec, coef_mat, omega = 1 / 2,
                       H_tau = .H_tau, df_penalty = NULL, coef_bound = 50,
                       remove_bad_tau = TRUE, min_tau = 11L) {
  response <- .validate_response(Y, delta)
  Y <- response$Y
  delta <- response$delta
  n <- length(Y)
  Xfit <- .numeric_matrix(Xfit, "Xfit", nrow_expected = n, allow_empty = TRUE)
  tau_vec <- .numeric_vector(tau_vec, "tau_vec", allow_empty = TRUE)
  if (any(tau_vec <= 0 | tau_vec >= 1) || any(diff(tau_vec) <= 0)) {
    stop("tau_vec must be strictly increasing and strictly between 0 and 1.", call. = FALSE)
  }
  if (!is.matrix(coef_mat) || !is.numeric(coef_mat) || is.complex(coef_mat)) {
    stop("coef_mat must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(coef_mat) != ncol(Xfit) + 1L || ncol(coef_mat) != length(tau_vec)) {
    stop("coef_mat must have ncol(Xfit) + 1 rows and length(tau_vec) columns.", call. = FALSE)
  }
  omega <- .scalar_number(omega, "omega")
  coef_bound <- .scalar_number(coef_bound, "coef_bound", lower = 0)
  if (coef_bound == 0) stop("coef_bound must be strictly positive.", call. = FALSE)
  remove_bad_tau <- .logical_scalar(remove_bad_tau, "remove_bad_tau")
  min_tau <- .positive_integer(min_tau, "min_tau")
  if (!is.function(H_tau)) stop("H_tau must be a function.", call. = FALSE)
  if (is.null(df_penalty)) df_penalty <- ncol(Xfit)
  df_penalty <- .scalar_number(df_penalty, "df_penalty", lower = 0)

  removed <- integer(0)
  if (remove_bad_tau) {
    removed <- which(colSums(!is.finite(coef_mat)) > 0 |
                       colSums(abs(coef_mat)) > coef_bound)
    if (length(removed)) {
      coef_mat <- coef_mat[, -removed, drop = FALSE]
      tau_vec <- tau_vec[-removed]
    }
  } else if (any(!is.finite(coef_mat))) {
    stop("coef_mat must be finite when remove_bad_tau = FALSE.", call. = FALSE)
  }
  L <- length(tau_vec)
  diagnostics <- list(tau_vec = tau_vec, removed_tau_indices = removed,
                      n_tau = L, min_tau = min_tau, df_penalty = df_penalty,
                      omega = omega)
  if (L < min_tau) {
    return(c(list(IAIC = Inf, IBIC = Inf, IGACV = Inf,
                  message = paste0("Only ", L, " quantile levels remain; at least ",
                                   min_tau, " are required.")), diagnostics))
  }

  X_full <- cbind(1, Xfit)
  eta_mat <- X_full %*% coef_mat
  if (any(!is.finite(eta_mat))) {
    stop("Fitted log quantiles overflow; rescale Xfit or coefficient values.", call. = FALSE)
  }
  qhat_mat <- exp(eta_mat)
  Y_mat <- matrix(Y, nrow = n, ncol = L)
  delta_mat <- matrix(delta, nrow = n, ncol = L)
  N_eval_mat <- (delta_mat == 1) * (Y_mat <= qhat_mat)

  previous_tau <- c(0, tau_vec[-L])
  H_vals <- .numeric_vector(H_tau(tau_vec), "H_tau(tau_vec)")
  H_prev <- .numeric_vector(H_tau(previous_tau), "H_tau(previous_tau)")
  if (length(H_vals) != L || length(H_prev) != L) {
    stop("H_tau must return one value per supplied quantile level.", call. = FALSE)
  }
  dH <- H_vals - H_prev
  if (any(!is.finite(dH)) || any(dH < 0)) {
    stop("H_tau increments must be finite and nonnegative.", call. = FALSE)
  }

  qhat_with_tau0 <- cbind(0, qhat_mat)
  I_ge_mat <- qhat_with_tau0[, seq_len(L), drop = FALSE] <= Y_mat
  cum_integral_mat <- matrix(0, nrow = n, ncol = L)
  cumulative <- numeric(n)
  for (j in seq_len(L)) {
    cumulative <- cumulative + I_ge_mat[, j] * dH[j]
    cum_integral_mat[, j] <- cumulative
  }
  logY_mat <- matrix(log(Y), nrow = n, ncol = L)
  pointwise_numerators <- colSums(
    -(logY_mat - eta_mat) * (N_eval_mat - cum_integral_mat)
  )
  if (any(!is.finite(pointwise_numerators))) {
    stop("Criterion numerators overflow; rescale the input values.", call. = FALSE)
  }
  denominators <- c(AIC = n - df_penalty,
                    BIC = n - log(n) * df_penalty,
                    GACV = n - log(n)^omega * df_penalty)
  invalid <- !is.finite(denominators) | denominators <= 0
  if (any(invalid)) {
    warning("Nonpositive or nonfinite denominator for ",
            paste(names(denominators)[invalid], collapse = ", "),
            "; returning Inf for the affected criteria.", call. = FALSE)
  }
  criterion <- function(i) {
    if (invalid[i]) rep(Inf, L) else pointwise_numerators / denominators[i]
  }
  pointwise_aic <- criterion(1L)
  pointwise_bic <- criterion(2L)
  pointwise_gacv <- criterion(3L)
  c(list(eta_mat = eta_mat, qhat_mat = qhat_mat, N_eval_mat = N_eval_mat,
         dH = dH, cum_integral_mat = cum_integral_mat,
         pointwise_numerators = pointwise_numerators,
         pointwise_aic = pointwise_aic, pointwise_bic = pointwise_bic,
         pointwise_gacv = pointwise_gacv,
         IAIC = mean(pointwise_aic), IBIC = mean(pointwise_bic),
         IGACV = mean(pointwise_gacv), denominators = denominators,
         message = NULL), diagnostics)
}

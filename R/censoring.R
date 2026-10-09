.beran_inputs <- function(Y, delta, X, x0, h) {
  response <- .validate_response(Y, delta)
  X <- .numeric_matrix(X, "X", nrow_expected = length(response$Y))
  x0 <- .numeric_vector(x0, "x0")
  if (length(x0) != ncol(X)) {
    stop("x0 must contain one value per covariate column in X.", call. = FALSE)
  }
  h <- .scalar_number(h, "h", lower = 0)
  if (h == 0) stop("Bandwidth h must be strictly positive.", call. = FALSE)
  list(Y = response$Y, delta = response$delta, X = X, x0 = x0, h = h)
}

.beran_curve <- function(Y, delta, X, x0, times, h) {
  # The common normalizing factor cancels from every risk-set ratio. Center
  # squared distances at their minimum before exponentiating, so at least one
  # kernel weight is one even when the unscaled Gaussian weights underflow.
  if (is.null(dim(X))) X <- matrix(X, ncol = 1L)
  difference <- sweep(X, 2, x0, "-")
  distance <- sqrt(rowSums(difference^2))
  if (any(!is.finite(distance))) {
    stop("Covariate distances overflow; rescale X and x0.", call. = FALSE)
  }
  distance_min <- min(distance)
  exponent <- ((distance - distance_min) / h) *
    (distance / h + distance_min / h)
  exponent[distance == distance_min] <- 0
  w <- exp(-0.5 * exponent)
  ord <- order(Y)
  Ys <- Y[ord]
  ws <- w[ord]
  cens <- 1 - delta[ord]
  unique_times <- sort(unique(Ys))
  groups <- match(Ys, unique_times)
  mass <- as.numeric(rowsum(ws, groups, reorder = FALSE))
  events <- as.numeric(rowsum(ws * cens, groups, reorder = FALSE))
  risk <- rev(cumsum(rev(mass)))
  factor <- rep(1, length(risk))
  eligible <- risk > 0
  factor[eligible] <- pmax(0, pmin(1, 1 - events[eligible] / risk[eligible]))
  survival <- cumprod(factor)
  index <- findInterval(times, unique_times)
  c(1, survival)[index + 1L]
}

#' Beran conditional censoring survival estimator
#'
#' Estimate the survival function of censoring time conditional on a single
#' covariate vector using an isotropic Gaussian product kernel and a weighted
#' product limit.
#'
#' @param Y Positive, finite observed times.
#' @param delta Event indicators: one for an observed response event and zero
#'   for a censored response.
#' @param X A numeric covariate vector or an observations-by-covariates matrix.
#' @param x0 Covariate values at which to estimate survival; one per column of X.
#' @param t A finite time at which to evaluate survival.
#' @param h A strictly positive Gaussian kernel bandwidth.
#' @return A scalar censoring survival estimate in \eqn{[0,1]}.
#' @details Survival is evaluated at \eqn{t}, including censoring events at
#'   exactly \eqn{t}; this is the right continuous convention in the supplied
#'   research code. Kernel weights are rescaled to prevent complete numerical
#'   underflow. Response event indicators are used as provided.
#' @export
#' @examples
#' beran_Ghat(c(1, 2, 3, 4), c(1, 0, 1, 1),
#'            X = c(0, 0.1, 0.2, 0.3), x0 = 0.15, t = 2, h = 0.3)
beran_Ghat <- function(Y, delta, X, x0, t, h) {
  input <- .beran_inputs(Y, delta, X, x0, h)
  t <- .scalar_number(t, "t")
  .beran_curve(input$Y, input$delta, input$X, input$x0, t, input$h)
}

#' Evaluate a Beran censoring survival curve
#'
#' @inheritParams beran_Ghat
#' @param times A finite numeric vector of evaluation times, in any order.
#' @return A numeric vector of censoring survival estimates in the same order
#'   as \code{times}.
#' @seealso \code{\link{beran_Ghat}}
#' @export
#' @examples
#' beran_Ghat_curve(c(1, 2, 3, 4), c(1, 0, 1, 1),
#'                  X = c(0, 0.1, 0.2, 0.3), x0 = 0.15,
#'                  times = c(3, 1, 2), h = 0.3)
beran_Ghat_curve <- function(Y, delta, X, x0, times, h) {
  input <- .beran_inputs(Y, delta, X, x0, h)
  times <- .numeric_vector(times, "times", allow_empty = TRUE)
  .beran_curve(input$Y, input$delta, input$X, input$x0, times, input$h)
}

.censor_model <- function(expr, method) {
  tryCatch(expr, error = function(e) {
    stop("The ", method, " censoring model could not be fitted: ",
         conditionMessage(e), call. = FALSE)
  })
}

#' Inverse probability of censoring weights
#'
#' Compute \eqn{\delta_i / \widehat G(Y_i \mid X_i)} from an estimated
#' censoring survival model.
#'
#' @param Y Positive, finite observed times.
#' @param X A finite numeric covariate matrix, data frame, or vector. A matrix
#'   with zero columns requests a marginal censoring model.
#' @param delta Event indicators: one for an observed response event and zero
#'   for a censored response.
#' @param method One of \code{"marginal"}, \code{"lognormal"},
#'   \code{"loglogistic"}, \code{"beran"}, or \code{"cox"}. A \code{NULL}
#'   value selects \code{"loglogistic"}. The legacy spelling \code{"bernan"}
#'   is accepted as an alias for \code{"beran"}.
#' @param bandwidth Positive bandwidth for \code{"beran"}; defaults to
#'   \eqn{n^{-1/(p+4)}}, where \eqn{p} is the number of columns in \code{X}.
#' @param min_survival Finite lower bound on estimated survival, strictly
#'   positive and at most one.
#' @param force_last_event Whether to locally set the indicator for the first
#'   observation attaining the maximum observed time to one. Defaults to
#'   \code{TRUE} to preserve the supplied research code.
#' @param beran_covariates Use \code{"first"} for the original fitting script's
#'   first-column conditioning, or \code{"all"} for the multivariate Gaussian
#'   kernel used in Section S4 of the manuscript.
#' @return A finite numeric vector of nonnegative weights in observation order.
#' @details All methods evaluate the right continuous survival function at
#'   \eqn{Y_i}, including a censoring event at the evaluation time. If
#'   \code{force_last_event = TRUE}, its modified event indicator is used both
#'   to fit censoring survival and in the weight numerator; the caller's input
#'   is not modified. With no censored observations all weights equal one.
#'
#'   By default the Beran method uses the first covariate only. Its bandwidth
#'   retains the original code's heuristic based on the full number of
#'   covariates; specify \code{bandwidth} to control this directly. Parametric
#'   and Cox models use every covariate column. These modeling choices and the
#'   last event override should be assessed for the application.
#' @export
#' @examples
#' Y <- c(1, 2, 3, 4)
#' delta <- c(1, 0, 1, 1)
#' get_weights(Y, matrix(numeric(0), 4, 0), delta, method = "marginal")
#' get_weights(Y, c(0, 0.1, 0.2, 0.3), delta, method = "beran")
get_weights <- function(Y, X, delta, method = NULL, bandwidth = NULL,
                        min_survival = 1e-4, force_last_event = TRUE,
                        beran_covariates = c("first", "all")) {
  response <- .validate_response(Y, delta)
  Y <- response$Y
  delta <- response$delta
  n <- length(Y)
  X <- .numeric_matrix(X, "X", nrow_expected = n, allow_empty = TRUE)
  if (is.null(method)) method <- "loglogistic"
  choices <- c("marginal", "lognormal", "loglogistic", "beran", "cox")
  if (!is.character(method) || length(method) != 1L || is.na(method)) {
    stop("method must be one of: ", paste(choices, collapse = ", "), ".", call. = FALSE)
  }
  if (method == "bernan") method <- "beran"
  if (!method %in% choices) {
    stop("method must be one of: ", paste(choices, collapse = ", "), ".", call. = FALSE)
  }
  min_survival <- .scalar_number(min_survival, "min_survival", lower = 0, upper = 1)
  if (min_survival == 0) stop("min_survival must be strictly positive.", call. = FALSE)
  force_last_event <- .logical_scalar(force_last_event, "force_last_event")
  beran_covariates <- match.arg(beran_covariates)
  if (!is.null(bandwidth)) {
    bandwidth <- .scalar_number(bandwidth, "bandwidth", lower = 0)
    if (bandwidth == 0) stop("bandwidth must be strictly positive.", call. = FALSE)
  }
  if (all(delta == 1)) return(rep(1, n))
  if (force_last_event) delta[which.max(Y)] <- 1
  if (all(delta == 1)) return(rep(1, n))
  if (ncol(X) == 0L) method <- "marginal"

  if (method == "marginal") {
    fit <- .censor_model(survival::survfit(survival::Surv(Y, 1 - delta) ~ 1), method)
    index <- findInterval(Y, fit$time)
    Ghat <- c(1, as.numeric(fit$surv))[index + 1L]
  } else if (method == "beran") {
    if (is.null(bandwidth)) bandwidth <- n^(-1 / (ncol(X) + 4))
    X_kernel <- if (beran_covariates == "all") X else X[, 1L, drop = FALSE]
    Ghat <- vapply(seq_len(n), function(i) {
      .beran_curve(Y, delta, X_kernel, X_kernel[i, ], Y[i], bandwidth)
    }, numeric(1))
  } else {
    dat <- data.frame(Y = Y, cens = 1 - delta, X)
    names(dat) <- c("Y", "cens", paste0("X", seq_len(ncol(X))))
    formula <- survival::Surv(Y, cens) ~ .
    if (method %in% c("lognormal", "loglogistic")) {
      fit <- .censor_model(survival::survreg(formula, data = dat, dist = method), method)
      coefficients <- stats::coef(fit)
      if (length(coefficients) != ncol(X) + 1L || any(!is.finite(coefficients)) ||
          length(fit$scale) != 1L || !is.finite(fit$scale) || fit$scale <= 0) {
        stop("The ", method, " censoring model has nonfinite coefficients or an invalid scale; ",
             "check collinearity and the number of censoring events.", call. = FALSE)
      }
      lp <- as.numeric(cbind(1, X) %*% coefficients)
      z <- (log(Y) - lp) / fit$scale
      Ghat <- if (method == "lognormal") {
        stats::pnorm(z, lower.tail = FALSE)
      } else {
        stats::plogis(z, lower.tail = FALSE)
      }
    } else {
      fit <- .censor_model(survival::coxph(formula, data = dat), method)
      coefficients <- stats::coef(fit)
      if (length(coefficients) != ncol(X) || any(!is.finite(coefficients))) {
        stop("The Cox censoring model has nonfinite coefficients; check ",
             "collinearity and the number of censoring events.", call. = FALSE)
      }
      # Centered baseline hazard and linear predictors use the same reference,
      # avoiding unnecessarily large exp(lp) from an uncentered baseline.
      baseline <- .censor_model(survival::basehaz(fit, centered = TRUE), method)
      index <- findInterval(Y, baseline$time)
      H0 <- c(0, baseline$hazard)[index + 1L]
      lp <- as.numeric(stats::predict(fit, type = "lp"))
      cumulative_hazard <- numeric(n)
      positive <- H0 > 0
      cumulative_hazard[positive] <- exp(log(H0[positive]) + lp[positive])
      Ghat <- exp(-cumulative_hazard)
    }
  }
  if (length(Ghat) != n || any(!is.finite(Ghat)) || any(Ghat < 0 | Ghat > 1)) {
    stop("The ", method, " censoring model produced invalid survival estimates.", call. = FALSE)
  }
  as.numeric(delta / pmax(Ghat, min_survival))
}

.fit_crq <- function(design, logY, delta, grid = NULL) {
  dat <- data.frame(.time = logY, .status = delta)
  if (ncol(design)) {
    predictors <- as.data.frame(design)
    names(predictors) <- paste0("V", seq_len(ncol(design)))
    dat <- cbind(dat, predictors)
    form <- stats::reformulate(names(predictors),
      response = "survival::Surv(.time, .status, type = 'right')")
  } else {
    form <- stats::as.formula("survival::Surv(.time, .status, type = 'right') ~ 1")
  }
  if (is.null(grid)) quantreg::crq(form, data = dat, method = "PengHuang")
  else quantreg::crq(form, data = dat, method = "PengHuang", grid = grid)
}

.crq_coef <- function(fit, tau, size) {
  b <- as.numeric(stats::coef(fit, taus = tau))
  if (length(b) != size) stop("Unexpected coefficient dimension from quantreg::crq.")
  b
}

.capture_fit <- function(expr) {
  warnings <- character(0)
  result <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) e)
  list(value = result, warnings = unique(warnings),
       error = if (inherits(result, "error")) conditionMessage(result) else NULL)
}

.image_dimensions <- function(image_dim, pixels, required = TRUE) {
  if (is.null(image_dim)) {
    side <- sqrt(pixels)
    if (side == as.integer(side)) image_dim <- c(side, side)
    else if (required) stop("Supply image_dim for a nonsquare image grid.")
    else return(NULL)
  }
  if (!is.numeric(image_dim) || length(image_dim) != 2L ||
      any(!is.finite(image_dim)) || any(image_dim < 2) ||
      any(image_dim != as.integer(image_dim)) || prod(image_dim) != pixels) {
    stop("image_dim must contain two integers >= 2 whose product matches ncol(Z).")
  }
  as.integer(image_dim)
}

#' Fit supervised functional censored quantile regression
#'
#' Fits the image-basis model supplied in the original research script. Selects
#' component count and smoothing parameter using an integrated criterion.
#' All candidates share a quantile grid; invalid candidates receive infinite
#' scores and their reasons are stored in `diagnostics`.
#'
#' @param data List with positive observed times `Y`, event indicators `delta`,
#'   scalar covariates `X` (n by p), vectorized images `Z` (n by pixels), and image
#'   basis `Bxy` (pixels by basis functions). `X` may be NULL or have zero columns.
#'   Pixels follow R's column-major matrix order.
#' @param tau Target quantile in (0, 1).
#' @param Mmax Maximum number of supervised components, capped at the available
#'   basis dimension and sample size.
#' @param lambda_list Nonnegative smoothing candidates.
#' @param image_dim Image rows and columns; inferred for square grids or read
#'   from `data$image_dim`.
#' @param censoring_method Censoring model passed to [get_weights()]. With no
#'   scalar covariates, marginal weights are used.
#' @param bandwidth,min_survival,force_last_event Passed to [get_weights()].
#' @param scale.X Standardize basis predictors during component extraction.
#' @param penalty Optional supplied positive semidefinite basis penalty matrix.
#' @param grid Optional strictly increasing quantile grid in (0, 1) spanning
#'   `tau`. By default, uses the grid rule of `quantreg` with `tau` included.
#' @param criterion Selection criterion: `"IBIC"` (the original fitting script),
#'   `"IGACV"` (the manuscript), or `"IAIC"`.
#' @param beran_covariates `"first"` retains the original fitting script;
#'   `"all"` conditions on every scalar covariate as in the manuscript.
#' @return An object of class `sfcqr` containing `coefficients`, the raw-design
#'   `intercept`, scalar coefficients `coef_X`, basis coefficients
#'   `coef_estimate`, selected `M` and `lambda`, the criterion matrix `tuning`,
#'   `diagnostics`, training scores, weights, basis, and fitted quantile model.
#'   Use [predict.sfcqr()] and [coef.sfcqr()] for predictions and image effects.
#' @rdname sfcqr-fit
#' @export
#' @examples
#' d <- simulate_sfcqr_data(n = 120, seed = 2026)
#' fit <- sfcqr(d, tau = 0.5, Mmax = 2, lambda_list = c(0, 1))
#' print(fit)
#' head(predict(fit))
sfcqr <- function(data, tau = 0.5, Mmax = 30,
                  lambda_list = c(0, 1, 10, 100, 1000, 1e6),
                  image_dim = NULL, censoring_method = "beran", bandwidth = NULL,
                  min_survival = 1e-4, force_last_event = TRUE, scale.X = TRUE,
                  penalty = NULL, grid = NULL,
                  criterion = c("IBIC", "IGACV", "IAIC"),
                  beran_covariates = c("first", "all")) {
  if (!is.list(data) || !all(c("Y", "delta", "Z", "Bxy") %in% names(data))) {
    stop("data must be a list containing Y, delta, Z, and Bxy.")
  }
  response <- .validate_response(data$Y, data$delta)
  Y <- response$Y
  delta <- response$delta
  n <- length(Y)
  .validate_tau(tau)
  criterion <- match.arg(criterion)
  beran_covariates <- match.arg(beran_covariates)
  .scalar_number(Mmax, "Mmax", lower = 1, upper = .Machine$integer.max)
  if (Mmax != as.integer(Mmax)) stop("Mmax must be an integer.")
  if (!is.numeric(lambda_list) || !length(lambda_list) ||
      any(!is.finite(lambda_list)) || any(lambda_list < 0)) {
    stop("lambda_list must contain finite nonnegative values.")
  }
  lambda_list <- unique(lambda_list)
  X <- if (is.null(data$X)) matrix(numeric(0), n, 0L) else
    .numeric_matrix(data$X, "data$X", nrow_expected = n, allow_empty = TRUE)
  Z <- .numeric_matrix(data$Z, "data$Z", nrow_expected = n)
  Bxy <- .numeric_matrix(data$Bxy, "data$Bxy", nrow_expected = ncol(Z))
  p <- ncol(X)
  q <- ncol(Bxy)
  if (qr(cbind(1, X))$rank < p + 1L) {
    stop("Scalar covariates must be linearly independent of each other and the intercept.")
  }
  Mmax <- min(as.integer(Mmax), q, n - p - 2L)
  if (Mmax < 1L) stop("Too few observations for scalar covariates and components.")
  if (is.null(image_dim)) image_dim <- data$image_dim
  image_dim <- .image_dimensions(image_dim, ncol(Z),
                                  required = is.null(penalty) && any(lambda_list > 0))
  if (is.null(penalty)) {
    penalty <- if (any(lambda_list > 0)) get_P_matrix(lapply(seq_len(q), function(j) {
      vec_to_image(Bxy[, j], image_dim[1], image_dim[2])
    })) else matrix(0, q, q)
  } else penalty <- .validate_penalty(penalty, q)
  if (is.null(grid)) {
    grid <- sort(unique(c(seq(1/n, 1 - 1/n,
                              by = min(0.01, 1/(2 * n^0.7))), tau)))
  }
  if (!is.numeric(grid) || length(grid) < 3L || any(!is.finite(grid)) ||
      any(grid <= 0 | grid >= 1) || any(diff(grid) <= 0) ||
      min(grid) >= tau || max(grid) <= tau) {
    stop("grid must be strictly increasing in (0, 1) and span tau.")
  }
  weight_result <- .capture_fit(get_weights(Y, X, delta,
    method = censoring_method, bandwidth = bandwidth, min_survival = min_survival,
    force_last_event = force_last_event, beran_covariates = beran_covariates))
  if (!is.null(weight_result$error)) stop(weight_result$error)
  weights <- weight_result$value
  ZB <- Z %*% Bxy
  tuning <- matrix(Inf, Mmax, length(lambda_list),
    dimnames = list(paste0("M", seq_len(Mmax)), paste0("lambda=", lambda_list)))
  diagnostics <- matrix("Not evaluated", Mmax, length(lambda_list),
                        dimnames = dimnames(tuning))
  best <- NULL
  best_ibic <- Inf
  for (j in seq_along(lambda_list)) {
    component_result <- .capture_fit(plsfit_cqcov(ZB, log(Y), Mmax,
      weights, delta, lambda = lambda_list[j], P = penalty, tau = tau,
      scale.X = scale.X, grid = grid))
    if (!is.null(component_result$error)) {
      diagnostics[, j] <- component_result$error
      next
    }
    components <- component_result$value
    if (components$ncomp < Mmax) {
      diagnostics[seq.int(components$ncomp + 1L, Mmax), j] <-
        "No further nonzero components"
    }
    for (i in seq_len(components$ncomp)) {
      design <- cbind(components$U[, seq_len(i), drop = FALSE], X)
      candidate <- .capture_fit({
        if (qr(cbind(1, design))$rank != ncol(design) + 1L) {
          stop("Candidate design is rank deficient.")
        }
        crqfit <- .fit_crq(design, log(Y), delta, grid)
        b <- .crq_coef(crqfit, tau, i + p + 1L)
        if (any(!is.finite(b))) stop("Target tau is outside this candidate's estimable range.")
        sol <- crqfit$sol
        keep <- which(sol[1, ] <= tau & sol[1, ] > 0)
        criteria <- gacv_fcrq(Y, delta, design, sol[1, keep],
          sol[seq.int(2L, ncol(design) + 2L), keep, drop = FALSE],
          H_tau = H_tau, remove_bad_tau = FALSE)
        list(crqfit = crqfit, coefficients = b, criteria = criteria)
      })
      if (!is.null(candidate$error)) {
        diagnostics[i, j] <- candidate$error
        next
      }
      ibic <- candidate$value$criteria[[criterion]]
      tuning[i, j] <- if (is.finite(ibic)) ibic else Inf
      notes <- unique(c(component_result$warnings, candidate$warnings))
      diagnostics[i, j] <- if (length(notes)) paste(notes, collapse = "; ")
        else if (is.finite(ibic)) "OK"
        else if (!is.null(candidate$value$criteria$message)) candidate$value$criteria$message
        else paste(criterion, "is not finite")
      if (is.finite(ibic) && ibic < best_ibic) {
        best_ibic <- ibic
        best <- list(M = i, lambda = lambda_list[j], components = components,
                     candidate = candidate$value)
      }
    }
  }
  if (is.null(best)) {
    reasons <- unique(as.vector(diagnostics))
    stop("No valid candidate model. Try another tau, fewer components, or more ",
         "observations. Diagnostics: ", paste(reasons, collapse = "; "))
  }
  M <- best$M
  b <- best$candidate$coefficients
  projection <- best$components$W[, seq_len(M), drop = FALSE]
  basis_coef <- as.numeric(projection %*% b[seq.int(2L, M + 1L)])
  scalar_coef <- b[M + 1L + seq_len(p)]
  intercept <- b[1] - sum(best$components$Xmean * basis_coef)
  coefficients <- c(intercept, scalar_coef, basis_coef)
  scalar_names <- if (p) paste0("X", seq_len(p)) else character(0)
  names(coefficients) <- c("(Intercept)", scalar_names,
                           paste0("B", seq_len(q)))
  names(scalar_coef) <- scalar_names
  names(basis_coef) <- paste0("B", seq_len(q))
  linear_predictor <- as.numeric(intercept + X %*% scalar_coef + ZB %*% basis_coef)
  structure(list(call = match.call(), tau = tau, M = M, lambda = best$lambda,
    tuning = tuning, diagnostics = diagnostics, weights = weights,
    weight_warnings = weight_result$warnings, crqfit = best$candidate$crqfit,
    criteria = best$candidate$criteria, criterion = criterion,
    beran_covariates = beran_covariates, coefficients = coefficients,
    intercept = unname(intercept), coef_X = scalar_coef, coef_estimate = basis_coef,
    Bxy = Bxy, image_dim = image_dim, projection = projection,
    Xmean = best$components$Xmean, Xscale = best$components$Xscale,
    scores = best$components$U[, seq_len(M), drop = FALSE],
    penalty = penalty, n_scalar = p, n = n,
    linear.predictors = linear_predictor, fitted.values = exp(linear_predictor),
    censoring_method = if (p == 0L) "marginal" else censoring_method,
    force_last_event = force_last_event), class = "sfcqr")
}

#' Compatibility wrapper for the original Sfcqr function
#' @param data,tau,Mmax See [sfcqr()].
#' @param fit_object Return an `sfcqr` fit when TRUE; otherwise return the legacy
#'   vector of component count, raw intercept, scalar and basis coefficients.
#' @param ... Additional arguments passed to [sfcqr()].
#' @return A fitted object or named numeric vector.
#' @rdname sfcqr-legacy
#' @export
Sfcqr <- function(data, tau, Mmax = 30, fit_object = FALSE, ...) {
  fit <- sfcqr(data, tau = tau, Mmax = Mmax, ...)
  if (isTRUE(fit_object)) return(fit)
  result <- c(M_num = fit$M, fit$coefficients)
  names(result)[2L] <- "X0"
  result
}

#' Predict conditional event-time quantiles
#' @param object A fitted [sfcqr()] object.
#' @param newdata List with `Z` and, if the model has scalar predictors, `X`.
#'   With NULL, returns training predictions.
#' @param type `"response"` for time quantiles or `"link"` for log-time quantiles.
#' @param ... Reserved for future extensions.
#' @return Numeric prediction vector at the model's fitted quantile level.
#' @export
predict.sfcqr <- function(object, newdata = NULL, type = c("response", "link"), ...) {
  type <- match.arg(type)
  if (is.null(newdata)) eta <- object$linear.predictors else {
    if (!is.list(newdata) || is.null(newdata$Z)) stop("newdata must be a list with Z.")
    Z <- .numeric_matrix(newdata$Z, "newdata$Z")
    if (ncol(Z) != nrow(object$Bxy)) stop("newdata$Z has the wrong pixel count.")
    X <- if (is.null(newdata$X) && object$n_scalar == 0L) {
      matrix(numeric(0), nrow(Z), 0)
    } else .numeric_matrix(newdata$X, "newdata$X", nrow_expected = nrow(Z),
                          allow_empty = TRUE)
    if (ncol(X) != object$n_scalar) stop("newdata$X has the wrong number of columns.")
    eta <- as.numeric(object$intercept + X %*% object$coef_X +
                        (Z %*% object$Bxy) %*% object$coef_estimate)
  }
  if (type == "link") eta else exp(eta)
}

#' Extract model coefficients
#' @param object A fitted [sfcqr()] object.
#' @param type `"all"` for the raw-design intercept, scalar and basis
#'   coefficients; `"basis"` for basis coefficients; `"image"` for the pixel
#'   coefficient surface.
#' @param ... Reserved for future extensions.
#' @return A named numeric vector, or an image matrix for `type = "image"`.
#' @export
coef.sfcqr <- function(object, type = c("all", "basis", "image"), ...) {
  type <- match.arg(type)
  switch(type, all = object$coefficients, basis = object$coef_estimate,
    image = {
      if (is.null(object$image_dim)) stop("Image dimensions were not supplied.")
      vec_to_image(as.numeric(object$Bxy %*% object$coef_estimate),
                   object$image_dim[1], object$image_dim[2])
    })
}

#' Print a supervised functional censored quantile fit
#' @param x A fitted [sfcqr()] object.
#' @param ... Reserved for future extensions.
#' @return The fit, invisibly.
#' @export
print.sfcqr <- function(x, ...) {
  cat("Supervised functional censored quantile regression\n")
  cat("tau:", x$tau, " | observations:", x$n, " | components:", x$M,
      " | lambda:", x$lambda, "\n")
  criterion <- if (is.null(x$criterion)) "IBIC" else x$criterion
  cat("Selected ", criterion, ": ",
      format(x$criteria[[criterion]], digits = 6), "\n", sep = "")
  invisible(x)
}

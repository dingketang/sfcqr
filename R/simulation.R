#' Generate a small illustrative censored image dataset
#'
#' Creates smooth orthonormal image basis functions and a log-linear event-time
#' model. Censoring times are independent lognormal draws; their location is
#' calibrated against the simulated event times to the requested censoring
#' fraction. The realized fraction is random. This is a new demonstration
#' generator, distinct from [simulate_sfcqr_main()] and [gene_data()], which
#' reproduce the supplied main simulation scenario.
#'
#' @param n Number of observations, at least 30.
#' @param image_dim Two image dimensions, each at least 2.
#' @param n_basis Number of image basis functions, at most the number of pixels.
#' @param n_scalar Number of scalar covariates (possibly zero).
#' @param censoring Desired approximate censoring proportion in `[0, 0.8]`.
#' @param seed Optional integer random seed. When supplied, the caller's random
#'   number state is restored on exit.
#' @return A list containing `Y`, `delta`, `X`, `Z`, `Bxy`, `image_dim`, latent
#'   `event_time`, `censoring_time`, and `truth` (log-model coefficients).
#' @export
#' @examples
#' d <- simulate_sfcqr_data(n = 100, seed = 2026)
#' dim(d$Z)
#' mean(d$delta == 0)
simulate_sfcqr_data <- function(n = 160, image_dim = c(8, 8), n_basis = 4,
                               n_scalar = 2, censoring = 0.2, seed = 2026) {
  .scalar_number(n, "n", lower = 30, upper = .Machine$integer.max)
  .scalar_number(n_basis, "n_basis", lower = 1, upper = .Machine$integer.max)
  .scalar_number(n_scalar, "n_scalar", lower = 0, upper = .Machine$integer.max)
  .scalar_number(censoring, "censoring", lower = 0, upper = 0.8)
  if (any(c(n, n_basis, n_scalar) != as.integer(c(n, n_basis, n_scalar)))) {
    stop("n, n_basis, and n_scalar must be integers.")
  }
  if (!is.numeric(image_dim) || length(image_dim) != 2L ||
      any(!is.finite(image_dim)) || any(image_dim < 2) ||
      any(image_dim != as.integer(image_dim))) {
    stop("image_dim must contain two integers >= 2.")
  }
  if (n_basis > prod(image_dim)) stop("n_basis exceeds the pixel count.")
  if (n_scalar + n_basis + 2 >= n) stop("n is too small for the requested dimensions.")
  if (!is.null(seed)) {
    .scalar_number(seed, "seed", lower = 0, upper = .Machine$integer.max)
    if (seed != as.integer(seed)) stop("seed must be an integer.")
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
  nr <- image_dim[1]
  nc <- image_dim[2]
  frequencies <- expand.grid(row = 0:(nr - 1), col = 0:(nc - 1))
  frequencies <- frequencies[order(frequencies$row + frequencies$col,
                                   frequencies$row, frequencies$col), ]
  Bxy <- vapply(seq_len(n_basis), function(j) {
    as.numeric(outer(cos(pi * ((seq_len(nr) - 0.5) / nr) * frequencies$row[j]),
                     cos(pi * ((seq_len(nc) - 0.5) / nc) * frequencies$col[j])))
  }, numeric(nr * nc))
  Bxy <- sweep(Bxy, 2, sqrt(colSums(Bxy^2)), "/")
  latent <- matrix(stats::rnorm(n * n_basis), n, n_basis)
  Z <- latent %*% t(Bxy) + matrix(stats::rnorm(n * nr * nc, sd = 0.03), n)
  X <- matrix(stats::rnorm(n * n_scalar), n, n_scalar)
  basis_coef <- c(0.6, -0.4, 0.25, rep(0, max(0, n_basis - 3)))[seq_len(n_basis)]
  scalar_coef <- rep(c(0.25, -0.2), length.out = n_scalar)
  intercept <- 1
  eta <- as.numeric(intercept + X %*% scalar_coef + (Z %*% Bxy) %*% basis_coef)
  event_time <- exp(eta + stats::rnorm(n, sd = 0.35))
  if (censoring == 0) censoring_time <- rep(Inf, n) else {
    log_event <- log(event_time)
    location <- stats::uniroot(function(mu) {
      mean(stats::pnorm(log_event, mean = mu, sd = 0.8)) - censoring
    }, interval = c(min(log_event) - 10, max(log_event) + 10))$root
    censoring_time <- exp(stats::rnorm(n, mean = location, sd = 0.8))
  }
  list(Y = pmin(event_time, censoring_time),
       delta = as.integer(event_time <= censoring_time), X = X, Z = Z, Bxy = Bxy,
       image_dim = as.integer(image_dim), event_time = event_time,
       censoring_time = censoring_time,
       truth = list(intercept = intercept, coef_X = scalar_coef,
                    coef_estimate = basis_coef, error_sd = 0.35))
}

#' Run the main-scenario three-quantile simulation workflow
#' @param seed Integer simulation seed.
#' @param n Number of simulated observations.
#' @param taus Quantiles to fit; defaults to 0.3, 0.5, and 0.7.
#' @param Mmax Maximum component count.
#' @param criterion Integrated criterion used for model selection.
#' @param beran_covariates Covariates used in Beran censoring estimation.
#' @param ... Additional fitting arguments passed to [Sfcqr()].
#' @return A list of legacy coefficient vectors, one per quantile.
#' @export
simu_one_SFCQR <- function(seed, n = 500, taus = c(0.3, 0.5, 0.7), Mmax = 30,
                           criterion = "IGACV", beran_covariates = "all", ...) {
  data <- gene_data(seed, n = n)
  lapply(taus, function(tau) Sfcqr(data, tau, Mmax = Mmax,
    criterion = criterion, beran_covariates = beran_covariates, ...))
}

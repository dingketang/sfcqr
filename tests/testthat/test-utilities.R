test_that("Beran survival handles ties and preserves evaluation order", {
  Y <- c(1, 2, 2, 3)
  delta <- c(1, 0, 1, 1)
  result <- beran_Ghat_curve(Y, delta, rep(0, 4), x0 = 0,
                             times = c(3, 1, 2, 0), h = 1)
  expect_equal(result, c(2 / 3, 1, 2 / 3, 1))
  expect_equal(beran_Ghat(Y, delta, matrix(0, 4, 1), 0, 2, 1), 2 / 3)
  expect_length(beran_Ghat_curve(Y, delta, rep(0, 4), 0, numeric(0), 1), 0)
})

test_that("Beran kernel scaling prevents complete weight underflow", {
  result <- beran_Ghat_curve(c(1, 2, 3), c(1, 1, 0), c(0, 1, 2),
                             x0 = 10000, times = c(0, 2, 3), h = 1e-200)
  expect_equal(result, c(1, 1, 0))
  expect_true(all(is.finite(result)))
})

test_that("marginal IPCW keeps observation order and exposes the last event rule", {
  empty <- matrix(numeric(0), 4, 0)
  weights <- get_weights(c(3, 1, 4, 2), empty, c(1, 1, 1, 0), method = "marginal")
  expect_equal(weights, c(1.5, 1, 1.5, 0))
  Y <- c(1, 2, 3)
  delta <- c(1, 1, 0)
  X <- matrix(numeric(0), 3, 0)
  expect_equal(get_weights(Y, X, delta, force_last_event = TRUE), c(1, 1, 0))
  expect_equal(get_weights(Y, X, delta, force_last_event = FALSE), c(1, 1, 0))
  expect_equal(delta, c(1, 1, 0))
})

test_that("all supported censoring models return finite weights", {
  set.seed(941)
  n <- 150
  X <- matrix(stats::rnorm(n * 3), nrow = n)
  T <- exp(1 + 0.2 * X[, 1] + stats::rnorm(n))
  C <- exp(1 + 0.1 * X[, 2] + stats::rnorm(n))
  Y <- pmin(T, C)
  delta <- as.numeric(T <= C)
  for (method in c("marginal", "lognormal", "loglogistic", "beran", "cox")) {
    weights <- get_weights(Y, X, delta, method = method)
    expect_length(weights, n)
    expect_true(all(is.finite(weights)))
    expect_true(all(weights >= 0))
    expect_equal(get_weights(Y, X, rep(1, n), method = method), rep(1, n))
  }
  expect_equal(get_weights(Y, X, delta, method = "bernan"),
               get_weights(Y, X, delta, method = "beran"))
})

test_that("image roughness is a symmetric PSD Gram matrix on rectangular images", {
  set.seed(12)
  basis <- list(matrix(1, 5, 9), matrix(stats::rnorm(45), 5, 9),
                matrix(stats::rnorm(45), 5, 9))
  P <- get_P_matrix(basis)
  expect_equal(dim(P), c(3L, 3L))
  expect_equal(P, t(P))
  expect_true(min(eigen(P, symmetric = TRUE, only.values = TRUE)$values) >= -1e-10)
  expect_equal(P[1, ], rep(0, 3))
  expect_equal(get_P_matrix(lapply(basis, function(b) 2 * b)), 4 * P)
  expect_error(get_P_matrix(list(matrix(1, 3, 4), matrix(1, 4, 3))), "same dimensions")
})

test_that("matrix inverses use SVD only when needed", {
  expect_equal(Minverse(diag(c(2, 4))), diag(c(0.5, 0.25)))
  X <- matrix(c(1, 0, 0, 1, 1, 1), nrow = 3)
  expect_equal(Minverse(X) %*% X, diag(2), tolerance = 1e-12)
  singular <- matrix(c(1, 2, 2, 4), 2, 2)
  expect_warning(inverse <- Minverse(singular), "effective rank 1")
  expect_equal(singular %*% inverse %*% singular, singular, tolerance = 1e-12)
})

test_that("GACV numerators agree with direct estimating-equation sums", {
  Y <- seq(0.8, 3, length.out = 30)
  delta <- rep(c(1, 1, 0), 10)
  X <- matrix(seq(-1, 1, length.out = 30), ncol = 1)
  tau <- seq(0.05, 0.8, length.out = 12)
  coefficients <- rbind(log(1 + tau), rep(0.1, length(tau)))
  result <- gacv_fcrq(Y, delta, X, tau, coefficients)
  eta <- cbind(1, X) %*% coefficients
  qhat <- exp(eta)
  dH <- H_tau(tau) - H_tau(c(0, head(tau, -1)))
  numerator <- numeric(length(tau))
  for (j in seq_along(tau)) {
    for (i in seq_along(Y)) {
      integral <- 0
      for (k in seq_len(j)) {
        previous_quantile <- if (k == 1) 0 else qhat[i, k - 1]
        integral <- integral + (Y[i] >= previous_quantile) * dH[k]
      }
      numerator[j] <- numerator[j] - (log(Y[i]) - eta[i, j]) *
        (delta[i] * (Y[i] <= qhat[i, j]) - integral)
    }
  }
  expect_equal(result$pointwise_numerators, numerator, tolerance = 1e-12)
  expect_equal(result$IAIC, mean(numerator / (30 - 1)))
  expect_equal(result$IBIC, mean(numerator / (30 - log(30))))
  expect_equal(result$IGACV, mean(numerator / (30 - sqrt(log(30)))))
})

test_that("GACV invalid denominators cannot improve a tuning criterion", {
  Y <- seq(1, 4, length.out = 30)
  tau <- seq(0.05, 0.8, length.out = 12)
  coefficients <- matrix(log(stats::quantile(Y, tau)), nrow = 1)
  X <- matrix(numeric(0), 30, 0)
  expect_warning(result <- gacv_fcrq(Y, rep(1, 30), X, tau, coefficients,
                                     df_penalty = 15), "BIC")
  expect_true(is.infinite(result$IBIC))
  expect_true(is.finite(result$IAIC))
  expect_true(is.finite(result$IGACV))
  expect_warning(result <- gacv_fcrq(Y, rep(1, 30), X, tau, coefficients,
                                     df_penalty = 1000), "denominator")
  expect_true(all(is.infinite(c(result$IAIC, result$IBIC, result$IGACV))))
})

test_that("GACV filtering and single-column dimensions are explicit", {
  tau <- seq(0.05, 0.8, length.out = 12)
  result <- gacv_fcrq(1:20, rep(1, 20), matrix(numeric(0), 20, 0),
                      tau, matrix(Inf, 1, 12))
  expect_equal(result$n_tau, 0L)
  expect_equal(result$removed_tau_indices, seq_len(12))
  expect_true(all(is.infinite(c(result$IAIC, result$IBIC, result$IGACV))))
  expect_match(result$message, "at least 11")
  one <- gacv_fcrq(1, 1, matrix(numeric(0), 1, 0), 0.5, matrix(0, 1, 1), min_tau = 1)
  expect_equal(dim(one$cum_integral_mat), c(1L, 1L))
  expect_equal(dim(one$eta_mat), c(1L, 1L))
  expect_equal(one$IBIC, 0)
})

test_that("public utilities reject invalid input rather than recycle it", {
  expect_error(vec_to_image(1:3, 2, 2), "length")
  expect_error(get_weights(c(1, NA), matrix(1, 2, 1), c(1, 0)), "finite")
  expect_error(get_weights(c(1, 2), matrix(1, 2, 1), c(1, 2)), "0 or 1")
  expect_error(beran_Ghat(1:3, c(1, 0, 1), 1:3, 1, 2, 0), "positive")
  expect_error(cq_cov(1:3, 1:3, 2, c(1, 1), 0.5), "same length")
  expect_error(q_cov(1:3, 1:3, 1), "between 0 and 1")
})

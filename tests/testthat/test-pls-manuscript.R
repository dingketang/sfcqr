test_that("the optimizer reaches the independent analytic optimum", {
  x <- rbind(c(-3, -1), c(-1, -3), c(1, 3), c(3, 1))
  y <- 1:4
  delta <- rep(1, 4)
  p <- diag(c(0, 3))
  fit <- plsfit_cqcov(x, y, 1, rep(1, 4), delta,
                     lambda = 1, P = p, scale.X = FALSE)
  # F(2) = 1/2 and n^(-1) X' [tau - I(Y <= 2)] = (1, 1).
  q <- c(1, 1)
  metric <- diag(c(1, 4))
  expected <- c(4, 1) / sqrt(20)
  w <- as.numeric(fit$directions_raw[, 1])
  expect_equal(w, expected, tolerance = 1e-12)
  expect_equal(drop(crossprod(w, metric %*% w)), 1, tolerance = 1e-12)
  expect_equal(sum(q * w)^2, 1.25, tolerance = 1e-12)
  expect_equal(fit$Y_quantile, 2)
  # The previous inverse-square-root direction differs on this example.
  old_direction <- c(1, .5) / sqrt(1.25)
  expect_gt(max(abs(w - old_direction)), .1)
})

test_that("standardization preserves the raw constrained solution and deflation", {
  set.seed(701)
  x <- sweep(matrix(rnorm(240), 60, 4), 2, c(.05, 2, 8, .6), "*")
  y <- x[, 1] + x[, 2] + .2 * x[, 3] + rnorm(60)
  delta <- as.numeric(seq_len(60) %% 5 != 0)
  weights <- delta * seq(1, 3, length.out = 60)
  p <- crossprod(matrix(c(2, 1, 0, 0, 1, 2, 1, 0,
                          0, 1, 2, 1, 0, 0, 1, 2), 4))
  raw <- plsfit_cqcov(x, y, 3, weights, delta, lambda = 2, P = p,
                     scale.X = FALSE, tau = .6)
  scaled <- plsfit_cqcov(x, y, 3, weights, delta, lambda = 2, P = p,
                        scale.X = TRUE, tau = .6)
  expect_equal(raw$W, scaled$W, tolerance = 1e-9)
  expect_equal(raw$U, scaled$U, tolerance = 1e-9)
  expect_equal(raw$Y_quantile, scaled$Y_quantile)
  expect_equal(unname(sweep(x, 2, scaled$Xmean, "-") %*% scaled$W),
               unname(scaled$U), tolerance = 1e-9)
  gram <- crossprod(scaled$U)
  expect_lt(max(abs(gram[row(gram) != col(gram)])) / max(diag(gram)), 1e-10)
  direction_norms <- diag(crossprod(scaled$directions_raw,
                                   scaled$metric %*% scaled$directions_raw))
  expect_equal(unname(direction_norms), rep(1, 3), tolerance = 1e-9)
  scaled_norms <- diag(crossprod(scaled$directions,
                                scaled$metric_scaled %*% scaled$directions))
  expect_equal(unname(scaled_norms), rep(1, 3), tolerance = 1e-9)
})

test_that("IPCW quantiles respect ties, CDF jumps and original censoring", {
  x <- cbind(c(-2, -1, 0, 1, 3), c(0, 1, 0, 2, -1))
  y <- c(1, 2, 2, 3, 10)
  delta <- c(1, 1, 1, 1, 0)
  fit <- plsfit_cqcov(x, y, 1, c(1, 2, 1, 1, 1e7), delta, tau = .6)
  reference <- plsfit_cqcov(x, y, 1, c(1, 2, 1, 1, 0), delta, tau = .6)
  expect_equal(fit$Y_quantile, 2)
  expect_equal(fit$W, reference$W)
  expect_equal(fit$ipcw_weights, c(1, 2, 1, 1, 0))
  expect_equal(.ipcw_marginal_quantile(1:4, c(1, 3, 2, 2), .5), 2)
  expect_equal(.ipcw_marginal_quantile(y, c(1, 2, 1, 1, 0), .9), 3)
  stable <- plsfit_cqcov(x, y, 1, c(1, 2, 1, 1, 0) * 1e300,
                        delta, tau = .6)
  expect_equal(stable$Y_quantile, fit$Y_quantile)
  expect_equal(stable$W, fit$W, tolerance = 1e-12)
})

test_that("undefined IPCW moments and invalid weights are rejected", {
  x <- cbind(1:4, c(3, 1, 4, 2))
  expect_error(plsfit_cqcov(x, 1:4, 1, rep(1, 4), rep(0, 4)),
               "observed event")
  expect_error(plsfit_cqcov(x, 1:4, 1, rep(0, 4), rep(1, 4)),
               "observed event")
  expect_error(plsfit_cqcov(x, 1:4, 1, c(1, -1, 1, 1), rep(1, 4)),
               "nonnegative")
  expect_error(plsfit_cqcov(x, 1:4, 1, c(1, Inf, 1, 1), rep(1, 4)),
               "finite")
})

test_that("a large penalty does not discard a well-defined direction", {
  x <- rbind(c(-3, -1), c(-1, -3), c(1, 3), c(3, 1))
  fit <- plsfit_cqcov(x, 1:4, 1, rep(1, 4), rep(1, 4),
                     lambda = 1e20, P = diag(2), scale.X = FALSE)
  expect_equal(fit$ncomp, 1L)
  norm <- drop(crossprod(fit$directions_raw,
                        fit$metric %*% fit$directions_raw))
  expect_equal(norm, 1, tolerance = 1e-12)
  expect_true(all(is.finite(fit$U)))
  expect_gt(max(abs(fit$U)), 0)
})

test_that("the model defaults use manuscript tuning and event indicators", {
  d <- simulate_sfcqr_data(n = 100, censoring = .15, seed = 206)
  fit <- sfcqr(d, Mmax = 2, lambda_list = c(0, 1))
  explicit <- sfcqr(d, Mmax = 2, lambda_list = c(0, 1),
                    criterion = "IGACV", beran_covariates = "all",
                    force_last_event = FALSE)
  expect_equal(fit$criterion, "IGACV")
  expect_equal(fit$beran_covariates, "all")
  expect_false(fit$force_last_event)
  expect_equal(fit$weights[d$delta == 0], rep(0, sum(d$delta == 0)))
  expect_equal(predict(fit), predict(explicit), tolerance = 1e-10)
})


test_that("numerical cancellation cannot create extra supervised components", {
  set.seed(740)
  v <- c(rep(-1, 4), rep(1, 4))
  q <- qr.Q(qr(cbind(rep(1, 8), v, matrix(rnorm(32), 8, 4))))
  x <- q[, 2:5]
  y <- as.numeric(v > 0)
  for (scaled in c(TRUE, FALSE)) {
    fit <- plsfit_cqcov(x, y, 4, rep(1, 8), rep(1, 8), scale.X = scaled)
    expect_equal(fit$ncomp, 1L)
    expect_equal(unname(sweep(x, 2, fit$Xmean, "-") %*% fit$W),
                 unname(fit$U), tolerance = 1e-12)
  }
})

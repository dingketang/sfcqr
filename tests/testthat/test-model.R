test_that("README fit and new-data prediction reconstruct the crq scores", {
  d <- simulate_sfcqr_data(seed = 2026)
  fit <- sfcqr(d, Mmax = 3, lambda_list = c(0, 1))
  expect_s3_class(fit, "sfcqr")
  expect_true(all(is.finite(coef(fit))))
  b <- as.numeric(stats::coef(fit$crqfit, taus = .5))
  score_eta <- as.numeric(cbind(1, fit$scores, d$X) %*% b)
  expect_equal(predict(fit, type = "link"), score_eta, tolerance = 1e-8)
  new <- list(X = d$X[1:5, , drop = FALSE], Z = d$Z[1:5, , drop = FALSE])
  pred <- predict(fit, new)
  expect_equal(pred, predict(fit)[1:5], tolerance = 1e-10)
  image <- coef(fit, type = "image")
  pixel_eta <- as.numeric(fit$intercept + new$X %*% fit$coef_X +
                           new$Z %*% as.numeric(image))
  expect_equal(log(pred), pixel_eta, tolerance = 1e-8)
  expect_true(all(pred > 0))
  expect_equal(dim(image), c(8L, 8L))
})

test_that("zero and one scalar covariate and rectangular images fit", {
  for (p in 0:1) {
    d <- simulate_sfcqr_data(n = 130, image_dim = c(6, 8), n_scalar = p,
                            censoring = .1, seed = 40 + p)
    fit <- sfcqr(d, Mmax = 2, lambda_list = c(0, 1))
    expect_length(fit$coef_X, p)
    expect_equal(dim(coef(fit, type = "image")), c(6L, 8L))
    expect_equal(predict(fit, list(X = d$X, Z = d$Z)), predict(fit))
    if (p == 0L) expect_equal(predict(fit, list(Z = d$Z)), predict(fit))
  }
})

test_that("identical images are rejected before basis projection", {
  d <- simulate_sfcqr_data(n = 60, censoring = 0, seed = 1)
  for (value in c(0, 1)) {
    d$Z[] <- value
    expect_error(sfcqr(d, Mmax = 1, lambda_list = 0),
                 "no between-subject variation")
  }
  # A shared spatial pattern also has no between-subject information.
  pattern <- seq_len(ncol(d$Z)) / ncol(d$Z)
  d$Z <- matrix(rep(pattern, each = nrow(d$Z)), nrow = nrow(d$Z))
  expect_error(sfcqr(d, Mmax = 1, lambda_list = 0),
               "no between-subject variation")
})

test_that("integrated tuning explains an insufficient quantile grid", {
  d <- simulate_sfcqr_data(n = 100, censoring = 0, seed = 1)
  expect_error(
    sfcqr(d, Mmax = 1, lambda_list = 0,
          grid = c(.1, .25, .5, .75, .9)),
    "at least 11 grid levels at or below tau"
  )
  expect_error(
    sfcqr(d, tau = .05, Mmax = 1, lambda_list = 0),
    "at least 11 grid levels at or below tau"
  )
})

test_that("a dense supplied grid supports a small target quantile", {
  d <- simulate_sfcqr_data(n = 100, censoring = 0, seed = 1)
  low_grid <- sort(unique(c(seq(.001, .05, length.out = 21),
                            seq(.06, .99, by = .01))))
  fit <- sfcqr(d, tau = .05, Mmax = 1, lambda_list = 0, grid = low_grid)
  expect_s3_class(fit, "sfcqr")
  expect_true(all(is.finite(coef(fit))))
  expect_true(is.finite(fit$criteria$IGACV))
  expect_gte(fit$criteria$n_tau, 11L)
  expect_true(all(fit$criteria$tau_vec <= .05))
})

test_that("input validation explains malformed data", {
  d <- simulate_sfcqr_data(n = 100, seed = 7)
  expect_error(sfcqr(d, tau = 1), "tau")
  d$Y[1] <- 0
  expect_error(sfcqr(d), "positive")
  d$Y[1] <- 1
  d$X <- cbind(d$X[, 1], d$X[, 1])
  expect_error(sfcqr(d), "linearly independent")
})

test_that("legacy wrapper returns a dynamically sized vector", {
  d <- simulate_sfcqr_data(n = 100, n_scalar = 1, censoring = 0, seed = 4)
  vec <- Sfcqr(d, .5, Mmax = 1, lambda_list = 0)
  expect_length(vec, 2 + 1 + 4)
  expect_equal(names(vec)[1:3], c("M_num", "X0", "X1"))
})

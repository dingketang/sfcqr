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

test_that("unusable candidates fail instead of choosing an infinite score", {
  d <- simulate_sfcqr_data(n = 60, censoring = 0, seed = 1)
  expect_error(sfcqr(d, tau = .05, Mmax = 1, lambda_list = 0),
               "No valid candidate")
  d$Z[] <- 1
  expect_error(sfcqr(d, Mmax = 1, lambda_list = 0), "No valid candidate")
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

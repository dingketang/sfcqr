test_that("deflation projections reproduce orthogonal scores on raw inputs", {
  d <- simulate_sfcqr_data(n = 130, n_basis = 4, censoring = 0, seed = 18)
  x <- d$Z %*% d$Bxy
  A <- matrix(c(2, 1, 0, 0, 1, 2, 1, 0, 0, 1, 2, 1, 0, 0, 1, 2), 4)
  P <- crossprod(A)
  for (scaled in c(TRUE, FALSE)) {
    fit <- plsfit_cqcov(x, log(d$Y), 3, rep(1, nrow(x)), d$delta,
                       lambda = 2, P = P, scale.X = scaled)
    expect_equal(ncol(fit$U), 3L)
    projected <- sweep(x, 2, fit$Xmean, "-") %*% fit$W
    expect_equal(unname(projected), unname(fit$U), tolerance = 1e-9)
    gram <- crossprod(fit$U)
    expect_lt(max(abs(gram[row(gram) != col(gram)])), 1e-8)
  }
})

test_that("unpenalized first component preserves quantile covariance extraction", {
  d <- simulate_sfcqr_data(n = 100, censoring = 0, seed = 9)
  x <- d$Z %*% d$Bxy
  fit <- plsfit_cqcov(x, log(d$Y), 2, rep(1, 100), d$delta)
  x0 <- sweep(x, 2, fit$Xmean, "-")
  manual <- as.numeric(cq_cov(x0, log(d$Y), fit$Y_quantile, rep(1, 100), .5))
  manual <- manual / sqrt(sum(manual^2))
  expect_equal(unname(fit$W[, 1]), manual, tolerance = 1e-10)
})

test_that("rank truncation and scalar penalization are safe", {
  d <- simulate_sfcqr_data(n = 100, censoring = 0, seed = 3)
  x <- d$Z %*% d$Bxy
  duplicate <- cbind(x[, 1], x[, 1])
  fit <- plsfit_cqcov(duplicate, log(d$Y), 3, rep(1, 100), d$delta)
  expect_equal(fit$ncomp, 1L)
  one <- plsfit_cqcov(x[, 1, drop = FALSE], log(d$Y), 1,
    rep(1, 100), d$delta, lambda = 1, P = matrix(1, 1, 1))
  expect_true(all(is.finite(one$W)))
  expect_error(plsfit_cqcov(x, log(d$Y), 1, rep(1, 100), d$delta,
                           lambda = 1, P = c(0, 0, 1)), "square matrix")
})

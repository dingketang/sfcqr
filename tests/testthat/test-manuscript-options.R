test_that("multivariate Beran uses the Gaussian product kernel", {
  Y <- c(1, 2, 3, 4)
  delta <- c(1, 0, 1, 1)
  X <- cbind(rep(0, 4), c(-1, 0, 1, 2))
  x0 <- c(0, 1)
  h <- .7
  kernel <- apply(sweep(X, 2, x0, "-") / h, 1, function(row) {
    prod(stats::dnorm(row))
  })
  expected <- 1 - kernel[2] / sum(kernel[Y >= 2])
  expect_equal(beran_Ghat(Y, delta, X, x0, 3, h), expected)
  all_weights <- get_weights(Y, X, delta, "beran", bandwidth = h,
                            beran_covariates = "all", force_last_event = FALSE)
  first_weights <- get_weights(Y, X, delta, "beran", bandwidth = h,
                              beran_covariates = "first", force_last_event = FALSE)
  expect_equal(all_weights[3], 1 / expected)
  expect_false(isTRUE(all.equal(all_weights, first_weights)))
  expect_error(beran_Ghat(Y, delta, X, 0, 3, h), "one value")
})

test_that("model selection uses the explicitly requested integrated criterion", {
  d <- simulate_sfcqr_data(n = 110, censoring = .1, seed = 81)
  for (criterion in c("IGACV", "IAIC", "IBIC")) {
    fit <- sfcqr(d, Mmax = 2, lambda_list = c(0, 1), criterion = criterion,
                 beran_covariates = "all")
    expect_equal(fit$criterion, criterion)
    expect_equal(fit$beran_covariates, "all")
    expect_equal(min(fit$tuning), unname(fit$criteria[[criterion]]))
    expect_equal(fit$weights, get_weights(d$Y, d$X, d$delta, "beran",
                                         beran_covariates = "all"))
  }
})

test_that("IPCW retains observed indicators when the maximum time is censored", {
  Y <- c(1, 2, 3, 4)
  delta <- c(1, 0, 1, 0)
  delta_before <- delta
  for (method in c("marginal", "beran")) {
    X <- matrix(0, nrow = length(Y), ncol = 1L)
    default <- get_weights(Y, X, delta, method, bandwidth = 1)
    explicit <- get_weights(Y, X, delta, method, bandwidth = 1,
                            force_last_event = FALSE)
    expect_equal(default, c(1, 0, 1.5, 0))
    expect_equal(default, explicit)
    forced <- get_weights(Y, X, delta, method, bandwidth = 1,
                          force_last_event = TRUE)
    expect_equal(forced[delta == 0], c(0, 0))
    expect_equal(forced, c(1, 0, 1.5, 0))
  }
  expect_identical(delta, delta_before)
})

test_that("a tail-only censoring override leaves its supervised weight zero", {
  Y <- c(1, 2, 3, 4)
  X <- c(-1, 0, 1, 2)
  delta <- c(1, 1, 1, 0)
  for (method in c("marginal", "beran", "lognormal", "loglogistic", "cox")) {
    # The override removes all censoring events in the G fit, where G = 1.
    expect_equal(get_weights(Y, X, delta, method, force_last_event = TRUE),
                 c(1, 1, 1, 0))
  }
  for (method in c("marginal", "beran")) {
    expect_equal(get_weights(Y, X, delta, method, bandwidth = 1),
                 c(1, 1, 1, 0))
  }
})

test_that("the optional override changes G only, including at tied tail times", {
  Y <- c(1, 2, 4, 4)
  X <- matrix(0, nrow = length(Y), ncol = 2L)
  delta <- c(1, 0, 0, 1)
  for (method in c("marginal", "beran")) {
    observed <- get_weights(Y, X, delta, method, bandwidth = 1,
                            beran_covariates = "all")
    forced <- get_weights(Y, X, delta, method, bandwidth = 1,
                          force_last_event = TRUE, beran_covariates = "all")
    # At time 2, G = 2/3. Without the override, censoring at time 4
    # halves it again to 1/3; the original event at time 4 has weight 3.
    expect_equal(observed, c(1, 0, 0, 3))
    expect_equal(forced, c(1, 0, 0, 1.5))
    expect_equal(forced[delta == 0], c(0, 0))
  }
  expect_identical(delta, c(1, 0, 0, 1))
})

test_that("all observed and all censored samples have exact IPCW limits", {
  Y <- c(1, 2, 3, 4)
  X <- c(-1, 0, 1, 2)
  for (method in c("marginal", "beran", "lognormal", "loglogistic", "cox")) {
    for (force in c(FALSE, TRUE)) {
      expect_equal(get_weights(Y, X, rep(1, 4), method,
                               force_last_event = force), rep(1, 4))
      expect_equal(get_weights(Y, X, rep(0, 4), method,
                               force_last_event = force), rep(0, 4))
    }
  }
})

test_that("zero covariates select the marginal model with original numerators", {
  Y <- c(1, 2, 4, 4)
  delta <- c(1, 0, 0, 1)
  X <- matrix(numeric(0), nrow = length(Y), ncol = 0L)
  expect_equal(get_weights(Y, X, delta, method = "beran"), c(1, 0, 0, 3))
  expect_equal(get_weights(Y, X, delta, method = "beran",
                           force_last_event = TRUE), c(1, 0, 0, 1.5))
})

test_that("Beran IPCW conditions on every supplied covariate by default", {
  Y <- c(1, 2, 3, 4)
  delta <- c(1, 0, 1, 1)
  X <- cbind(rep(0, 4), c(-1, 0, 1, 2))
  h <- .7
  w <- get_weights(Y, X, delta, "beran", bandwidth = h)
  kernel_at_third <- exp(-.5 * ((X[, 2] - X[3, 2]) / h)^2)
  G_at_third <- 1 - kernel_at_third[2] / sum(kernel_at_third[Y >= 2])
  expect_equal(w[3], 1 / G_at_third)
  expect_equal(w, get_weights(Y, X, delta, "beran", bandwidth = h,
                              beran_covariates = "all"))
  expect_false(isTRUE(all.equal(w, get_weights(Y, X, delta, "beran",
                    bandwidth = h, beran_covariates = "first"))))
})

test_that("fitted censoring models retain zero weights with the optional override", {
  set.seed(614)
  n <- 120L
  X <- cbind(stats::rnorm(n), stats::rnorm(n))
  event_time <- exp(.2 * X[, 1] + stats::rnorm(n))
  censor_time <- exp(.2 * X[, 2] + stats::rnorm(n))
  Y <- pmin(event_time, censor_time)
  delta <- as.numeric(event_time <= censor_time)
  # Ensure that the separate tail convention is exercised in each model.
  delta[which.max(Y)] <- 0
  for (method in c("marginal", "beran", "lognormal", "loglogistic", "cox")) {
    for (force in c(FALSE, TRUE)) {
      w <- get_weights(Y, X, delta, method, force_last_event = force)
      expect_equal(w[delta == 0], rep(0, sum(delta == 0)))
      expect_true(all(is.finite(w[delta == 1]) & w[delta == 1] >= 1))
    }
  }
})

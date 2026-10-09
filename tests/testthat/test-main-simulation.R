test_that("main simulation equals the supplied research generator", {
  # Literal main-branch reference protects the research random draw order,
  # tensor ordering, normalization, pixel inner products, and censoring law.
  reference <- function(seed, n) {
    set.seed(seed)
    X1 <- runif(n)
    X2 <- runif(n)
    B <- splines::bs(1:100, df = 6, degree = 3, intercept = TRUE)
    Borth <- qr.Q(qr(B))
    Bxy <- matrix(0, 10000, 36)
    for (i in 0:5) {
      for (j in 0:5) {
        Bxy[, i * 6 + j + 1] <- c(matrix(Borth[, i + 1]) %*%
                                   t(matrix(Borth[, j + 1])))
      }
    }
    Bxy <- Bxy %*% diag(1 / sqrt(diag(t(Bxy) %*% Bxy)))
    A <- matrix(runif(n * 36) * 4, ncol = 36, nrow = n)
    coefz <- diag(sqrt(sqrt(1 / (1:36))))
    Z <- A %*% coefz %*% t(Bxy)
    C1 <- Bxy[, 1] + Bxy[, 8]
    C2 <- 0.8 * Bxy[, 22] + 1.2 * Bxy[, 29]
    epsilon <- runif(n) * 2 - 1
    var_epsilon <- var(epsilon)
    logTi <- Z %*% C2 + Z %*% C1 * epsilon + X1 * epsilon + X2 * 0.5
    Ti <- pmax(exp(logTi), 1e-10)
    u <- rlogis(n, location = 0, scale = 4)
    C <- exp(4 * X1 + 4 * X2 + u - 2)
    delta <- C > Ti
    Y <- Ti
    Y[!delta] <- C[!delta]
    list(Y = Y, X = cbind(X1, X2), delta = delta, Z = Z, Bxy = Bxy,
         Ti = Ti, A = A, var_epsilon = var_epsilon)
  }
  expected <- reference(17, 40)
  d <- simulate_sfcqr_main(n = 40, seed = 17)
  expect_identical(d[names(expected)], expected)
  expect_identical(gene_data(seed = 17, n = 40), d)
})

test_that("main data reconstructs the image and conditional event-time law", {
  d <- simulate_sfcqr_main(n = 40, seed = 17)
  expect_identical(dim(d$X), c(40L, 2L))
  expect_identical(dim(d$Z), c(40L, 10000L))
  expect_identical(dim(d$Bxy), c(10000L, 36L))
  expect_identical(d$image_dim, c(100L, 100L))
  expect_equal(crossprod(d$Bxy), diag(36), tolerance = 1e-12)
  expect_equal(d$Z %*% d$Bxy, d$latent_scores, tolerance = 1e-11)
  expect_equal(which(d$truth$C1 != 0), c(B1 = 1L, B8 = 8L))
  expect_equal(unname(d$truth$C1[c(1, 8)]), c(1, 1))
  expect_equal(which(d$truth$C2 != 0), c(B22 = 22L, B29 = 29L))
  expect_equal(unname(d$truth$C2[c(22, 29)]), c(0.8, 1.2))
  reconstructed_log_T <- as.numeric(((d$Z %*% d$Bxy) %*% d$truth$C1 +
                                      d$X[, 1]) * d$error +
                                     (d$Z %*% d$Bxy) %*% d$truth$C2 +
                                     0.5 * d$X[, 2])
  expect_equal(d$log_event_time, reconstructed_log_T, tolerance = 1e-11)
  expect_identical(d$event_time, d$Ti)
  expect_equal(as.numeric(d$event_time), pmax(exp(d$log_event_time), 1e-10))
  expect_equal(d$censoring_time, exp(d$log_censoring_time))
  expect_equal(d$Y, pmin(d$event_time, d$censoring_time))
  expect_identical(d$delta, d$censoring_time > d$event_time)
  expect_identical(d$var_epsilon, var(d$error))
  expect_true(all(d$X >= 0 & d$X <= 1))
  expect_true(all(d$error >= -1 & d$error <= 1))
  expect_true(all(d$A >= 0 & d$A <= 4))
})

test_that("main basis retains the supplied spline and tensor ordering", {
  d <- simulate_sfcqr_main(n = 30, seed = 2)
  Q <- qr.Q(qr(splines::bs(1:100, df = 6, degree = 3, intercept = TRUE)))
  expected_B2 <- as.numeric(outer(Q[, 1], Q[, 2]))
  expected_B7 <- as.numeric(outer(Q[, 2], Q[, 1]))
  expect_equal(d$Bxy[, 2], expected_B2 / sqrt(sum(expected_B2^2)))
  expect_equal(d$Bxy[, 7], expected_B7 / sqrt(sum(expected_B7^2)))
  expect_equal(d$setting$tensor_index,
               data.frame(row = rep(1:6, each = 6), col = rep(1:6, 6)))
  expect_equal(unname(d$setting$knots), c(34, 67))
  expect_equal(d$setting$boundary_knots, c(1, 100))
  expect_identical(d$setting$censor_shift, -2)
})

test_that("main conditional quantiles have the stated coefficient truth", {
  d <- simulate_sfcqr_main(n = 30, seed = 11)
  expect_equal(d$truth$taus, c(0.3, 0.5, 0.7))
  expect_equal(unname(d$truth$coef_X_tau), rbind(c(-0.4, 0, 0.4), rep(0.5, 3)))
  expect_equal(unname(d$truth$intercept_tau), rep(0, 3))
  scores <- d$Z %*% d$Bxy
  for (j in seq_along(d$truth$taus)) {
    tau <- d$truth$taus[j]
    expected_quantile <- as.numeric((scores %*% d$truth$C1 + d$X[, 1]) *
      (2 * tau - 1) + scores %*% d$truth$C2 + 0.5 * d$X[, 2])
    coefficient_quantile <- as.numeric(scores %*% d$truth$coef_tau[, j] +
      d$X %*% d$truth$coef_X_tau[, j])
    expect_equal(coefficient_quantile, expected_quantile, tolerance = 1e-12)
  }
})

test_that("main simulation preserves RNG and allows a fixed censoring shift", {
  set.seed(42)
  previous_rng <- .Random.seed
  a <- simulate_sfcqr_main(n = 30, seed = 19)
  expect_identical(.Random.seed, previous_rng)
  b <- simulate_sfcqr_main(n = 30, seed = 19)
  expect_identical(a, b)
  shifted <- simulate_sfcqr_main(n = 30, seed = 19, censor_shift = 1)
  expect_identical(shifted$event_time, a$event_time)
  expect_equal(shifted$log_censoring_time, a$log_censoring_time + 3)
  expect_identical(.Random.seed, previous_rng)
})

test_that("main simulation validates inputs", {
  expect_error(simulate_sfcqr_main(n = 29), "n")
  expect_error(simulate_sfcqr_main(n = 30.5), "integer")
  expect_error(simulate_sfcqr_main(seed = 1.5), "integer")
  expect_error(simulate_sfcqr_main(censor_shift = NA_real_), "finite")
})

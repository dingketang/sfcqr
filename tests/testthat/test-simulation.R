test_that("simulation is reproducible and preserves caller RNG", {
  set.seed(17)
  rng <- .Random.seed
  a <- simulate_sfcqr_data(n = 80, seed = 9)
  expect_identical(.Random.seed, rng)
  b <- simulate_sfcqr_data(n = 80, seed = 9)
  expect_identical(a, b)
  expect_equal(a$Y, pmin(a$event_time, a$censoring_time))
  expect_equal(a$delta, as.integer(a$event_time <= a$censoring_time))
  expect_equal(crossprod(a$Bxy), diag(4), tolerance = 1e-12)
})

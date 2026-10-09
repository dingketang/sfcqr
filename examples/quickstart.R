library(sfcqr)

dat <- simulate_sfcqr_main(n = 500, seed = 2026, censor_shift = -2)
taus <- c(0.3, 0.5, 0.7)

stopifnot(
  identical(dim(dat$X), c(500L, 2L)),
  identical(dim(dat$Z), c(500L, 10000L)),
  identical(dim(dat$Bxy), c(10000L, 36L)),
  all(dat$delta == as.integer(dat$censoring_time > dat$event_time)),
  isTRUE(all.equal(dat$Y, pmin(dat$event_time, dat$censoring_time)))
)
cat(sprintf("Observed censoring: %.1f%%\n", 100 * mean(dat$delta == 0)))

fits <- setNames(
  lapply(taus, function(tau) sfcqr(
    data = dat,
    tau = tau,
    Mmax = 30,
    lambda_list = c(0, 1, 10, 100, 1000, 1e6),
    image_dim = c(100, 100),
    criterion = "IGACV",
    censoring_method = "beran",
    beran_covariates = "all",
    bandwidth = nrow(dat$X)^(-1/6),
    force_last_event = FALSE
  )),
  paste0("tau=", taus)
)

newdat <- list(
  X = dat$X[1:5, , drop = FALSE],
  Z = dat$Z[1:5, , drop = FALSE]
)

for (fit in fits) {
  pred <- predict(fit, newdata = newdat, type = "response")
  reconstructed <- as.numeric(exp(
    fit$intercept +
      newdat$X %*% fit$coef_X +
      (newdat$Z %*% dat$Bxy) %*% fit$coef_estimate
  ))
  stopifnot(
    all(is.finite(coef(fit))),
    length(pred) == 5L,
    all(is.finite(pred)),
    all(pred > 0),
    isTRUE(all.equal(pred, reconstructed, tolerance = 1e-8))
  )
}

fit <- fits[["tau=0.5"]]
print(fit)
cat("Main-setting example passed at tau = 0.3, 0.5, and 0.7.\n")

# sfcqr

`sfcqr` is an R package for supervised functional censored quantile regression with right-censored outcomes, scalar covariates, and two-dimensional image predictors. The package provides documented fitting functions, coefficient reconstruction, prediction, reproducible simulations, tests, and GitHub Actions.

The model uses an image basis expansion, inverse censoring probability weights, penalized quantile covariance component extraction, and `quantreg::crq(..., method = "PengHuang")` for the final quantile regression. Version 0.1.2 defaults to IGACV tuning, Beran conditioning on all scalar covariates, and the original event indicators. The example below uses the main simulation mechanism and these settings.

## Associated manuscript

Tang, D., He, X., and Jiang, S. (2026). *Supervised functional censored quantile regression with application to mammogram data*. Manuscript under revision.

The section references below refer to this manuscript and its supplementary material. Use `citation("sfcqr")` to obtain the manuscript citation from R.

## Installation

Install the current development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("dingketang/sfcqr", ref = "main")
library(sfcqr)
```

For a fixed source snapshot, replace `"main"` with the full Git commit SHA of the version used for the analysis.

## Reproducible example: the main simulation setting

This example follows the main simulation in Section 5.1 of the manuscript and the Beran configuration in Supplementary Section S4. It uses the supplied main-branch `gene_data()` mechanism, with these settings:

| Setting | Value |
| --- | --- |
| Sample size | `n = 500`, one of the manuscript's sample sizes |
| Scalar covariates | Two independent `Unif(0, 1)` covariates |
| Image grid | `100 x 100` pixels |
| Working basis | Six cubic B-splines per spatial dimension, orthonormalized by QR; 36 tensor-product functions |
| Image scores | Independent `Unif(0, 4)` scores multiplied by `ell^(-1/4)` |
| Noise | `Unif(-1, 1)` |
| Censoring | Covariate-dependent logistic censoring with scale 4 and fixed shift -2, as in the supplied main branch |
| Target quantiles | `0.3`, `0.5`, and `0.7` |
| Component search | `Mmax = 30` |
| Penalty candidates | `c(0, 1, 10, 100, 1000, 1e6)` |
| Tuning criterion | IGACV |
| Censoring weights | Gaussian Beran estimator conditional on both scalar covariates; `h = n^(-1/6)` |

Writing `B_ell` for the tensor-product basis functions, the image and coefficient surfaces are

$$
Z_i(s) = \sum_{\ell=1}^{36} A_{i\ell}\ell^{-1/4}B_\ell(s),
\qquad A_{i\ell} \sim \mathrm{Unif}(0,4),
$$

$$
C_1(s) = B_1(s) + B_8(s),
\qquad C_2(s) = 0.8B_{22}(s) + 1.2B_{29}(s).
$$

The main-branch generator uses

$$
\log T_i = \{\langle Z_i,C_1\rangle + X_{i1}\}\epsilon_i
            + \langle Z_i,C_2\rangle + 0.5X_{i2},
\qquad \epsilon_i \sim \mathrm{Unif}(-1,1),
$$

$$
\log C_i = 4X_{i1} + 4X_{i2} + U_i - 2,
\qquad U_i \sim \mathrm{Logistic}(0,4).
$$

The implementation retains `T_i = pmax(exp(log(T_i)), 1e-10)`, the observed time `Y_i = min(T_i, C_i)`, and the original event indicator `C_i > T_i`. The fixed censoring shift targets the manuscript's approximately 50% censoring setting; the observed percentage varies across random samples.

The true quantile coefficients are `beta(tau) = c(2 * tau - 1, 0.5)` and `alpha(s, tau) = (2 * tau - 1) * C1(s) + C2(s)`. The generator returns these in `dat$truth` for comparisons with fitted coefficients.

```r
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
```

After installation, run the same example from the source repository:

```sh
Rscript examples/quickstart.R
```

A successful run ends with `Main-setting example passed at tau = 0.3, 0.5, and 0.7.` The checks verify the main-setting data dimensions, observed-data construction, finite fitted coefficients, positive predictions, and agreement between `predict()` and direct reconstruction from the original-scale coefficients.

Reference output for the example above under R 4.4.1:

```text
Observed censoring: 54.8%
Supervised functional censored quantile regression
tau: 0.5  | observations: 500  | components: 3  | lambda: 1
Selected IGACV: 0.400344
Main-setting example passed at tau = 0.3, 0.5, and 0.7.
```

This example demonstrates model fitting, coefficient reconstruction, and prediction for one dataset from the main simulation setting. The manuscript's Monte Carlo study also considers `n = 1000` and `n = 2000` with 1,000 replications.

For exact correspondence with the supplied generator, the one-dimensional basis is `splines::bs(1:100, df = 6, degree = 3, intercept = TRUE)`, followed by `qr.Q(qr(...))`. Tensor columns use the original ordering with the first basis index as the outer loop and the second as the inner loop; columns are normalized. Functional inner products are discrete pixel sums, matching the supplied code.

## Data format and model interface

The `data` argument is a list containing:

| Element | Type and dimensions | Meaning |
| --- | --- | --- |
| `Y` | Numeric vector of length `n` | Strictly positive observed time `min(T, C)` |
| `delta` | Binary vector of length `n` | `1` for an observed event; `0` for right censoring |
| `X` | Numeric `n x p` matrix | Scalar covariates; use `matrix(numeric(0), n, 0)` when there are none |
| `Z` | Numeric `n x pixels` matrix | One vectorized image per subject |
| `Bxy` | Numeric `pixels x q` matrix | Basis values in the same pixel order as `Z` |

To preserve the supplied generator's original output, `simulate_sfcqr_main()` returns `Y`, `Ti`, and `delta` as single-column matrices. The fitter accepts these directly, as well as ordinary response vectors.

Supply complete, finite data. The package does not impute missing values. The pixel ordering of `Z` and `Bxy` must agree: basis projection is `Z %*% Bxy`, with no additional integration weights. For externally supplied images, use the same vectorization convention for the image and basis; R's column-major matrix order is the package convention.

```r
sfcqr(
  data,
  tau = 0.5,
  Mmax = 30,
  lambda_list = c(0, 1, 10, 100, 1000, 1e6),
  image_dim = NULL,
  censoring_method = "beran",
  bandwidth = NULL,
  min_survival = 1e-4,
  force_last_event = FALSE,
  scale.X = TRUE,
  penalty = NULL,
  grid = NULL,
  criterion = "IGACV",
  beran_covariates = "all"
)
```

`tau` is the target quantile, `Mmax` is the maximum number of supervised components, and `lambda_list` contains the candidate smoothing penalties. `image_dim = c(nrow, ncol)` specifies the image grid; its product must equal `ncol(Z)`. A supplied `penalty` must be a `q x q` matrix. `scale.X` controls numerical standardization of basis projections during supervised component extraction; both settings preserve the constraint defined in the original basis coordinates. `grid` optionally supplies the quantile grid used to compare candidate models.

Integrated model selection requires at least 11 usable positive quantile levels at or below `tau`. Supplied and default grids are checked before candidate fitting. For a small target quantile, supply a denser grid below `tau` and retain at least one level above it. Censoring can further limit the estimable range.

For example, an uncensored low-quantile fit can use:

```r
d <- simulate_sfcqr_data(n = 100, censoring = 0, seed = 1)
low_grid <- sort(unique(c(
  seq(0.001, 0.05, length.out = 21),
  seq(0.06, 0.99, by = 0.01)
)))
low_fit <- sfcqr(d, tau = 0.05, Mmax = 1, lambda_list = 0, grid = low_grid)
```

`criterion` accepts `"IGACV"`, `"IAIC"`, or `"IBIC"`, with `"IGACV"` as the default. To use the earlier information-criterion choice, specify `criterion = "IBIC"` explicitly.

Censoring weights support `"beran"`, `"marginal"`, `"lognormal"`, `"loglogistic"`, and `"cox"`. For Beran weights, `beran_covariates = "all"` conditions on all columns of `X`, as specified in Supplementary Section S4. The default is `"all"`; `"first"` remains available as an explicit option. The default bandwidth is `n^(-1 / (p + 4))`, where `p = ncol(X)`; the example supplies `n^(-1/6)` explicitly because there are two covariates. With no scalar covariates, the package uses marginal weights. The original spelling `"bernan"` remains a compatibility alias.

`min_survival` bounds the estimated censoring survival probability below, and weights use the right-continuous survival probability at the observed time, `G(Y)`, including censoring events tied at that time. The default `force_last_event = FALSE` fits that distribution using the original event indicators. The numerator of each IPCW weight is always the original `delta`, so censored observations have zero supervision weight. Setting `force_last_event = TRUE` explicitly applies the legacy largest-time adjustment to the censoring-distribution fit only; it does not change the weight numerator or the final quantile regression event indicators.

## Inspecting a fitted model

The fitted object records the selected component count and smoothing parameter. Coefficient and prediction methods return results on the original input scale:

```r
fit$M
fit$lambda
alpha_hat <- coef(fit, type = "image")
dim(alpha_hat)
predict(fit, newdata = newdat, type = "response")
predict(fit, newdata = newdat, type = "link")
```

`type = "response"` returns the fitted event-time quantile; `type = "link"` returns its logarithm. See `?sfcqr`, `?coef.sfcqr`, and `?predict.sfcqr` for the available fields and arguments.

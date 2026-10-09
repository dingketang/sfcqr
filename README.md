# sfcqr

`sfcqr` is an R package for supervised functional censored quantile regression with right-censored outcomes, scalar covariates, and two-dimensional image predictors. The package provides documented fitting functions, coefficient reconstruction, prediction, reproducible simulations, tests, and GitHub Actions.

The model uses an image basis expansion, inverse censoring probability weights, penalized quantile covariance component extraction, and `quantreg::crq(..., method = "PengHuang")` for the final quantile regression. Version 0.1.2 defaults to IGACV tuning, Beran conditioning on all scalar covariates, and the original event indicators. The example below uses the main simulation mechanism and these settings.

## Installation

Install the runtime dependencies in R:

```r
install.packages(c("quantreg", "survival"))
```

From the directory containing the delivered source archive:

```r
install.packages("sfcqr_0.1.2.tar.gz", repos = NULL, type = "source")
library(sfcqr)
```

Alternatively, from the parent of the package source directory, run:

```sh
R CMD INSTALL sfcqr
```

After the repository is published, replace `USERNAME` with its owner:

```r
install.packages("remotes")
remotes::install_github("USERNAME/sfcqr", ref = "v0.1.2")
```

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
\qquad A_{i\ell} \sim \operatorname{Unif}(0,4),
$$

$$
C_1(s) = B_1(s) + B_8(s),
\qquad C_2(s) = 0.8B_{22}(s) + 1.2B_{29}(s).
$$

The main-branch generator uses

$$
\log T_i = \{\langle Z_i,C_1\rangle + X_{i1}\}\epsilon_i
            + \langle Z_i,C_2\rangle + 0.5X_{i2},
\qquad \epsilon_i \sim \operatorname{Unif}(-1,1),
$$

$$
\log C_i = 4X_{i1} + 4X_{i2} + U_i - 2,
\qquad U_i \sim \operatorname{Logistic}(0,4).
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

This example fits all three quantiles to one simulated dataset. The manuscript also considers `n = 1000` and `n = 2000` and summarizes 1,000 Monte Carlo replications. The example and software tests validate the current implementation; they do not establish that earlier manuscript tables were generated with this release.

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

`criterion` accepts `"IGACV"`, `"IAIC"`, or `"IBIC"`, with `"IGACV"` as the default. To use the earlier information-criterion choice, specify `criterion = "IBIC"` explicitly.

Censoring weights support `"beran"`, `"marginal"`, `"lognormal"`, `"loglogistic"`, and `"cox"`. For Beran weights, `beran_covariates = "all"` conditions on all columns of `X`, as specified in Supplementary Section S4. The default is `"all"`; `"first"` remains available as an explicit option. The default bandwidth is `n^(-1 / (p + 4))`, where `p = ncol(X)`; the example supplies `n^(-1/6)` explicitly because there are two covariates. With no scalar covariates, the package uses marginal weights. The original spelling `"bernan"` remains a compatibility alias.

`min_survival` bounds the estimated censoring survival probability below, and weights use the right-continuous survival probability at the observed time, `G(Y)`, including censoring events tied at that time. The default `force_last_event = FALSE` fits that distribution using the original event indicators. The numerator of each IPCW weight is always the original `delta`, so censored observations have zero supervision weight. Setting `force_last_event = TRUE` explicitly applies the legacy largest-time adjustment to the censoring-distribution fit only; it does not change the weight numerator or the final quantile regression event indicators.

## Component extraction and changes in version 0.1.2

For a censored quantile covariance vector `q` and a basis penalty `P`, extraction solves

$$
A_\lambda = I + \lambda P,\qquad
A_\lambda a = q,\qquad
w = \frac{a}{\sqrt{a^\top A_\lambda a}}.
$$

This gives `t(w) %*% A_lambda %*% w = 1` in the original residual-predictor coordinates. When `scale.X = TRUE`, the code uses an equivalent coordinate transformation, including the Euclidean part of the metric, then maps the scores and fitted coefficients back to the original basis. `plsfit_cqcov()` returns these directions in `directions_raw` and the metric in `metric`; the final projection matrix `W` also includes the deflation mappings.

The marginal event-time threshold is the first observed time at which the normalized IPCW empirical CDF reaches `tau`. It uses the original event indicators and does not interpolate between observed times or fit an intercept-only `crq()` model. The `grid` argument continues to control candidate quantile regressions and their integrated criteria, but is unused by this empirical threshold.

Version 0.1.2 corrects the extraction direction and threshold and changes the defaults listed above. Fitted values and selected models can therefore differ from version 0.1.1. Explicitly specifying earlier tuning or censoring options does not restore the earlier component algorithm. See [NEWS.md](NEWS.md) for the release changes.

## Results and prediction

The fitted object has class `sfcqr`:

| Interface or field | Result |
| --- | --- |
| `fit$M`, `fit$lambda` | Selected component count and smoothing penalty |
| `fit$tuning` | Candidate scores for the selected integrated criterion; unavailable candidates receive `Inf` |
| `fit$diagnostics` | Warnings and reasons for unavailable candidate fits |
| `coef(fit)` | Original-scale intercept, scalar coefficients, and basis coefficients |
| `coef(fit, type = "basis")` | Basis coefficient vector |
| `coef(fit, type = "image")` | Reconstructed image coefficient on the specified grid |
| `predict(fit, newdata, type = "response")` | Positive time prediction at the target quantile |
| `predict(fit, newdata, type = "link")` | Prediction on the log-time scale |

New data must supply `X` and `Z` with the training column order. There is no need to supply `Bxy` again. Use `drop = FALSE` when subsetting a single column or row to preserve matrix dimensions.

For the main simulator, `dat$truth$taus` contains the reference quantiles, `dat$truth$coef_tau` contains the corresponding 36-dimensional basis coefficients, and `dat$truth$coef_X_tau` contains the two scalar coefficients. The true log-time intercepts are in `dat$truth$intercept_tau`; `dat$setting` records the simulation configuration.

See `?sfcqr`, `?simulate_sfcqr_main`, and `?get_weights` for details. `gene_data(seed, n = 1000)` provides a compatibility entry point for the supplied main-branch generator. The separate `simulate_sfcqr_data()` function remains available for generic software examples; it is not the data generator used in the README's main-setting example.

## Changes from the research script

- Removed `rm(list = ls())`, global `library()` calls, and unused dependencies so that functions do not alter the user's workspace.
- Replaced hard-coded image dimensions, penalty candidates, and scalar coefficient counts with validated arguments and dynamic dimensions.
- Added checks for incompatible dimensions and handling for dropped matrix dimensions, constant columns, zero-variance components, singular matrices, and unavailable quantile fits.
- Unified centering, scaling, and component mappings between fitting and prediction, and converted fitted coefficients back to the original input scale. The reconstruction checks in the example verify this mapping.
- Implemented the original `imager::imgradient()` default `scheme = 3` rotation-invariant `3 x 3` difference stencil and replicated boundaries in base R, applying the stencil twice in each direction to construct the smoothing penalty. The implementation was checked numerically against the original gradient scheme, allowing removal of the `imager` dependency. Existing penalty matrices can also be supplied directly.
- Added the supplied main-branch data generator, IGACV selection, and Beran conditioning on all scalar covariates; version 0.1.2 uses these fitting options by default.
- Aligned the component optimizer and IPCW empirical threshold with the stated mathematical procedure and preserved the original event indicators in the supervision weights.

## Development and GitHub publication

Install `testthat` before running the checks. To run the test suite from the source directory:

```sh
Rscript -e 'testthat::test_local(".")'
```

From the parent of the source directory:

```sh
R CMD build sfcqr
R CMD check --no-manual sfcqr_0.1.2.tar.gz
```

The repository includes a GitHub Actions workflow that checks the package on Linux, macOS, and Windows. The author and maintainer is Dingke Tang (`dtang@uottawa.ca`), and the source is distributed under the MIT license.

Follow the [GitHub publication guide](GITHUB_SUBMISSION.md) for repository creation, authentication, pushing the source, and publishing a tagged release.

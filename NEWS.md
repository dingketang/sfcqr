# sfcqr 0.1.2

- Component extraction now solves `(I + lambda * P) a = q` and normalizes
  the direction in the same penalized metric.
- Predictor standardization is an equivalent coordinate transformation of
  the original basis metric; changing `scale.X` preserves the component space.
- The marginal event-time quantile is obtained by inverting the normalized
  IPCW empirical CDF, using positive event weights and no interpolation.
- Censored observations have zero supervision weight. `get_weights()` always
  uses the original event indicator in its numerator; the optional largest-time
  adjustment affects the censoring-distribution fit only.
- Fitting defaults are now `criterion = "IGACV"`, `beran_covariates = "all"`,
  and `force_last_event = FALSE`. The earlier criterion and censoring options
  remain available when explicitly selected.
- `plsfit_cqcov()` exposes original-scale extraction directions, their metric,
  and the effective IPCW weights for numerical verification.
- Added regression checks for the penalized optimizer, IPCW threshold,
  event-indicator preservation, scaling, deflation, and coefficient mapping.
- These corrections can change fitted values and selected models relative
  to version 0.1.1; the previous extraction algorithm is not retained.

# sfcqr 0.1.1

- English README and GitHub submission guide.
- A reproducible generator following the supplied main simulation scenario:
  100 x 100 images, 36 tensor B-splines, uniform heteroscedastic errors, and
  covariate-dependent logistic censoring with shift -2.
- A README example using n = 500, all three manuscript quantiles, the full
  tuning grids, multivariate Beran weights, and IGACV selection.
- Optional multivariate Beran conditioning and selectable integrated criteria;
  the original fitting defaults remain available.

# sfcqr 0.1.0

- Initial package release of supervised functional censored quantile regression.
- Configurable image dimensions, component counts, penalty candidates, and censoring estimation.
- Raw-scale coefficient, image reconstruction, prediction, and print methods.
- Reproducible illustrative simulation and README smoke example.
- Input validation, numerical safeguards, and base R image roughness penalties.
- GitHub submission instructions and cross-platform R package checks.

unsurv 0.6.0 (2026-07-22)
-------------------------
- Added `unsurv_compare()` to summarize and compare multiple cluster partitions
  (e.g., an `unsurv` curve-based partition against scalar-risk or covariate-PCA
  baselines) against observed survival outcomes, reporting cluster-size balance,
  Adjusted Rand Index agreement against a reference partition, and per-cluster
  Kaplan-Meier medians.
- Added `autoplot()`/`plot()` methods for `"unsurv_compare"` objects, producing
  Kaplan-Meier curves faceted by comparison method.

unsurv 0.5.0 (2026-03-12)
-------------------------
- Added introductory vignette covering fitting, visualization, prediction, and stability checks.
- Expanded test suite for weighting options, auto-K selection, monotonic enforcement, and prediction validation.
- Added citation metadata and NEWS file for release tracking.
- Set up GitHub Actions workflow for R CMD check.

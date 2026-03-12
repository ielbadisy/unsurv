## Test environments
- Local: Ubuntu 24.04.3 LTS, R 4.5.1 (2025-06-13)

## R CMD check results
- 0 errors | 0 warnings | 1 note

### Notes
- "unable to verify current time" arises from the sandbox clock; it does not affect package behavior.

## Additional validation
- Vignettes build successfully (`rmarkdown::html_vignette`).
- Test suite (`testthat`) runs clean during `R CMD check --as-cran`.

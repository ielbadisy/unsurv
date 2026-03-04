## Test environments
- Local: Ubuntu 24.04.3 LTS, R 4.5.1 (2025-06-13)

## R CMD check results
- 0 errors | 0 warnings | 2 notes

### Notes
1) CRAN incoming reports:
   - "Package has a VignetteBuilder field but no prebuilt vignette index."  
     The `R CMD build` tarball includes `inst/doc/index.html` plus the prebuilt vignette outputs (`unsurv-intro.html`, `unsurv-intro.Rmd`, `unsurv-intro.R`). Please let me know if an additional index format is required.
   - URL checks fail while offline (github.com, ielbadisy.github.io) because the check environment lacks internet access.
2) "unable to verify current time" arises from the sandbox clock; it does not affect package behavior.

## Additional validation
- Vignettes build successfully (`rmarkdown::html_vignette`).
- Test suite (`testthat`) runs clean during `R CMD check --as-cran`.

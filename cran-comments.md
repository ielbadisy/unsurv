## Resubmission

This is a resubmission (version 0.7.2). Changes since the previous submission:

- The `Description` field now points to the published method reference,
  El Badisy (2026) <doi:10.1093/bioadv/vbag218>, and `inst/CITATION` was
  switched to the corresponding article entry.
- New exported function `unsurv_compare()` (with `print()`/`autoplot()`/`plot()`
  methods) for comparing cluster partitions against observed survival outcomes.
- The introductory vignette gained the paper's worked example, run from a small
  pre-computed data set in `inst/extdata/` (no additional packages needed at
  build time).
- Dropped the `dplyr` and `tidyr` dependencies: the two plot helpers'
  internal reshape step now uses base R only.

## Test environments

- Local: Ubuntu 24.04.3 LTS, R 4.5.1 (2025-06-13)

## R CMD check results

`R CMD check --as-cran` on the source tarball: 0 errors | 0 warnings | 1 note.

### Note

- "New maintainer" — the Maintainer field's family-name casing changed from
  "EL BADISY" to "El Badisy" (title case). Same person, same email address,
  no change in maintainership.

## Additional validation

- Vignette builds and re-builds cleanly (`rmarkdown::html_vignette`).
- Test suite (`testthat`, edition 3) runs clean under `R CMD check --as-cran`.
- The DOI in `DESCRIPTION` resolves.

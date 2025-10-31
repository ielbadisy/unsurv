predict_unsurv <- function(fit, S_new, clamp = TRUE) {
  stopifnot(inherits(fit, "unsurv"))
  S_new <- as.matrix(S_new)
  if (ncol(S_new) != length(fit$times))
    stop("S_new must have the same number of columns/time grid as the fit.")
  if (clamp) S_new <- .clamp01(S_new)
  if (isTRUE(fit$enforce_monotone)) S_new <- .enforce_monotone(S_new)
  if (is.numeric(fit$smooth_median_width) &&
      fit$smooth_median_width >= 3 &&
      fit$smooth_median_width %% 2 == 1) {
    S_new <- .smooth_median(S_new, fit$smooth_median_width)
  }

  # recreate weighted (and optionally standardized) feature space
  X_new <- .weight_features(S_new, fit$weights, fit$distance)
  X_med <- .weight_features(fit$medoids, fit$weights, fit$distance)
  X_new <- as.matrix(X_new); X_med <- as.matrix(X_med)

  if (isTRUE(fit$standardize_cols)) {
    mu <- fit$center; sdv <- fit$scale
    if (is.null(mu) || is.null(sdv)) stop("fit was standardized but center/scale not stored.")
    X_new <- sweep(X_new, 2, mu, `-`); X_new <- sweep(X_new, 2, sdv, `/`)
    X_med <- sweep(X_med, 2, mu, `-`); X_med <- sweep(X_med, 2, sdv, `/`)
  }

  if (fit$distance == "L2") {
    rn <- rowSums(X_new^2)
    rm <- rowSums(X_med^2)
    d2 <- outer(rn, rm, `+`) - 2 * (X_new %*% t(X_med))
    d2[d2 < 0] <- 0
    d <- sqrt(d2)
  } else {
    d <- matrix(NA_real_, nrow(X_new), nrow(X_med))
    for (k in seq_len(nrow(X_med))) {
      d[, k] <- rowSums(abs(X_new - matrix(X_med[k, ], nrow(X_new), byrow = TRUE)))
    }
  }
  max.col(-d)  # K cluster id per row
}

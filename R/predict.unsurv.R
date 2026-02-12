
## ---- prediction: S3 method (export as predict.unsurv in package) ----
predict.unsurv <- function(object, newdata, clamp = TRUE, ...) {
  stopifnot(inherits(object, "unsurv"))
  S_new <- as.matrix(newdata)
  if (ncol(S_new) != length(object$times))
    stop("newdata must have the same number of columns/time grid as the fit.")
  if (isTRUE(clamp)) S_new <- .clamp01(S_new)
  if (isTRUE(object$enforce_monotone)) S_new <- .enforce_monotone(S_new)
  if (is.numeric(object$smooth_median_width) &&
      object$smooth_median_width >= 3 &&
      object$smooth_median_width %% 2 == 1) {
    S_new <- .smooth_median(S_new, object$smooth_median_width)
  }
  
  # recreate weighted (and optionally standardized) feature space
  X_new <- .weight_features(S_new, object$weights, object$distance)
  X_med <- .weight_features(object$medoids, object$weights, object$distance)
  X_new <- as.matrix(X_new); X_med <- as.matrix(X_med)
  
  if (isTRUE(object$standardize_cols)) {
    mu <- object$center; sdv <- object$scale
    if (is.null(mu) || is.null(sdv)) stop("fit was standardized but center/scale not stored.")
    X_new <- sweep(X_new, 2, mu, `-`); X_new <- sweep(X_new, 2, sdv, `/`)
    X_med <- sweep(X_med, 2, mu, `-`); X_med <- sweep(X_med, 2, sdv, `/`)
  }
  
  if (object$distance == "L2") {
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



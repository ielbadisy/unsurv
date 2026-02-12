
## ---- printing & summary: S3 methods (export in package) ----
print.unsurv <- function(x, ...) {
  cat("unsurv (PAM) fit\n",
      "  K:", x$K, "\n",
      "  distance:", x$distance, " silhouette_mean:", sprintf("%.3f", x$silhouette_mean), "\n",
      "  n:", length(x$clusters), " Q:", length(x$times), "\n", sep = "")
  invisible(x)
}

summary.unsurv <- function(object, ...) {
  out <- list(
    K = object$K,
    silhouette_mean = object$silhouette_mean,
    size = as.integer(table(object$clusters))
  )
  class(out) <- "summary.unsurv"
  out
}

print.summary.unsurv <- function(x, ...) {
  cat("K:\n"); print(x$K)
  cat("\nSilhouette mean:\n"); print(x$silhouette_mean)
  cat("\nCluster sizes:\n"); print(x$size)
  invisible(x)
}


## ---- computing stability (export in package as unsurv_stability) ----
unsurv_stability <- function(
    S, times, fit,
    B = 30, frac = 0.5, mode = c("bootstrap","subsample"),
    jitter_sd = 0.01,
    weight_perturb = 0.30,
    eps_jitter = 0.02,
    return_distribution = TRUE
){
  mode <- match.arg(mode)
  stopifnot(inherits(fit, "unsurv"))
  if (ncol(S) != length(times)) stop("S/times dimension mismatch.")
  
  .cluster_stability(
    S, times,
    B = B, frac = frac, mode = mode,
    jitter_sd = jitter_sd, weight_perturb = weight_perturb,
    K = fit$K,
    distance = fit$distance,
    standardize_cols = isTRUE(fit$standardize_cols),
    enforce_monotone = isTRUE(fit$enforce_monotone),
    smooth_median_width = ifelse(is.null(fit$smooth_median_width), 0, fit$smooth_median_width),
    eps_jitter = eps_jitter,
    seed = NULL,
    return_distribution = return_distribution
  )
}






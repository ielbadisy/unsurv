#' Print an unsurv model
#'
#' @param x An object of class \code{"unsurv"}.
#' @param ... Unused.
#'
#' @return Invisibly returns \code{x}.
#' @export
print.unsurv <- function(x, ...) {
  cat("unsurv (PAM) fit\n",
      "  K:", x$K, "\n",
      "  distance:", x$distance,
      " silhouette_mean:", sprintf("%.3f", x$silhouette_mean), "\n",
      "  n:", length(x$clusters),
      " Q:", length(x$times), "\n", sep = "")
  invisible(x)
}

#' Summarize an unsurv model
#'
#' @param object An object of class \code{"unsurv"}.
#' @param ... Unused.
#'
#' @return An object of class \code{"summary.unsurv"} with elements:
#'   \itemize{
#'     \item \code{K}: number of clusters
#'     \item \code{silhouette_mean}: mean silhouette width
#'     \item \code{size}: cluster sizes
#'   }
#' @export
summary.unsurv <- function(object, ...) {
  out <- list(
    K = object$K,
    silhouette_mean = object$silhouette_mean,
    size = as.integer(table(object$clusters))
  )
  class(out) <- "summary.unsurv"
  out
}

#' Print a summary of an unsurv model
#'
#' @param x An object of class \code{"summary.unsurv"}.
#' @param ... Unused.
#'
#' @return Invisibly returns \code{x}.
#' @export
print.summary.unsurv <- function(x, ...) {
  cat("K:\n"); print(x$K)
  cat("\nSilhouette mean:\n"); print(x$silhouette_mean)
  cat("\nCluster sizes:\n"); print(x$size)
  invisible(x)
}

#' Stability assessment for an unsurv clustering
#'
#' Computes a resampling-based stability score for a fitted \code{unsurv} model
#' using the Adjusted Rand Index (ARI) computed on overlap sets across resamples.
#'
#' @param S Numeric matrix of survival probabilities used for stability assessment
#'   (\eqn{n \times m}), with columns matching \code{times}.
#' @param times Numeric vector of time grid points (length \eqn{m}).
#' @param fit An object of class \code{"unsurv"}, typically returned by \code{\link{unsurv}}.
#' @param B Integer; number of resamples.
#' @param frac Numeric in (0, 1]; fraction of rows sampled per resample.
#' @param mode Resampling mode: \code{"bootstrap"} (with replacement) or
#'   \code{"subsample"} (without replacement).
#' @param jitter_sd Nonnegative numeric; curve-space noise level applied before clamping/monotone enforcement.
#' @param weight_perturb Numeric in \eqn{[0,1]}; blends trapezoidal weights with a random simplex
#'   to perturb the weighting scheme.
#' @param eps_jitter Nonnegative numeric; feature-space jitter used inside the clustering during resamples.
#' @param return_distribution Logical; if \code{TRUE}, returns the full ARI distribution, else returns only the mean.
#'
#' @return If \code{return_distribution = TRUE}, a list with:
#'   \itemize{
#'     \item \code{mean}: mean ARI across resample-pair overlaps
#'     \item \code{aris}: numeric vector of ARIs
#'   }
#' Otherwise, returns a single numeric mean ARI.
#'
#' @examples
#' if (requireNamespace("cluster", quietly = TRUE)) {
#'   set.seed(2025)
#'   n <- 120; Q <- 60
#'   times <- seq(0, 5, length.out = Q)
#'   rates <- c(0.12, 0.38, 0.8)
#'   grp <- sample(1:3, n, TRUE, c(0.4, 0.4, 0.2))
#'   S <- t(vapply(1:n, function(i)
#'     pmin(pmax(exp(-rates[grp[i]] * times) + rnorm(Q, 0, 0.01), 0), 1),
#'     numeric(Q)
#'   ))
#'
#'   fit <- unsurv(S, times, K = NULL, K_max = 6, distance = "L2",
#'                enforce_monotone = TRUE, standardize_cols = FALSE,
#'                eps_jitter = 0, seed = NULL)
#'
#'   stab <- unsurv_stability(S, times, fit, B = 30, frac = 0.55, mode = "bootstrap",
#'                            jitter_sd = 0.3, weight_perturb = 0.0, eps_jitter = 0.3,
#'                            return_distribution = TRUE)
#'   stab$mean
#' }
#' @export
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
  if (ncol(S) != length(times)) stop("S/times dimension mismatch.", call. = FALSE)

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

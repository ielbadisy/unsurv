#' Print an unsurv model
#'
#' @param x An object of class \code{"unsurv"}.
#' @param ... Unused.
#'
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' if (requireNamespace("cluster", quietly = TRUE)) {
#'   set.seed(1)
#'   times <- seq(0, 4, length.out = 20)
#'   grp <- rep(1:2, each = 6)
#'   rates <- c(0.2, 0.6)
#'   S <- sapply(times, function(t) exp(-rates[grp] * t))
#'   fit <- unsurv(S, times, K = 2)
#'   print(fit)
#' }
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
#'
#' @examples
#' if (requireNamespace("cluster", quietly = TRUE)) {
#'   set.seed(1)
#'   times <- seq(0, 4, length.out = 20)
#'   grp <- rep(1:2, each = 6)
#'   rates <- c(0.2, 0.6)
#'   S <- sapply(times, function(t) exp(-rates[grp] * t))
#'   fit <- unsurv(S, times, K = 2)
#'   summary(fit)
#' }
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
#'
#' @examples
#' if (requireNamespace("cluster", quietly = TRUE)) {
#'   set.seed(1)
#'   times <- seq(0, 4, length.out = 20)
#'   grp <- rep(1:2, each = 6)
#'   rates <- c(0.2, 0.6)
#'   S <- sapply(times, function(t) exp(-rates[grp] * t))
#'   fit <- unsurv(S, times, K = 2)
#'   s <- summary(fit)
#'   print(s)
#' }
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
#' @param frac Numeric in (0, 1]; size of each resample as a fraction of \code{nrow(S)}.
#' @param mode Resampling mode: \code{"bootstrap"} (with replacement, default) or
#'   \code{"subsample"} (without replacement).
#' @param jitter_sd Nonnegative numeric; sd of Gaussian noise added to the curves before
#'   clamping/monotone enforcement. Default 0 (no curve noise).
#' @param weight_perturb Numeric in \eqn{[0,1]}; blends the time weights of \code{fit} with a
#'   random simplex to perturb the weighting scheme. Default 0 (weights of \code{fit}).
#' @param eps_jitter Nonnegative numeric; feature-space tie-breaking jitter used inside each
#'   refit. Defaults to the value used by \code{fit}.
#' @param return_distribution Logical; if \code{TRUE}, returns the full ARI distribution, else returns only the mean.
#' @param seed Optional integer seed for the resampling; \code{NULL} (default) uses the current
#'   random number stream.
#'
#' @details Each resample is refitted with the settings of \code{fit}, including its time weights
#'   \code{fit$weights} and its tie-breaking jitter \code{fit$eps_jitter}, and the partitions of
#'   two refits are compared by the ARI on the patients they share (Ben-Hur et al., 2002).
#'
#'   By default only the patients are resampled: a nonparametric bootstrap of size \eqn{n}, the
#'   scheme recommended by Hennig (2007). The index then measures how reproducible the partition is
#'   under sampling of patients. Noise added to the curves (\code{jitter_sd}), to the weights
#'   (\code{weight_perturb}) or to the features (\code{eps_jitter} above the value of the fit) is
#'   available as a sensitivity analysis but is off by default: its scale is arbitrary, feature
#'   noise adds \eqn{2 m \sigma^2} to every squared distance and so grows with the number of grid
#'   points \eqn{m}, and in simulations it lowered stability most for the true number of clusters
#'   and pulled the most stable \eqn{K} towards 2. Before version 0.7.5 the defaults were
#'   \code{frac = 0.5}, \code{jitter_sd = 0.01}, \code{weight_perturb = 0.3} and
#'   \code{eps_jitter = 0.02}.
#'
#' @references
#' Ben-Hur, A., Elisseeff, A. and Guyon, I. (2002). A stability based method for discovering
#' structure in clustered data. \emph{Pacific Symposium on Biocomputing}, 7, 6--17.
#'
#' Hennig, C. (2007). Cluster-wise assessment of cluster stability. \emph{Computational
#' Statistics & Data Analysis}, 52(1), 258--271. \doi{10.1016/j.csda.2006.11.025}
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
#'   n <- 60; Q <- 40
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
#'   stab <- unsurv_stability(S, times, fit, B = 8, frac = 0.55, mode = "bootstrap",
#'                            jitter_sd = 0.3, weight_perturb = 0.0, eps_jitter = 0.3,
#'                            return_distribution = TRUE)
#'   stab$mean
#' }
#' @export
unsurv_stability <- function(
    S, times, fit,
    B = 30, frac = 1, mode = c("bootstrap","subsample"),
    jitter_sd = 0,
    weight_perturb = 0,
    eps_jitter = fit$eps_jitter,
    return_distribution = TRUE,
    seed = NULL
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
    weights = fit$weights,
    seed = seed,
    return_distribution = return_distribution
  )
}

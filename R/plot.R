#' Unsupervised clustering of individualized survival curves
#'
#' Clusters individuals using their survival-probability curves evaluated on a
#' common time grid. The method computes a weighted feature representation of
#' the curves and applies PAM (Partitioning Around Medoids) on the resulting
#' dissimilarity matrix. If \code{K} is not provided, it is selected by maximizing
#' the mean silhouette width over \code{K = 2, ..., K_max}.
#'
#' @param S Numeric matrix of survival probabilities with shape \eqn{n \times m}.
#'   Rows are individuals; columns correspond to the time grid \code{times}.
#'   Values are clamped to \eqn{[0,1]} internally.
#' @param times Numeric vector of length \eqn{m} giving the time grid (must be
#'   strictly increasing).
#' @param K Optional integer number of clusters. If \code{NULL}, selected via
#'   silhouette width.
#' @param K_max Maximum number of clusters considered when \code{K} is \code{NULL}.
#'   Effective maximum is \code{min(K_max, nrow(S))}.
#' @param distance Distance type: \code{"L2"} (euclidean) or \code{"L1"} (manhattan).
#' @param weights Optional nonnegative numeric vector of length \eqn{m} giving
#'   time-point weights. If \code{NULL}, trapezoidal weights are used.
#' @param enforce_monotone Logical; if \code{TRUE}, enforces non-increasing curves
#'   over time.
#' @param smooth_median_width Integer; if \eqn{\ge 3} and odd, applies median
#'   smoothing along the time axis.
#' @param standardize_cols Logical; if \code{TRUE}, standardizes feature columns
#'   prior to clustering (center/scale stored for prediction).
#' @param eps_jitter Nonnegative numeric; standard deviation of small Gaussian
#'   noise added in feature space to break ties.
#' @param seed Optional integer seed for reproducibility.
#'
#' @details
#' This function requires the \pkg{cluster} package for PAM clustering and
#' silhouette widths.
#'
#' The returned object stores medoid curves and metadata required for prediction
#' on new curves via \code{\link[stats]{predict}} (method \code{predict.unsurv}).
#'
#' @return An object of class \code{"unsurv"} with components including:
#' \itemize{
#'   \item \code{clusters}: integer vector of cluster assignments
#'   \item \code{K}: number of clusters
#'   \item \code{times}: time grid
#'   \item \code{medoids}: medoid survival curves (one per cluster)
#'   \item \code{silhouette_mean}: mean silhouette width
#'   \item plus preprocessing/settings fields used for prediction
#' }
#'
#' @examples
#' if (requireNamespace("cluster", quietly = TRUE)) {
#'   set.seed(2025)
#'   n <- 40; Q <- 30
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
#'   print(fit)
#'   summary(fit)
#'   plot(fit)
#'
#'   pred <- predict(fit, S[1:5, ])
#'   pred
#' }
#' @export
unsurv <- function(
    S, times,
    K = NULL, K_max = 10,
    distance = c("L2", "L1"),
    weights = NULL,
    enforce_monotone = TRUE,
    smooth_median_width = 0,
    standardize_cols = FALSE,
    eps_jitter = 0.001,
    seed = NULL
){
  if (!requireNamespace("cluster", quietly = TRUE)) {
    stop("Install 'cluster' package.", call. = FALSE)
  }

  distance <- match.arg(distance)
  if (!is.null(seed)) set.seed(seed)

  S <- .check_inputs(S, times)
  S <- .clamp01(S)
  if (isTRUE(enforce_monotone)) S <- .enforce_monotone(S)
  if (smooth_median_width >= 3 && smooth_median_width %% 2 == 1) {
    S <- .smooth_median(S, smooth_median_width)
  }

  # weights (default trapezoid)
  w <- if (is.null(weights)) .trap_weights(times) else {
    if (length(weights) != ncol(S)) stop("'weights' length must equal number of time points.", call. = FALSE)
    if (any(weights < 0)) stop("'weights' must be nonnegative.", call. = FALSE)
    sw <- sum(weights)
    if (sw <= 0) stop("'weights' must sum to a positive value.", call. = FALSE)
    weights / sw
  }

  # weighted features
  X <- .weight_features(S, w, distance)

  # optional column standardization
  center <- NULL; scalev <- NULL
  if (isTRUE(standardize_cols)) {
    X <- scale(X)
    center <- attr(X, "scaled:center")
    scalev <- attr(X, "scaled:scale"); scalev[scalev == 0] <- 1
    X <- as.matrix(X)
  } else {
    X <- as.matrix(X)
  }

  # tiny jitter to break ties / determinism
  if (eps_jitter > 0) {
    X <- X + matrix(stats::rnorm(length(X), 0, eps_jitter), nrow(X))
  }

  D <- .build_dist(X, distance)
  n <- nrow(S)

  # choose K by silhouette if not provided
  if (is.null(K)) {
    K_max <- min(as.integer(K_max), n)
    if (K_max < 2) stop("'K_max' must be >= 2 (and <= nrow(S)) when K is NULL.", call. = FALSE)
    rng <- 2:K_max
    sil <- rep(NA_real_, length(rng))
    for (i in seq_along(rng)) {
      pm <- cluster::pam(D, k = rng[i], diss = TRUE)
      sil[i] <- mean(cluster::silhouette(pm$clustering, D)[, 3])
    }
    K <- rng[which.max(sil)]
  } else {
    K <- as.integer(K)
    if (K < 2) stop("'K' must be >= 2.", call. = FALSE)
    if (K > n) stop("'K' cannot exceed nrow(S).", call. = FALSE)
  }

  pam_fit <- cluster::pam(D, k = K, diss = TRUE)
  clusters <- pam_fit$clustering
  medoid_indices <- pam_fit$id.med
  medoids <- S[medoid_indices, , drop = FALSE]
  sil_mean <- mean(cluster::silhouette(clusters, D)[, 3])

  out <- list(
    clusters = as.integer(clusters),
    K = as.integer(K),
    times = as.numeric(times),
    medoid_indices = as.integer(medoid_indices),
    medoids = medoids,
    silhouette_mean = as.numeric(sil_mean),
    weights = w,
    distance = distance,
    standardize_cols = isTRUE(standardize_cols),
    center = center,
    scale = scalev,
    enforce_monotone = isTRUE(enforce_monotone),
    smooth_median_width = smooth_median_width,
    eps_jitter = eps_jitter,
    seed = seed
  )
  class(out) <- "unsurv"
  out
}

#' Plot medoid survival curves from an unsurv fit
#'
#' Produces a base R plot of the cluster medoid survival curves stored in the
#' fitted object.
#'
#' @param x An object of class \code{"unsurv"}.
#' @param ... Additional arguments passed to \code{\link[graphics]{matplot}}
#'   (e.g., \code{lwd}).
#'
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' if (requireNamespace("cluster", quietly = TRUE)) {
#'   set.seed(1)
#'   times <- seq(0, 4, length.out = 30)
#'   grp <- rep(1:2, each = 8)
#'   rates <- c(0.2, 0.55)
#'   S <- sapply(times, function(t) exp(-rates[grp] * t))
#'   fit <- unsurv(S, times, K = 2)
#'   plot(fit)
#' }
#' @export
plot.unsurv <- function(x, ...) {
  stopifnot(inherits(x, "unsurv"))
  times <- x$times
  med <- x$medoids
  if (is.null(med) || nrow(med) < 1L) stop("No medoids found in 'x$medoids'.", call. = FALSE)
  if (length(times) != ncol(med)) stop("Dimension mismatch: length(x$times) must equal ncol(x$medoids).", call. = FALSE)

  K <- nrow(med)
  graphics::matplot(times, t(med), type = "l", lty = 1,
                    xlab = "Time", ylab = "Survival", ...)
  graphics::legend("topright", legend = paste0("Cluster ", seq_len(K)),
                   lty = 1, bty = "n")
  invisible(x)
}

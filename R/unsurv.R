#' Unsupervised clustering of survival curves
#'
#' Fits an unsupervised clustering model on survival-probability curves evaluated
#' on a common time grid. Clustering is performed using PAM (Partitioning Around
#' Medoids) on a weighted feature representation of the curves.
#'
#' If \code{K} is \code{NULL}, the number of clusters is selected by maximizing
#' the mean silhouette width over \code{K = 2, ..., K_max}.
#'
#' @param S Numeric matrix of survival probabilities with shape \eqn{n \times m}.
#'   Rows are subjects, columns correspond to \code{times}. Values are clamped to \eqn{[0,1]}.
#' @param times Numeric vector of length \eqn{m} (strictly increasing time grid).
#' @param K Optional integer number of clusters. If \code{NULL}, selected by silhouette.
#' @param K_max Maximum \code{K} considered when \code{K} is \code{NULL}.
#' @param distance Distance type: \code{"L2"} (euclidean) or \code{"L1"} (manhattan).
#' @param weights Optional nonnegative vector of length \eqn{m} for time-point weights.
#'   If \code{NULL}, trapezoidal weights are used.
#' @param enforce_monotone Logical; enforce non-increasing survival curves over time.
#' @param smooth_median_width Integer; if \eqn{\ge 3} and odd, apply median smoothing along time.
#' @param standardize_cols Logical; standardize feature columns before clustering.
#' @param eps_jitter Nonnegative numeric; feature-space Gaussian jitter sd to break ties.
#' @param seed Optional integer seed.
#'
#' @details
#' Requires the \pkg{cluster} package (recommended in \code{Suggests}).
#'
#' @return An object of class \code{"unsurv"}.
#'
#' @examples
#' if (requireNamespace("cluster", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 80
#'   times <- seq(0, 5, length.out = 60)
#'   grp <- sample(1:2, n, TRUE)
#'   rates <- ifelse(grp == 1, 0.2, 0.6)
#'   S <- sapply(times, function(t) exp(-rates * t))
#'   S <- S + matrix(stats::rnorm(n * length(times), 0, 0.02), nrow = n)
#'   fit <- unsurv(S, times, K = NULL, K_max = 6, seed = 123)
#'   table(fit$clusters, grp)
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
) {
  .fit_unsurv(
    S = S, times = times,
    K = K, K_max = K_max,
    distance = distance,
    weights = weights,
    enforce_monotone = enforce_monotone,
    smooth_median_width = smooth_median_width,
    standardize_cols = standardize_cols,
    eps_jitter = eps_jitter,
    seed = seed
  )
}

#' @keywords internal
.fit_unsurv <- function(
  S, times,
  K = NULL, K_max = 10,
  distance = c("L2", "L1"),
  weights = NULL,
  enforce_monotone = TRUE,
  smooth_median_width = 0,
  standardize_cols = FALSE,
  eps_jitter = 0.001,
  seed = NULL
) {
  if (!requireNamespace("cluster", quietly = TRUE)) {
    stop("Package 'cluster' is required. Please install it.", call. = FALSE)
  }

  distance <- match.arg(distance)
  if (!is.null(seed)) set.seed(seed)

  S <- .check_inputs(S, times)
  S <- .clamp01(S)

  if (isTRUE(enforce_monotone)) S <- .enforce_monotone(S)
  if (is.numeric(smooth_median_width) && smooth_median_width >= 3 && smooth_median_width %% 2 == 1) {
    S <- .smooth_median(S, smooth_median_width)
  }

  # weights (default trapezoid)
  w <- if (is.null(weights)) .trap_weights(times) else {
    weights <- as.numeric(weights)
    if (length(weights) != ncol(S)) stop("'weights' length must equal number of time points.", call. = FALSE)
    if (anyNA(weights)) stop("'weights' contains NA.", call. = FALSE)
    if (any(weights < 0)) stop("'weights' must be nonnegative.", call. = FALSE)
    sw <- sum(weights)
    if (sw <= 0) stop("'weights' must sum to a positive value.", call. = FALSE)
    weights / sw
  }

  # weighted features
  X <- .weight_features(S, w, distance)

  # optional column standardization
  center <- NULL; scalev <- NULL
  X <- as.matrix(X)
  if (isTRUE(standardize_cols)) {
    Xs <- scale(X)
    center <- attr(Xs, "scaled:center")
    scalev <- attr(Xs, "scaled:scale"); scalev[scalev == 0] <- 1
    X <- as.matrix(Xs)
  }

  # tiny jitter to break ties
  eps_jitter <- as.numeric(eps_jitter)
  if (!is.finite(eps_jitter) || eps_jitter < 0) stop("'eps_jitter' must be >= 0.", call. = FALSE)
  if (eps_jitter > 0) {
    X <- X + matrix(stats::rnorm(length(X), 0, eps_jitter), nrow(X))
  }

  D <- .build_dist(X, distance)

  n <- nrow(S)
  if (!is.null(K)) {
    K <- as.integer(K)
    if (K < 2) stop("'K' must be >= 2.", call. = FALSE)
    if (K > n) stop("'K' cannot exceed nrow(S).", call. = FALSE)
  } else {
    K_max <- as.integer(K_max)
    if (K_max < 2) stop("'K_max' must be >= 2 when 'K' is NULL.", call. = FALSE)
    K_max <- min(K_max, n)
    rng <- 2:K_max
    sil <- rep(NA_real_, length(rng))
    for (i in seq_along(rng)) {
      pm <- cluster::pam(D, k = rng[i], diss = TRUE)
      sil[i] <- mean(cluster::silhouette(pm$clustering, D)[, 3])
    }
    K <- rng[which.max(sil)]
  }

  pam_fit <- cluster::pam(D, k = K, diss = TRUE)
  clusters <- pam_fit$clustering
  medoid_indices <- pam_fit$id.med
  medoids <- S[medoid_indices, , drop = FALSE]
  sil_mean <- mean(cluster::silhouette(clusters, D)[, 3])

  structure(
    list(
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
      seed = seed,
      call = match.call()
    ),
    class = "unsurv"
  )
}

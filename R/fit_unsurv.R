fit_unsurv <- function(
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
  if (!requireNamespace("cluster", quietly = TRUE)) stop("Install 'cluster' package.")
  distance <- match.arg(distance)
  if (!is.null(seed)) set.seed(seed)

  S <- .check_inputs(S, times)
  S <- .clamp01(S)
  if (enforce_monotone) S <- .enforce_monotone(S)
  if (smooth_median_width >= 3 && smooth_median_width %% 2 == 1)
    S <- .smooth_median(S, smooth_median_width)

  # weights (default trapezoid)
  w <- if (is.null(weights)) .trap_weights(times) else {
    if (length(weights) != ncol(S)) stop("'weights' length must equal number of time points.")
    if (any(weights < 0)) stop("'weights' must be nonnegative.")
    weights / sum(weights)
  }

  # weighted features
  X <- .weight_features(S, w, distance)

  # optional column standardization
  center <- NULL; scalev <- NULL
  if (standardize_cols) {
    X <- scale(X)
    center <- attr(X, "scaled:center")
    scalev <- attr(X, "scaled:scale"); scalev[scalev == 0] <- 1
    X <- as.matrix(X)
  } else {
    X <- as.matrix(X)
  }

  # tiny jitter to break ties / determinism
  if (eps_jitter > 0) {
    X <- X + matrix(rnorm(length(X), 0, eps_jitter), nrow(X))
  }

  D <- .build_dist(X, distance)

  # choose K by silhouette if not provided
  if (is.null(K)) {
    rng <- 2:max(2, K_max)
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

  out <- list(
    clusters = clusters,
    K = K,
    times = times,
    medoid_indices = medoid_indices,
    medoids = medoids,
    silhouette_mean = sil_mean,
    weights = w,
    distance = distance,
    standardize_cols = standardize_cols,
    center = center,
    scale = scalev,
    enforce_monotone = enforce_monotone,
    smooth_median_width = smooth_median_width,
    eps_jitter = eps_jitter,
    seed = seed
  )
  class(out) <- "unsurv"
  out
}


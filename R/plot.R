## ---- main user function (export in package as 'unsurv') ----
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
  if (!requireNamespace("cluster", quietly = TRUE)) stop("Install 'cluster' package.")
  distance <- match.arg(distance)
  if (!is.null(seed)) set.seed(seed)
  
  S <- .check_inputs(S, times)
  S <- .clamp01(S)
  if (isTRUE(enforce_monotone)) S <- .enforce_monotone(S)
  if (smooth_median_width >= 3 && smooth_median_width %% 2 == 1)
    S <- .smooth_median(S, smooth_median_width)
  
  # weights (default trapezoid)
  w <- if (is.null(weights)) .trap_weights(times) else {
    if (length(weights) != ncol(S)) stop("'weights' length must equal number of time points.")
    if (any(weights < 0)) stop("'weights' must be nonnegative.")
    sw <- sum(weights)
    if (sw <= 0) stop("'weights' must sum to a positive value.")
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
    if (K_max < 2) stop("'K_max' must be >= 2 (and <= nrow(S)) when K is NULL.")
    rng <- 2:K_max
    sil <- rep(NA_real_, length(rng))
    for (i in seq_along(rng)) {
      pm <- cluster::pam(D, k = rng[i], diss = TRUE)
      sil[i] <- mean(cluster::silhouette(pm$clustering, D)[, 3])
    }
    K <- rng[which.max(sil)]
  } else {
    K <- as.integer(K)
    if (K < 2) stop("'K' must be >= 2.")
    if (K > n) stop("'K' cannot exceed nrow(S).")
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


## ---- base plotting: S3 method (export as plot.unsurv in package) ----
plot.unsurv <- function(x, ...) {
  stopifnot(inherits(x, "unsurv"))
  times <- x$times
  med <- x$medoids
  if (is.null(med) || nrow(med) < 1L) stop("No medoids found in 'x$medoids'.")
  if (length(times) != ncol(med)) stop("Dimension mismatch: length(x$times) must equal ncol(x$medoids).")
  
  K <- nrow(med)
  graphics::matplot(times, t(med), type = "l", lty = 1,
                    xlab = "Time", ylab = "Survival", ...)
  graphics::legend("topright", legend = paste0("Cluster ", seq_len(K)),
                   lty = 1, bty = "n")
  invisible(x)
}





#' @keywords internal
.clamp01 <- function(M) {
  M[M < 0] <- 0
  M[M > 1] <- 1
  M
}

#' Validate inputs for unsurv
#' @keywords internal
.check_inputs <- function(S, times) {
  S <- as.matrix(S)
  if (!is.numeric(times) || any(!is.finite(times))) stop("'times' must be finite numeric.", call. = FALSE)
  if (length(times) < 2) stop("'times' must have length >= 2.", call. = FALSE)
  if (is.unsorted(times, strictly = TRUE)) stop("'times' must be strictly increasing.", call. = FALSE)
  if (ncol(S) != length(times)) stop("ncol(S) must equal length(times).", call. = FALSE)
  if (nrow(S) < 2) stop("Need at least 2 rows in 'S' to cluster.", call. = FALSE)
  if (any(!is.finite(S))) stop("S has non-finite values.", call. = FALSE)
  S
}

#' Enforce non-increasing survival curves over time
#' @keywords internal
.enforce_monotone <- function(S) {
  S <- as.matrix(S)
  n <- nrow(S); Q <- ncol(S)
  for (i in seq_len(n)) for (q in 2:Q) S[i, q] <- min(S[i, q], S[i, q - 1])
  S
}

#' Median smooth survival curves along the time axis
#' @keywords internal
.smooth_median <- function(S, width) {
  if (width < 3 || width %% 2 != 1) return(S)
  S <- as.matrix(S); n <- nrow(S); Q <- ncol(S)
  hw <- (width - 1) / 2; Sm <- S
  for (i in seq_len(n)) for (q in seq_len(Q)) {
    lo <- max(1, q - hw); hi <- min(Q, q + hw)
    Sm[i, q] <- stats::median(S[i, lo:hi])
  }
  Sm
}

#' Trapezoidal weights on a time grid
#' @keywords internal
.trap_weights <- function(times) {
  Q <- length(times)
  if (Q < 2) return(rep(1 / Q, Q))
  dt <- diff(times)
  w <- numeric(Q)
  w[1] <- dt[1]
  if (Q > 2) w[2:(Q - 1)] <- (dt[-1] + dt[-length(dt)]) / 2
  w[Q] <- dt[length(dt)]
  w / sum(w)
}

#' Apply time-point weights to curve features
#' @keywords internal
.weight_features <- function(S, w, distance) {
  if (distance == "L2") sweep(S, 2, sqrt(w), `*`) else sweep(S, 2, w, `*`)
}

#' Build a distance object from features
#' @keywords internal
.build_dist <- function(X, distance) {
  stats::dist(X, method = ifelse(distance == "L2", "euclidean", "manhattan"))
}


#' Plot medoid curves with ggplot2 (internal helper)
#' @keywords internal
plot_surv_medoids <- function(fit) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("tidyr", quietly = TRUE) ||
      !requireNamespace("dplyr", quietly = TRUE) ||
      !requireNamespace("scales", quietly = TRUE)) {
    stop("Install ggplot2, tidyr, dplyr, scales for plotting.", call. = FALSE)
  }
  med <- as.data.frame(fit$medoids)
  colnames(med) <- paste0("t", seq_along(fit$times))
  med$cluster <- factor(seq_len(nrow(med)))

  long <- tidyr::pivot_longer(
    med, dplyr::starts_with("t"),
    names_to = "gid", values_to = "S"
  ) |>
    dplyr::mutate(t = fit$times[as.integer(sub("t", "", gid))])

  ggplot2::ggplot(long, ggplot2::aes(t, S, color = cluster)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_y_continuous(labels = scales::percent_format()) +
    ggplot2::labs(
      title = paste0("PAM medoid curves (K = ", fit$K, ")"),
      x = "Time", y = "Survival", color = "Cluster"
    ) +
    ggplot2::theme_minimal()
}

#' Plot individual curves with ggplot2 (internal helper)
#' @keywords internal
plot_surv_samples <- function(S, times, clusters = NULL, alpha = 0.2) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("tidyr", quietly = TRUE) ||
      !requireNamespace("dplyr", quietly = TRUE)) {
    stop("Install ggplot2, tidyr, dplyr for plotting.", call. = FALSE)
  }
  S <- as.data.frame(S)
  colnames(S) <- paste0("t", seq_along(times))
  S$id <- seq_len(nrow(S))
  if (!is.null(clusters)) S$cluster <- factor(clusters)

  long <- tidyr::pivot_longer(
    S, dplyr::starts_with("t"),
    names_to = "gid", values_to = "S"
  ) |>
    dplyr::mutate(t = times[as.integer(sub("t", "", gid))])

  ggplot2::ggplot(long, ggplot2::aes(
    t, S, group = id,
    color = if (!is.null(clusters)) cluster else NULL
  )) +
    ggplot2::geom_line(alpha = alpha) +
    ggplot2::labs(x = "Time", y = "Survival", color = "Cluster") +
    ggplot2::theme_minimal()
}

#' Adjusted Rand Index (ARI)
#' @keywords internal
.ari <- function(a, b) {
  tab <- table(a, b); n <- sum(tab)
  if (n <= 1) return(NA_real_)
  c2 <- function(x) sum(choose(x, 2))
  sum_all  <- c2(tab)
  sum_rows <- c2(rowSums(tab))
  sum_cols <- c2(colSums(tab))
  expected <- (sum_rows * sum_cols) / choose(n, 2)
  max_idx  <- 0.5 * (sum_rows + sum_cols)
  denom <- max_idx - expected
  if (denom == 0) return(0)
  (sum_all - expected) / denom
}


#' Resampling-based stability assessment using ARI
#' @keywords internal
.cluster_stability <- function(
    S, times,
    B = 30, frac = 0.7,
    mode = c("subsample", "bootstrap"),
    jitter_sd = 0.001,
    weight_perturb = 0.001,
    seed = NULL,
    return_distribution = TRUE,
    ...   # forwarded to unsurv()
) {
  mode <- match.arg(mode)
  if (!is.null(seed)) set.seed(seed)

  n <- nrow(S)
  labs_list <- vector("list", B)
  idx_list  <- vector("list", B)

  pam_formals <- names(formals(unsurv))
  base_args <- list(...)
  base_args <- base_args[names(base_args) %in% pam_formals]
  base_args$seed <- NULL

  for (b in seq_len(B)) {
    m <- max(2, floor(frac * n))
    idx <- if (mode == "subsample") sort(sample.int(n, m, replace = FALSE))
    else sort(sample.int(n, m, replace = TRUE))
    S_b <- S[idx, , drop = FALSE]

    if (jitter_sd > 0) {
      S_b <- .clamp01(S_b + matrix(stats::rnorm(length(S_b), 0, jitter_sd), nrow(S_b)))
      S_b <- .enforce_monotone(S_b)
    }

    args_b <- base_args
    if (weight_perturb > 0) {
      w0 <- .trap_weights(times)
      z  <- stats::rexp(length(w0), 1); z <- z / sum(z)
      args_b$weights <- (1 - weight_perturb) * w0 + weight_perturb * z
    }

    fit_b <- do.call(unsurv, c(list(S = S_b, times = times), args_b))
    labs_list[[b]] <- fit_b$clusters
    idx_list[[b]]  <- idx
  }

  aris <- c()
  for (i in 1:(B - 1)) for (j in (i + 1):B) {
    overlap <- intersect(idx_list[[i]], idx_list[[j]])
    if (length(overlap) < 2) next
    a <- labs_list[[i]][match(overlap, idx_list[[i]])]
    b <- labs_list[[j]][match(overlap, idx_list[[j]])]
    aris <- c(aris, .ari(a, b))
  }

  if (!length(aris)) {
    if (return_distribution) return(list(mean = NA_real_, aris = numeric(0)))
    return(NA_real_)
  }
  if (return_distribution) list(mean = mean(aris, na.rm = TRUE), aris = aris)
  else mean(aris, na.rm = TRUE)
}


#' Plot stability distribution (internal helper)
#' @keywords internal
plot_stability <- function(stab) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install 'ggplot2'.", call. = FALSE)

  if (is.list(stab) && !is.null(stab$aris)) {
    aris <- stab$aris
  } else if (is.numeric(stab)) {
    aris <- stab
  } else {
    stop("Pass the list returned by unsurv_stability(..., return_distribution = TRUE).", call. = FALSE)
  }

  df <- data.frame(ARI = aris)
  ggplot2::ggplot(df, ggplot2::aes(ARI)) +
    ggplot2::geom_histogram(bins = 20) +
    ggplot2::geom_vline(xintercept = mean(aris), linetype = 2) +
    ggplot2::labs(
      title = sprintf("Stability (mean ARI = %.3f)", mean(aris)),
      x = "Adjusted Rand Index across resamples",
      y = "Count"
    ) +
    ggplot2::theme_minimal()
}

#' ggplot2 autoplot for unsurv objects
#'
#' This method requires the \pkg{ggplot2} package.
#'
#' @param object An object of class \code{"unsurv"}.
#' @param ... Unused.
#'
#' @return A ggplot object.
#' @export
autoplot.unsurv <- function(object, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install 'ggplot2' to use autoplot().", call. = FALSE)
  }

  times <- object$times
  med <- object$medoids
  K <- nrow(med)

  df <- data.frame(
    t = rep(times, times = K),
    S = as.vector(t(med)),
    cluster = factor(rep(seq_len(K), each = length(times)))
  )

  ggplot2::ggplot(df, ggplot2::aes(x = t, y = S, color = cluster)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::labs(
      title = paste0("unsurv medoid curves (K = ", object$K, ")"),
      x = "Time", y = "Survival", color = "Cluster"
    ) +
    ggplot2::theme_minimal()
}

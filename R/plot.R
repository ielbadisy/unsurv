plot_surv_medoids <- function(fit) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("tidyr", quietly = TRUE) ||
      !requireNamespace("dplyr", quietly = TRUE) ||
      !requireNamespace("scales", quietly = TRUE)) {
    stop("Install ggplot2, tidyr, dplyr, scales for plotting.")
  }
  med <- as.data.frame(fit$medoids)
  colnames(med) <- paste0("t", seq_along(fit$times))
  med$cluster <- factor(seq_len(nrow(med)))
  long <- tidyr::pivot_longer(med, dplyr::starts_with("t"),
                              names_to = "gid", values_to = "S") |>
    dplyr::mutate(t = fit$times[as.integer(sub("t", "", gid))])
  ggplot2::ggplot(long, ggplot2::aes(t, S, color = cluster)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_y_continuous(labels = scales::percent_format()) +
    ggplot2::labs(title = paste0("PAM medoid curves (K = ", fit$K, ")"),
                  x = "Time", y = "Survival", color = "Cluster") +
    ggplot2::theme_minimal()
}



plot_surv_samples <- function(S, times, clusters = NULL, alpha = 0.2) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("tidyr", quietly = TRUE) ||
      !requireNamespace("dplyr", quietly = TRUE)) {
    stop("Install ggplot2, tidyr, dplyr for plotting.")
  }
  S <- as.data.frame(S)
  colnames(S) <- paste0("t", seq_along(times))
  S$id <- seq_len(nrow(S))
  if (!is.null(clusters)) S$cluster <- factor(clusters)
  long <- tidyr::pivot_longer(S, dplyr::starts_with("t"),
                              names_to = "gid", values_to = "S") |>
    dplyr::mutate(t = times[as.integer(sub("t", "", gid))])
  ggplot2::ggplot(long, ggplot2::aes(t, S, group = id,
                                     color = if (!is.null(clusters)) cluster else NULL)) +
    ggplot2::geom_line(alpha = alpha) +
    ggplot2::labs(x = "Time", y = "Survival", color = "Cluster") +
    ggplot2::theme_minimal()
}


plot_stability <- function(stab){
  #if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install 'ggplot2'.")
  if (is.list(stab) && !is.null(stab$aris)) {
    aris <- stab$aris
  } else if (is.numeric(stab)) {
    aris <- stab
  } else stop("Pass the list returned by cluster_stability(..., return_distribution=TRUE).")
  df <- data.frame(ARI = aris)
  ggplot2::ggplot(df, ggplot2::aes(ARI)) +
    ggplot2::geom_histogram(bins = 20) +
    ggplot2::geom_vline(xintercept = mean(aris), linetype = 2) +
    ggplot2::labs(title = sprintf("Stability (mean ARI = %.3f)", mean(aris)),
                  x = "Adjusted Rand Index across resamples", y = "Count") +
    ggplot2::theme_minimal()
}



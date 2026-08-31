#' Per-cluster Kaplan-Meier median survival
#' @keywords internal
.km_cluster_medians <- function(lab, time, status) {
  lab <- factor(lab)
  sf <- survival::survfit(survival::Surv(time, status) ~ lab)
  tbl <- summary(sf)$table
  if (is.null(dim(tbl))) tbl <- matrix(tbl, nrow = 1, dimnames = list(levels(lab), names(tbl)))
  data.frame(
    cluster = levels(lab),
    n = as.integer(tbl[, "records"]),
    events = as.integer(tbl[, "events"]),
    median_survival = as.numeric(tbl[, "median"]),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

#' Compare cluster partitions against observed survival outcomes
#'
#' Summarizes and compares one or more cluster-label partitions of the same
#' individuals against observed time-to-event outcomes. This is intended for
#' comparing an \code{unsurv} curve-based partition against baseline
#' partitions (e.g., PAM on a scalar risk summary, or PAM on covariate PCA
#' scores), or for comparing the same partitioning rule applied to different
#' patient sets (e.g., a partition-defining set and a held-out validation
#' set) to check that survival separation generalizes.
#'
#' @param labels A named list of cluster-label vectors (integer or factor),
#'   each of the same length as \code{time}/\code{status}. Names are used as
#'   method labels; unnamed elements are labeled \code{"method1"},
#'   \code{"method2"}, etc.
#' @param time Numeric vector of observed follow-up times.
#' @param status Numeric/integer vector of event indicators (\code{1} = event,
#'   \code{0} = censored).
#' @param reference Name or integer index of the element of \code{labels}
#'   used as the reference partition for Adjusted Rand Index (ARI)
#'   agreement. Defaults to the first element.
#'
#' @details
#' Requires the \pkg{survival} package for Kaplan-Meier medians.
#'
#' For each partition, the Adjusted Rand Index quantifies agreement with the
#' reference partition. A log-rank test is deliberately not reported: when a
#' partition is itself fit to separate the curves (as \code{unsurv} and the
#' baselines are), a log-rank test against those same labels is circular and
#' close to guaranteed to be "significant," so it is not a fair basis for
#' comparing methods. Instead, per-cluster Kaplan-Meier medians are reported,
#' which is useful for checking whether the ordering of clusters by survival
#' (e.g., "cluster 2 has better survival than clusters 1 and 3") is preserved
#' across sets, such as a partition-defining set and an independent
#' validation set.
#'
#' @return An object of class \code{"unsurv_compare"} with elements:
#' \itemize{
#'   \item \code{summary}: one row per method with \code{K}, cluster-size
#'     range, and ARI against the reference partition.
#'   \item \code{cluster_summary}: one row per method/cluster with size,
#'     event count, and Kaplan-Meier median survival.
#'   \item \code{labels}, \code{time}, \code{status}, \code{reference}: the
#'     inputs, stored for plotting.
#' }
#'
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 120
#'   time <- stats::rexp(n, 0.1)
#'   status <- sample(0:1, n, TRUE)
#'   labs <- list(
#'     unsurv_curve = sample(1:3, n, TRUE),
#'     scalar_risk = sample(1:3, n, TRUE)
#'   )
#'   cmp <- unsurv_compare(labs, time, status)
#'   print(cmp)
#' }
#' @export
unsurv_compare <- function(labels, time, status, reference = 1) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Install the 'survival' package to use unsurv_compare().", call. = FALSE)
  }
  if (!is.list(labels) || length(labels) < 1) {
    stop("'labels' must be a non-empty list of cluster-label vectors.", call. = FALSE)
  }
  n <- length(time)
  if (length(status) != n) stop("'time' and 'status' must have the same length.", call. = FALSE)
  for (l in labels) {
    if (length(l) != n) stop("Every element of 'labels' must have length(time) entries.", call. = FALSE)
  }

  nm <- names(labels)
  if (is.null(nm) || any(nm == "")) nm <- paste0("method", seq_along(labels))
  names(labels) <- nm

  ref_idx <- if (is.character(reference)) match(reference, nm) else as.integer(reference)
  if (is.na(ref_idx) || ref_idx < 1 || ref_idx > length(labels)) {
    stop("'reference' must identify an element of 'labels'.", call. = FALSE)
  }
  ref_lab <- labels[[ref_idx]]

  rows <- vector("list", length(labels))
  cluster_rows <- vector("list", length(labels))
  for (i in seq_along(labels)) {
    lab <- factor(labels[[i]])
    tab <- table(lab)

    rows[[i]] <- data.frame(
      method = nm[i],
      K = length(tab),
      min_size = as.integer(min(tab)),
      max_size = as.integer(max(tab)),
      ari_ref = .ari(labels[[i]], ref_lab),
      stringsAsFactors = FALSE
    )

    cs <- .km_cluster_medians(lab, time, status)
    cs$method <- nm[i]
    cluster_rows[[i]] <- cs
  }

  structure(
    list(
      summary = do.call(rbind, rows),
      cluster_summary = do.call(rbind, cluster_rows),
      labels = labels,
      time = as.numeric(time),
      status = as.numeric(status),
      reference = nm[ref_idx]
    ),
    class = "unsurv_compare"
  )
}

#' Print a partition comparison
#'
#' @param x An object of class \code{"unsurv_compare"}.
#' @param ... Unused.
#'
#' @return Invisibly returns \code{x}.
#' @export
print.unsurv_compare <- function(x, ...) {
  cat("unsurv_compare: partition comparison (reference = ", x$reference, ")\n\n", sep = "")
  print(x$summary, row.names = FALSE)
  invisible(x)
}

#' Kaplan-Meier plot for a partition comparison
#'
#' Plots Kaplan-Meier survival curves for each cluster, faceted by
#' comparison method, for an object returned by \code{\link{unsurv_compare}}.
#'
#' @param object An object of class \code{"unsurv_compare"}.
#' @param ... Unused.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' if (requireNamespace("survival", quietly = TRUE) &&
#'     requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 120
#'   time <- stats::rexp(n, 0.1)
#'   status <- sample(0:1, n, TRUE)
#'   labs <- list(
#'     unsurv_curve = sample(1:3, n, TRUE),
#'     scalar_risk = sample(1:3, n, TRUE)
#'   )
#'   cmp <- unsurv_compare(labs, time, status)
#'   ggplot2::autoplot(cmp)
#' }
#' @export
#' @method autoplot unsurv_compare
autoplot.unsurv_compare <- function(object, ...) {
  if (!requireNamespace("survival", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install 'survival' and 'ggplot2' to use autoplot().", call. = FALSE)
  }

  dfs <- lapply(names(object$labels), function(m) {
    lab <- factor(object$labels[[m]])
    sf <- survival::survfit(survival::Surv(object$time, object$status) ~ lab)
    data.frame(
      time = sf$time,
      survival = sf$surv,
      cluster = sub("^lab=", "", rep(names(sf$strata), sf$strata)),
      method = m,
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, dfs)

  ggplot2::ggplot(df, ggplot2::aes(time, survival, color = cluster)) +
    ggplot2::geom_step(linewidth = 0.8) +
    ggplot2::facet_wrap(~method) +
    ggplot2::labs(x = "Time", y = "Kaplan-Meier survival", color = "Cluster") +
    ggplot2::theme_minimal()
}

#' @rdname autoplot.unsurv_compare
#' @param x An object of class \code{"unsurv_compare"}.
#' @export
#' @method plot unsurv_compare
plot.unsurv_compare <- function(x, ...) {
  p <- autoplot.unsurv_compare(x, ...)
  print(p)
  invisible(x)
}

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
  print(x$K)
  print(x$silhouette_mean)
  cat("$size\n\n")
  print(x$size)
  invisible(x)
}

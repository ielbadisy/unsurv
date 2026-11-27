#' Plot medoid survival curves from an unsurv fit
#'
#' Produces a base R plot of the cluster medoid survival curves stored in the
#' fitted object.
#'
#' @param x An object of class \code{"unsurv"}.
#' @param col Line colours, recycled over clusters; also used in the legend.
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
plot.unsurv <- function(x, col = seq_len(nrow(x$medoids)), ...) {
  stopifnot(inherits(x, "unsurv"))
  times <- x$times
  med <- x$medoids
  if (is.null(med) || nrow(med) < 1L) stop("No medoids found in 'x$medoids'.", call. = FALSE)
  if (length(times) != ncol(med)) stop("Dimension mismatch: length(x$times) must equal ncol(x$medoids).", call. = FALSE)

  K <- nrow(med)
  col <- rep_len(col, K)
  graphics::matplot(times, t(med), type = "l", lty = 1, col = col,
                    xlab = "Time", ylab = "Survival", ...)
  graphics::legend("topright", legend = paste0("Cluster ", seq_len(K)),
                   lty = 1, col = col, bty = "n")
  invisible(x)
}

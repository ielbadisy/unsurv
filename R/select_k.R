#' Compare numbers of clusters
#'
#' Fits \code{\link{unsurv}} for each candidate \code{K} and returns, side by side, the mean
#' silhouette width, the resampling stability of \code{\link{unsurv_stability}}, the
#' prediction strength of Tibshirani and Walther (2005) and the smallest cluster size.
#'
#' @param S Numeric matrix of survival probabilities (\eqn{n \times m}), columns matching
#'   \code{times}.
#' @param times Numeric vector of time grid points (length \eqn{m}).
#' @param K Integer vector of candidate numbers of clusters (each \eqn{\ge 2}).
#' @param B Integer; number of resamples for the stability index.
#' @param splits Integer; number of random half splits for the prediction strength.
#' @param seed Optional integer seed; \code{NULL} uses the current random number stream.
#' @param ... Further arguments passed to \code{\link{unsurv}} (for example \code{distance},
#'   \code{weights} or \code{eps_jitter}).
#'
#' @details The prediction strength of a split into halves \eqn{A} and \eqn{B} clusters both
#'   halves, assigns the patients of \eqn{B} to the nearest medoid of the clustering of
#'   \eqn{A} (\code{\link{predict.unsurv}}), and takes, over the clusters of \eqn{B}, the smallest
#'   share of pairs of patients of the same cluster that the medoids of \eqn{A} also put together.
#'   The column reports its mean over \code{splits} random splits. Tibshirani and Walther suggest
#'   choosing the largest \code{K} whose prediction strength exceeds about 0.8. The silhouette
#'   and the stability index tend to favour small \code{K}, and none of the columns can say that
#'   there is no cluster structure; read the table, not only its maxima.
#'
#' @return A data frame with columns \code{K}, \code{silhouette}, \code{stability},
#'   \code{pred_strength} and \code{min_size}.
#'
#' @references
#' Tibshirani, R. and Walther, G. (2005). Cluster validation by prediction strength.
#' \emph{Journal of Computational and Graphical Statistics}, 14(3), 511--528.
#' \doi{10.1198/106186005X59243}
#'
#' @examples
#' if (requireNamespace("cluster", quietly = TRUE)) {
#'   set.seed(1)
#'   times <- seq(0, 4, length.out = 15)
#'   grp <- sample(1:3, 60, TRUE)
#'   S <- t(sapply(c(0.15, 0.5, 1.2)[grp] * exp(rnorm(60, 0, 0.1)), function(r) exp(-r * times)))
#'   unsurv_select_k(S, times, K = 2:4, B = 5, splits = 3, seed = 1)
#' }
#' @export
unsurv_select_k <- function(S, times, K = 2:6, B = 20, splits = 10, seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  fits <- lapply(K, function(k) unsurv(S, times, K = k, ...))
  data.frame(
    K = as.integer(K),
    silhouette = vapply(fits, `[[`, numeric(1), "silhouette_mean"),
    stability = vapply(fits, function(f)
      unsurv_stability(S, times, f, B = B, return_distribution = FALSE), numeric(1)),
    pred_strength = vapply(K, function(k)
      mean(replicate(splits, .pred_strength(S, times, k, ...))), numeric(1)),
    min_size = vapply(fits, function(f) min(tabulate(f$clusters)), integer(1))
  )
}

#' Prediction strength of one random half split
#' @keywords internal
.pred_strength <- function(S, times, K, ...) {
  i <- sample.int(nrow(S), nrow(S) %/% 2)
  a <- unsurv(S[i, , drop = FALSE], times, K = K, ...)
  b <- unsurv(S[-i, , drop = FALSE], times, K = K, ...)$clusters
  p <- predict(a, S[-i, , drop = FALSE])
  min(vapply(seq_len(K), function(k) {
    j <- which(b == k); nj <- length(j)
    if (nj < 2) return(1)
    (sum(outer(p[j], p[j], "==")) - nj) / (nj * (nj - 1))
  }, numeric(1)))
}

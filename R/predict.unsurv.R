#' Predict cluster membership for new survival curves
#'
#' Assigns new survival-probability curves to clusters using the medoids from a
#' fitted \code{unsurv} object. New curves are preprocessed using the same
#' weighting, optional monotonic enforcement, smoothing, and standardization
#' parameters as the fitted model.
#'
#' @param object An object of class \code{"unsurv"}, returned by \code{\link{unsurv}}.
#' @param newdata Numeric matrix of survival probabilities with shape
#'   \eqn{n_{new} \times m}, where columns correspond to the same time grid
#'   used during fitting.
#' @param clamp Logical; if \code{TRUE}, clamps values to \eqn{[0,1]} before
#'   preprocessing.
#' @param ... Unused. Included for compatibility with the generic.
#'
#' @details
#' Cluster assignment is performed by computing distances between the new curves
#' and the stored medoid curves in the weighted feature space defined during
#' fitting. The distance metric (\code{"L1"} or \code{"L2"}) and any
#' standardization parameters are reused from the fitted model.
#'
#' @return An integer vector of cluster labels of length \code{nrow(newdata)},
#'   taking values in \code{1, ..., object$K}.
#'
#' @examples
#' if (requireNamespace("cluster", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 60; Q <- 40
#'   times <- seq(0, 5, length.out = Q)
#'   grp <- sample(1:2, n, TRUE)
#'   rates <- c(0.2, 0.6)
#'   S <- sapply(times, function(t) exp(-rates[grp] * t))
#'   S <- S + matrix(stats::rnorm(n * Q, 0, 0.02), nrow = n)
#'
#'   fit <- unsurv(S, times, K = 2)
#'
#'   # predict cluster membership for first 5 curves
#'   predict(fit, S[1:5, ])
#' }
#'
#' @export
#' @method predict unsurv
#' @importFrom stats predict
predict.unsurv <- function(object, newdata, clamp = TRUE, ...) {

  stopifnot(inherits(object, "unsurv"))

  S_new <- as.matrix(newdata)

  if (ncol(S_new) != length(object$times)) {
    stop(
      "newdata must have the same number of columns/time grid as the fit.",
      call. = FALSE
    )
  }

  if (any(!is.finite(S_new))) {
    stop("newdata contains non-finite values.", call. = FALSE)
  }

  # preprocessing identical to fitting
  if (isTRUE(clamp)) {
    S_new <- .clamp01(S_new)
  }

  if (isTRUE(object$enforce_monotone)) {
    S_new <- .enforce_monotone(S_new)
  }

  if (is.numeric(object$smooth_median_width) &&
      object$smooth_median_width >= 3 &&
      object$smooth_median_width %% 2 == 1) {

    S_new <- .smooth_median(S_new, object$smooth_median_width)
  }

  # recreate weighted feature space
  X_new <- .weight_features(S_new, object$weights, object$distance)
  X_med <- .weight_features(object$medoids, object$weights, object$distance)

  X_new <- as.matrix(X_new)
  X_med <- as.matrix(X_med)

  # apply stored standardization if used during fitting
  if (isTRUE(object$standardize_cols)) {

    mu <- object$center
    sdv <- object$scale

    if (is.null(mu) || is.null(sdv)) {
      stop(
        "fit was standardized but center/scale not stored.",
        call. = FALSE
      )
    }

    X_new <- sweep(X_new, 2, mu, "-")
    X_new <- sweep(X_new, 2, sdv, "/")

    X_med <- sweep(X_med, 2, mu, "-")
    X_med <- sweep(X_med, 2, sdv, "/")
  }

  # compute distances to medoids
  if (object$distance == "L2") {

    rn <- rowSums(X_new^2)
    rm <- rowSums(X_med^2)

    d2 <- outer(rn, rm, "+") - 2 * (X_new %*% t(X_med))
    d2[d2 < 0] <- 0

    d <- sqrt(d2)

  } else {

    d <- matrix(NA_real_, nrow(X_new), nrow(X_med))

    for (k in seq_len(nrow(X_med))) {

      d[, k] <- rowSums(
        abs(X_new - matrix(X_med[k, ], nrow(X_new), byrow = TRUE))
      )
    }
  }

  # assign cluster = nearest medoid
  cluster <- max.col(-d)

  as.integer(cluster)
}

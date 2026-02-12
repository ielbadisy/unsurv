#' unsurv: Unsupervised clustering of individualized survival curves
#'
#' The \pkg{unsurv} package provides tools for unsupervised clustering of
#' individualized survival curves using medoid-based clustering (PAM).
#'
#' It is designed for settings where each individual is represented by a survival
#' probability curve evaluated on a common time grid, such as predictions from:
#'
#' \itemize{
#'   \item Kaplan–Meier estimates,
#'   \item Cox models,
#'   \item parametric survival models,
#'   \item deep learning survival models,
#'   \item or other individualized survival predictors.
#' }
#'
#' Core features include:
#'
#' \itemize{
#'   \item PAM clustering using weighted L1 or L2 distances,
#'   \item automatic cluster selection via silhouette width,
#'   \item optional monotonicity enforcement,
#'   \item optional median smoothing,
#'   \item prediction of cluster membership for new curves,
#'   \item stability assessment via resampling and Adjusted Rand Index,
#'   \item base R and ggplot2 visualization methods.
#' }
#'
#' Main functions:
#'
#' \itemize{
#'   \item \code{\link{unsurv}} — fit clustering model
#'   \item \code{\link{predict.unsurv}} — predict cluster membership
#'   \item \code{\link{plot.unsurv}} — plot medoid curves
#'   \item \code{\link{summary.unsurv}} — summarize clustering
#'   \item \code{\link{unsurv_stability}} — assess stability
#' }
#'
#' @docType package
#' @name unsurv-package
#'
#' @author
#' Imad EL BADISY
#'
#' @references
#' Kaufman, L., & Rousseeuw, P. J. (1990).
#' \emph{Finding Groups in Data: An Introduction to Cluster Analysis}.
#' Wiley.
#'
#' @examples
#' if (requireNamespace("cluster", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 10
#'   times <- seq(0, 5, length.out = 40)
#'   rates <- sample(c(0.2, 0.6), n, TRUE)
#'   S <- sapply(times, function(t) exp(-rates * t))
#'
#'   fit <- unsurv(S, times, K = 2)
#'   plot(fit)
#' }
"_PACKAGE"

## ---- utilities ----
.clamp01 <- function(M) { M[M < 0] <- 0; M[M > 1] <- 1; M }

.check_inputs <- function(S, times) {
  S <- as.matrix(S)
  if (!is.numeric(times) || any(!is.finite(times))) stop("'times' must be finite numeric.")
  if (is.unsorted(times, strictly = TRUE)) stop("'times' must be strictly increasing.")
  if (ncol(S) != length(times)) stop("ncol(S) must equal length(times).")
  if (any(!is.finite(S))) stop("S has non-finite values.")
  S
}

.enforce_monotone <- function(S) {
  S <- as.matrix(S)
  n <- nrow(S); Q <- ncol(S)
  for (i in seq_len(n)) for (q in 2:Q) S[i, q] <- min(S[i, q], S[i, q - 1])
  S
}

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

.trap_weights <- function(times) {
  Q <- length(times)
  if (Q < 2) return(rep(1, Q))
  dt <- diff(times)
  w <- numeric(Q)
  w[1] <- dt[1]
  if (Q > 2) w[2:(Q - 1)] <- (dt[-1] + dt[-length(dt)]) / 2
  w[Q] <- dt[length(dt)]
  w / sum(w)
}

.weight_features <- function(S, w, distance) {
  if (distance == "L2") sweep(S, 2, sqrt(w), `*`) else sweep(S, 2, w, `*`)
}

.build_dist <- function(X, distance) {
  stats::dist(X, method = ifelse(distance == "L2", "euclidean", "manhattan"))
}

make_toy_curves <- function(n = 30, Q = 40, seed = 1) {
  set.seed(seed)
  times <- seq(0, 5, length.out = Q)
  grp <- sample(1:2, n, TRUE)
  rates <- c(0.2, 0.7)

  S <- sapply(times, function(t) exp(-rates[grp] * t))
  S <- S + matrix(stats::rnorm(n * Q, 0, 0.01), nrow = n)
  S[S < 0] <- 0
  S[S > 1] <- 1

  list(S = S, times = times, grp = grp)
}

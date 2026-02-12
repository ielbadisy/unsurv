test_that("unsurv fits and returns expected structure", {
  skip_if_not_installed("cluster")

  dat <- make_toy_curves(n = 30, Q = 30, seed = 1)
  fit <- unsurv(dat$S, dat$times, K = 2, eps_jitter = 0)

  expect_s3_class(fit, "unsurv")
  expect_true(is.integer(fit$clusters))
  expect_equal(length(fit$clusters), nrow(dat$S))
  expect_equal(fit$K, 2L)
  expect_equal(length(fit$times), ncol(dat$S))
  expect_equal(dim(fit$medoids), c(2, ncol(dat$S)))

  # sanity: silhouette in [-1, 1]
  expect_true(is.finite(fit$silhouette_mean))
  expect_gte(fit$silhouette_mean, -1)
  expect_lte(fit$silhouette_mean, 1)
})

test_that("unsurv enforces monotonicity when requested", {
  skip_if_not_installed("cluster")

  dat <- make_toy_curves(n = 20, Q = 25, seed = 2)
  # Introduce a monotonicity violation
  dat$S[1, 10] <- 0.1
  dat$S[1, 11] <- 0.9

  fit <- unsurv(dat$S, dat$times, K = 2, enforce_monotone = TRUE, eps_jitter = 0)

  # medoids should be non-increasing if enforce_monotone=TRUE was applied before PAM
  med <- fit$medoids
  diffs <- apply(med, 1, diff)
  expect_true(all(diffs <= 1e-12))
})

test_that("unsurv errors on bad inputs", {
  skip_if_not_installed("cluster")

  dat <- make_toy_curves(n = 10, Q = 10, seed = 3)

  expect_error(unsurv(dat$S, rev(dat$times), K = 2), "strictly increasing")
  expect_error(unsurv(dat$S[, -1, drop = FALSE], dat$times, K = 2), "ncol\\(S\\) must equal")
  expect_error(unsurv(dat$S[1, , drop = FALSE], dat$times, K = 2), "at least 2 rows")
  expect_error(unsurv(dat$S, dat$times, K = 0), "must be")
  expect_error(unsurv(dat$S, dat$times, K = 1000), "cannot exceed")
})

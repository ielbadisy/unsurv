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

test_that("unsurv can leave curves non-monotone when requested", {
  skip_if_not_installed("cluster")

  S <- rbind(
    c(1, 0.9, 0.8, 0.7, 0.6),
    c(1, 0.85, 0.95, 0.9, 0.75), # has upward bump
    c(1, 0.9, 0.78, 0.65, 0.55)
  )
  times <- seq_along(S[1, ])

  fit <- unsurv(S, times, K = 2, enforce_monotone = FALSE, eps_jitter = 0)
  expect_false(fit$enforce_monotone)

  # one medoid should retain the upward bump
  bump_sizes <- apply(fit$medoids, 1, function(x) max(diff(x)))
  expect_gt(max(bump_sizes), 0)
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

test_that("unsurv handles distance, weights, and automatic K selection", {
  skip_if_not_installed("cluster")

  dat <- make_toy_curves(n = 50, Q = 30, seed = 8)
  custom_w <- seq_len(ncol(dat$S))

  fit <- unsurv(
    dat$S, dat$times,
    K = NULL, K_max = 5,
    distance = "L1",
    weights = custom_w,
    standardize_cols = TRUE,
    eps_jitter = 0
  )

  expect_equal(fit$distance, "L1")
  expect_true(abs(sum(fit$weights) - 1) < 1e-8)
  expect_equal(length(fit$weights), ncol(dat$S))
  expect_true(fit$K >= 2 && fit$K <= 5)
  expect_true(is.finite(fit$silhouette_mean))
})

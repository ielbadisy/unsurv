test_that("unsurv_stability returns expected output", {
  skip_if_not_installed("cluster")

  dat <- make_toy_curves(n = 30, Q = 25, seed = 6)
  fit <- unsurv(dat$S, dat$times, K = 2, eps_jitter = 0)

  stab <- unsurv_stability(
    dat$S, dat$times, fit,
    B = 6, frac = 0.7, mode = "subsample",
    jitter_sd = 0.01, weight_perturb = 0, eps_jitter = 0,
    return_distribution = TRUE
  )

  expect_true(is.list(stab))
  expect_true(all(c("mean", "aris") %in% names(stab)))
  expect_true(is.numeric(stab$mean) && length(stab$mean) == 1)
  expect_true(is.numeric(stab$aris))
  expect_gte(stab$mean, -1)
  expect_lte(stab$mean, 1)
})

test_that("unsurv_stability errors on dimension mismatch", {
  skip_if_not_installed("cluster")

  dat <- make_toy_curves(n = 20, Q = 20, seed = 7)
  fit <- unsurv(dat$S, dat$times, K = 2, eps_jitter = 0)

  expect_error(
    unsurv_stability(dat$S[, -1, drop = FALSE], dat$times, fit, B = 3),
    "dimension mismatch"
  )
})

test_that("unsurv_stability refits with the weights of fit", {
  skip_if_not_installed("cluster")
  dat <- make_toy_curves(n = 20, Q = 10, seed = 8)
  w <- c(rep(0, 9), 1)
  fit <- unsurv(dat$S, dat$times, K = 2, weights = w, eps_jitter = 0)
  seen <- list()
  local_mocked_bindings(unsurv = function(S, times, K = NULL, weights = NULL, ...) {
    seen <<- c(seen, list(weights))
    list(clusters = rep(1:2, length.out = nrow(S)))
  })
  unsurv_stability(dat$S, dat$times, fit, B = 3, weight_perturb = 0)
  expect_length(seen, 3); for (s in seen) expect_equal(s, w)
})

test_that("unsurv_stability is reproducible with seed", {
  skip_if_not_installed("cluster")
  dat <- make_toy_curves(n = 30, Q = 15, seed = 9)
  fit <- unsurv(dat$S, dat$times, K = 2)
  a <- unsurv_stability(dat$S, dat$times, fit, B = 4, seed = 3)
  b <- unsurv_stability(dat$S, dat$times, fit, B = 4, seed = 3)
  expect_identical(a, b)
})

test_that("unsurv_stability resamples patients only by default", {
  skip_if_not_installed("cluster")
  dat <- make_toy_curves(n = 20, Q = 10, seed = 10)
  fit <- unsurv(dat$S, dat$times, K = 2, eps_jitter = 0.003)
  seen <- list()
  local_mocked_bindings(unsurv = function(S, times, K = NULL, weights = NULL, eps_jitter = NULL, ...) {
    seen <<- c(seen, list(list(S = S, w = weights, eps = eps_jitter)))
    list(clusters = rep(1:2, length.out = nrow(S)))
  })
  unsurv_stability(dat$S, dat$times, fit, B = 3)
  S0 <- unsurv:::.clamp01(dat$S)
  for (s in seen) {
    expect_equal(nrow(s$S), nrow(S0)); expect_equal(s$eps, 0.003); expect_equal(s$w, fit$weights)
    expect_true(all(apply(s$S, 1, function(r) any(colSums(abs(t(S0) - r)) == 0))))
  }
})

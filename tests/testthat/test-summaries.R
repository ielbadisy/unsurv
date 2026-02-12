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

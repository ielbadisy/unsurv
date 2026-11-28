test_that("unsurv_select_k returns one row per K with valid indices", {
  skip_if_not_installed("cluster")
  dat <- make_toy_curves(n = 40, Q = 15, seed = 11)
  tab <- unsurv_select_k(dat$S, dat$times, K = 2:4, B = 4, splits = 3, seed = 1)
  expect_equal(tab$K, 2:4)
  expect_equal(tab$silhouette, vapply(2:4, function(k)
    unsurv(dat$S, dat$times, K = k, eps_jitter = 0)$silhouette_mean, numeric(1)), tolerance = 0.05)
  expect_true(all(tab$pred_strength >= 0 & tab$pred_strength <= 1))
  expect_true(all(tab$stability <= 1)); expect_true(all(tab$min_size >= 1))
  expect_identical(tab, unsurv_select_k(dat$S, dat$times, K = 2:4, B = 4, splits = 3, seed = 1))
})

test_that("prediction strength is 1 for well separated groups at the true K", {
  skip_if_not_installed("cluster")
  set.seed(2); times <- seq(0, 4, length.out = 10); g <- rep(1:2, 20)
  S <- t(sapply(c(0.1, 2)[g], function(r) exp(-r * times)))
  expect_equal(mean(replicate(5, unsurv:::.pred_strength(S, times, 2))), 1)
})

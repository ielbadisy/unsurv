test_that("predict.unsurv returns valid cluster labels", {
  skip_if_not_installed("cluster")

  dat <- make_toy_curves(n = 30, Q = 35, seed = 4)
  fit <- unsurv(dat$S, dat$times, K = 2, eps_jitter = 0)

  pred <- predict(fit, dat$S[1:7, , drop = FALSE])
  expect_type(pred, "integer")
  expect_equal(length(pred), 7)
  expect_true(all(pred %in% 1:fit$K))

  # dimension mismatch should error
  expect_error(predict(fit, dat$S[1:7, -1, drop = FALSE]), "same number of columns")
})

test_that("predict.unsurv works with standardize_cols=TRUE", {
  skip_if_not_installed("cluster")

  dat <- make_toy_curves(n = 40, Q = 30, seed = 5)
  fit <- unsurv(dat$S, dat$times, K = 2, standardize_cols = TRUE, eps_jitter = 0)

  pred <- predict(fit, dat$S[1:10, , drop = FALSE])
  expect_true(all(pred %in% 1:fit$K))
})

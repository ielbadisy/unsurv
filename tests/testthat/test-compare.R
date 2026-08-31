test_that("unsurv_compare returns expected structure", {
  skip_if_not_installed("survival")

  set.seed(11)
  n <- 90
  time <- stats::rexp(n, 0.15)
  status <- sample(0:1, n, TRUE)
  labs <- list(
    unsurv_curve = sample(1:3, n, TRUE),
    scalar_risk = sample(1:3, n, TRUE)
  )

  cmp <- unsurv_compare(labs, time, status)

  expect_s3_class(cmp, "unsurv_compare")
  expect_equal(nrow(cmp$summary), 2)
  expect_true(all(c("method", "K", "min_size", "max_size", "ari_ref") %in% names(cmp$summary)))
  expect_false("logrank_p" %in% names(cmp$summary))
  expect_equal(cmp$summary$ari_ref[cmp$summary$method == "unsurv_curve"], 1)
  expect_true(all(cmp$cluster_summary$method %in% names(labs)))
  expect_equal(cmp$reference, "unsurv_curve")
})

test_that("unsurv_compare validates inputs", {
  skip_if_not_installed("survival")

  n <- 20
  time <- stats::rexp(n, 0.1)
  status <- sample(0:1, n, TRUE)

  expect_error(unsurv_compare(list(), time, status), "non-empty")
  expect_error(
    unsurv_compare(list(a = sample(1:2, n - 1, TRUE)), time, status),
    "length\\(time\\)"
  )
  expect_error(
    unsurv_compare(list(a = sample(1:2, n, TRUE)), time, status, reference = "missing"),
    "reference"
  )
})

test_that("autoplot.unsurv_compare returns a ggplot", {
  skip_if_not_installed("survival")
  skip_if_not_installed("ggplot2")

  set.seed(12)
  n <- 60
  time <- stats::rexp(n, 0.1)
  status <- sample(0:1, n, TRUE)
  labs <- list(a = sample(1:2, n, TRUE), b = sample(1:2, n, TRUE))
  cmp <- unsurv_compare(labs, time, status)

  p <- ggplot2::autoplot(cmp)
  expect_s3_class(p, "gg")
})

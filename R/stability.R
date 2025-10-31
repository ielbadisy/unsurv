.ari <- function(a, b) {
  tab <- table(a, b); n <- sum(tab)
  if (n <= 1) return(NA_real_)
  c2 <- function(x) sum(choose(x, 2))
  sum_all  <- c2(tab)
  sum_rows <- c2(rowSums(tab))
  sum_cols <- c2(colSums(tab))
  expected <- (sum_rows * sum_cols) / choose(n, 2)
  max_idx  <- 0.5 * (sum_rows + sum_cols)
  denom <- max_idx - expected
  if (denom == 0) return(0)
  (sum_all - expected) / denom
}

## ---- stability via ARI over resamples ----
cluster_stability <- function(
    S, times,
    B = 30, frac = 0.7,
    mode = c("subsample", "bootstrap"),
    jitter_sd = 0.001,              # curve-level additive noise before monotone clamp
    weight_perturb = 0.001,         # 0..1: blend base trapezoid weights with random simplex
    seed = NULL,
    return_distribution = TRUE,
    ...                         # forwarded to fit_unsurv
){
  mode <- match.arg(mode)
  if (!is.null(seed)) set.seed(seed)

  n <- nrow(S)
  labs_list <- vector("list", B)
  idx_list  <- vector("list", B)

  # keep only formal args for the fitter
  pam_formals <- names(formals(fit_unsurv))
  base_args <- list(...)
  base_args <- base_args[names(base_args) %in% pam_formals]
  base_args$seed <- NULL   # ensure variability across resamples

  for (b in seq_len(B)) {
    m <- max(2, floor(frac * n))
    idx <- if (mode == "subsample") sort(sample.int(n, m, replace = FALSE))
    else sort(sample.int(n, m, replace = TRUE))
    S_b <- S[idx, , drop = FALSE]

    # curve jitter
    if (jitter_sd > 0) {
      S_b <- .clamp01(S_b + matrix(rnorm(length(S_b), 0, jitter_sd), nrow(S_b)))
      S_b <- .enforce_monotone(S_b)
    }

    # weight perturbation
    args_b <- base_args
    if (weight_perturb > 0) {
      w0 <- .trap_weights(times)
      z  <- rexp(length(w0), 1); z <- z / sum(z)
      args_b$weights <- (1 - weight_perturb) * w0 + weight_perturb * z
    }

    fit_b <- do.call(fit_unsurv, c(list(S = S_b, times = times), args_b))
    labs_list[[b]] <- fit_b$clusters
    idx_list[[b]]  <- idx
  }

  ## ARI on overlaps across all pairs
  aris <- c()
  for (i in 1:(B - 1)) for (j in (i + 1):B) {
    overlap <- intersect(idx_list[[i]], idx_list[[j]])
    if (length(overlap) < 2) next
    a <- labs_list[[i]][match(overlap, idx_list[[i]])]
    b <- labs_list[[j]][match(overlap, idx_list[[j]])]
    aris <- c(aris, .ari(a, b))
  }

  if (!length(aris)) {
    if (return_distribution) {
      return(list(mean = NA_real_, aris = numeric(0)))
    } else {
      return(NA_real_)
    }
  }
  if (return_distribution) list(mean = mean(aris, na.rm = TRUE), aris = aris)
  else mean(aris, na.rm = TRUE)
}

## ---- computing stability ----
unsurv_stability <- function(
    S, times, fit,
    B = 30, frac = 0.5, mode = c("bootstrap","subsample"),
    jitter_sd = 0.01,          # curve-space noise
    weight_perturb = 0.30,     # 0..1 blend of base trapezoid weights
    eps_jitter = 0.02,         # feature-space jitter inside PAM
    return_distribution = TRUE
){
  mode <- match.arg(mode)
  stopifnot(inherits(fit, "unsurv"))
  if (ncol(S) != length(times)) stop("S/times dimension mismatch.")

  cluster_stability(
    S, times,
    B = B, frac = frac, mode = mode,
    jitter_sd = jitter_sd, weight_perturb = weight_perturb,
    K = fit$K,
    distance = fit$distance,
    standardize_cols = isTRUE(fit$standardize_cols),
    enforce_monotone = isTRUE(fit$enforce_monotone),
    smooth_median_width = ifelse(is.null(fit$smooth_median_width), 0, fit$smooth_median_width),
    eps_jitter = eps_jitter,
    seed = NULL,
    return_distribution = return_distribution
  )
}

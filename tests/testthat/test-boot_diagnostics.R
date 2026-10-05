data(pop)
data(samp)

# a rescaled elevation keeps the fits well conditioned
elev_m <- mean(pop$elev)
elev_s <- sd(pop$elev)
pop$elev_s <- (pop$elev - elev_m) / elev_s
samp$elev_s <- (samp$elev - elev_m) / elev_s

real_fit_zi <- fit_zi

# Replace fit_zi() so that chosen calls fail or warn. Call 1 is the original
# model fit; calls 2 to B + 1 are the bootstrap refits, in order.
mock_fit_zi <- function(fail_calls = integer(0), warn_calls = integer(0)) {
  n <- 0
  function(...) {
    n <<- n + 1
    if (n %in% fail_calls) stop("forced failure")
    out <- real_fit_zi(...)
    if (n %in% warn_calls) warning("forced convergence warning")
    out
  }
}

fit_diag <- function(seed = 1, ...) {
  set.seed(seed)
  suppressMessages(
    saeczi(samp,
           pop,
           lin_formula = DRYBIO_AG_TPA_live_ADJ ~ tcc16 + elev_s,
           log_formula = DRYBIO_AG_TPA_live_ADJ ~ elev_s,
           domain_level = "COUNTYFIPS",
           mse_est = TRUE,
           parallel = FALSE,
           ...)
  )
}

test_that("summarise_boot_mse matches sd / sqrt(K) and one-way ANOVA", {
  set.seed(3)
  K <- 8
  m <- 4
  pop_id <- rep(seq_len(K), each = m)
  sq_err <- matrix(rexp(3 * K * m), 3, K * m) *
    matrix(rep(rexp(K, 0.5), each = m), 3, K * m, byrow = TRUE)

  out <- summarise_boot_mse(sq_err, pop_id, K)

  pop_means <- t(apply(sq_err, 1, \(r) tapply(r, pop_id, mean)))
  expect_equal(out$by_domain$mse, rowMeans(sq_err))
  expect_equal(out$by_domain$n_boot_used, rep(K * m, 3))
  expect_equal(out$by_domain$mse_se, apply(pop_means, 1, sd) / sqrt(K))

  vc <- t(sapply(1:3, \(j) {
    a <- anova(lm(sq_err[j, ] ~ factor(pop_id)))
    ms_within <- a["Residuals", "Mean Sq"]
    c((a[1, "Mean Sq"] - ms_within) / m, ms_within)
  }))
  expect_equal(unname(as.matrix(out$var_comp)), unname(vc))
})

test_that("summarise_boot_mse gives no mse_se with one population and no components with one replicate each", {
  sq_err <- matrix(rexp(20), 2, 10)

  one_pop <- summarise_boot_mse(sq_err, rep(1L, 10), 1L)
  expect_true(all(is.na(one_pop$by_domain$mse_se)))
  expect_null(one_pop$var_comp)

  per_rep <- summarise_boot_mse(sq_err, 1:10, 10L)
  expect_true(all(!is.na(per_rep$by_domain$mse_se)))
  expect_null(per_rep$var_comp)
})

test_that("a failed first refit no longer makes every MSE NaN", {
  local_mocked_bindings(fit_zi = mock_fit_zi(fail_calls = 2))
  expect_warning(res <- fit_diag(B = 6L), "1 of 6 bootstrap refits failed")

  expect_true(all(!is.na(res$res$mse)))
  expect_true(all(res$res$n_boot_used == 5))
  expect_equal(res$boot_info$n_error, 1)
  expect_equal(res$boot_info$n_dropped, 1)
  expect_equal(res$boot_info$replicates$status[1], "error")
  expect_false(res$boot_info$replicates$used[1])
})

test_that("boot_failures = 'redraw' replaces failed refits", {
  local_mocked_bindings(fit_zi = mock_fit_zi(fail_calls = 2))
  expect_warning(res <- fit_diag(B = 6L, boot_failures = "redraw"), "0 still failed")

  expect_true(all(res$res$n_boot_used == 6))
  expect_equal(res$boot_info$n_dropped, 0)
  expect_equal(res$boot_info$n_pop_drawn, 2)
  expect_equal(res$boot_info$replicates$attempts[1], 2)
  expect_equal(res$boot_info$replicates$population[1], 2)
})

test_that("redrawing stops after the attempt limit", {
  local_mocked_bindings(fit_zi = mock_fit_zi(fail_calls = c(2, 8:20)))
  expect_warning(res <- fit_diag(B = 6L, boot_failures = "redraw"), "1 still failed")

  expect_equal(res$boot_info$replicates$attempts[1], 6)
  expect_equal(res$boot_info$n_dropped, 1)
  expect_true(all(res$res$n_boot_used == 5))
})

test_that("boot_warnings controls whether warned refits are used", {
  local_mocked_bindings(fit_zi = mock_fit_zi(warn_calls = c(3, 5)))
  expect_warning(kept <- fit_diag(B = 6L), "2 of 6 bootstrap refits gave warnings")
  expect_equal(kept$boot_info$n_warning, 2)
  expect_true(all(kept$res$n_boot_used == 6))

  local_mocked_bindings(fit_zi = mock_fit_zi(warn_calls = c(3, 5)))
  expect_warning(dropped <- fit_diag(B = 6L, boot_warnings = "fail"), "2 of 6 bootstrap refits failed")
  expect_true(all(dropped$res$n_boot_used == 4))
  expect_equal(dropped$boot_info$replicates$status[c(2, 4)], c("warning", "warning"))
})

test_that("invalid boot_warnings and boot_failures values are rejected", {
  expect_error(fit_diag(B = 4L, boot_warnings = "drop"), "Invalid boot_warnings")
  expect_error(fit_diag(B = 4L, boot_failures = "keep"), "Invalid boot_failures")
})

test_that("res gains n_boot_used and mse_se as its last columns", {
  res <- fit_diag(B = 6L, n_boot_pop = 3L)

  expect_equal(names(res$res), c("COUNTYFIPS", "mse", "est", "oos_flag", "n_boot_used", "mse_se"))
  expect_true(all(!is.na(res$res$mse_se)))
  expect_true(all(res$res$mse_se >= 0))
})

test_that("boot_info reports timing, costs, and the suggestion only when it can", {
  one <- fit_diag(B = 6L)
  expect_true(all(is.na(one$res$mse_se)))
  expect_null(one$boot_info$var_comp)
  expect_true(is.na(one$boot_info$suggested_n_boot_pop))
  expect_true(all(one$boot_info$time$phase >= 0))
  expect_true(all(one$boot_info$cost > 0))

  mid <- fit_diag(B = 6L, n_boot_pop = 3L)
  expect_s3_class(mid$boot_info$var_comp, "data.frame")
  expect_equal(nrow(mid$boot_info$var_comp), nrow(mid$res))
  expect_true(mid$boot_info$suggested_n_boot_pop %in% 1:6)

  expect_output(print(mid), "Bootstrap MSE")
})

test_that("the defaults leave the original columns unchanged", {
  set.seed(4)
  res <- suppressWarnings(suppressMessages(
    saeczi(samp, pop,
           lin_formula = DRYBIO_AG_TPA_live_ADJ ~ tcc16,
           domain_level = "COUNTYFIPS",
           mse_est = TRUE,
           B = 4L)
  ))
  expect_equal(names(res$res)[1:4], c("COUNTYFIPS", "mse", "est", "oos_flag"))
})

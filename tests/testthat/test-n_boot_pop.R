data(pop)
data(samp)

# a domain with auxiliary (pop) data but no sample data
oos_doms <- unique(pop$COUNTYFIPS)[1]
samp_oos <- samp[!(samp$COUNTYFIPS %in% oos_doms), ]

fit_k <- function(seed, ...) {
  set.seed(seed)
  suppressWarnings(suppressMessages(
    saeczi(samp_oos,
           pop,
           lin_formula = DRYBIO_AG_TPA_live_ADJ ~ tcc16,
           domain_level = "COUNTYFIPS",
           mse_est = TRUE,
           parallel = FALSE,
           ...)
  ))$res
}

test_that("boot_pop_index splits B replicates evenly across populations", {
  expect_equal(boot_pop_index(10L, 1L), rep(1L, 10))
  expect_equal(boot_pop_index(10L, 10L), 1:10)
  expect_equal(boot_pop_index(10L, 3L), c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L))
  expect_length(boot_pop_index(7L, 7L), 7)
})

test_that("n_boot_pop = 1L (the default) gives the same result as omitting it", {
  expect_identical(fit_k(11, B = 6L), fit_k(11, B = 6L, n_boot_pop = 1L))
})

test_that("n_boot_pop = B and uneven splits produce complete results", {
  res_b <- fit_k(12, B = 6L, n_boot_pop = 6L)
  res_4 <- fit_k(12, B = 6L, n_boot_pop = 4L)

  for (res in list(res_b, res_4)) {
    expect_equal(nrow(res), length(unique(pop$COUNTYFIPS)))
    expect_true(all(!is.na(res$mse)))
    expect_true(all(res$mse >= 0))
    expect_equal(res$COUNTYFIPS[res$oos_flag], oos_doms)
  }
})

test_that("n_boot_pop changes the MSE estimates but not the point estimates", {
  res_1 <- fit_k(13, B = 6L)
  res_b <- fit_k(13, B = 6L, n_boot_pop = 6L)

  expect_equal(res_1$est, res_b$est)
  expect_false(isTRUE(all.equal(res_1$mse, res_b$mse)))
})

test_that("invalid n_boot_pop values are rejected", {
  expect_error(fit_k(14, B = 6L, n_boot_pop = 0L), "n_boot_pop must be")
  expect_error(fit_k(14, B = 6L, n_boot_pop = 7L), "n_boot_pop must be")
  expect_error(fit_k(14, B = 6L, n_boot_pop = c(1L, 2L)), "n_boot_pop must be")
  expect_error(fit_k(14, B = 6L, n_boot_pop = NA_integer_), "n_boot_pop must be")
  expect_error(fit_k(14, B = 6L, n_boot_pop = 2), "class integer")
})

test_that("draw_boot_pop reproduces generate_boot_pop for the same seed", {
  out <- fit_zi(samp,
                DRYBIO_AG_TPA_live_ADJ ~ tcc16 + (1 | COUNTYFIPS),
                (DRYBIO_AG_TPA_live_ADJ != 0) ~ tcc16 + (1 | COUNTYFIPS),
                "COUNTYFIPS")
  setup <- boot_pop_setup(out, pop, "COUNTYFIPS", "tcc16", "tcc16")

  set.seed(15)
  pop_a <- generate_boot_pop(out, pop, "COUNTYFIPS", "tcc16", "tcc16")
  set.seed(15)
  pop_b <- draw_boot_pop(setup)

  expect_identical(pop_a, pop_b)
})

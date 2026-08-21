data(pop)
data(samp)

# a couple of domains that have auxiliary (pop) data but no sample data
oos_doms <- unique(pop$COUNTYFIPS)[1:2]
samp_oos <- samp[!(samp$COUNTYFIPS %in% oos_doms), ]

test_that("predict_oos = TRUE (default) produces synthetic estimates and mse for oos domains", {

  set.seed(42)
  suppressWarnings(suppressMessages(
    result <- saeczi(samp_oos,
                     pop,
                     lin_formula = DRYBIO_AG_TPA_live_ADJ ~ tcc16 + elev,
                     domain_level = "COUNTYFIPS",
                     mse_est = TRUE,
                     B = 10L,
                     parallel = FALSE)
  ))

  expect_equal(nrow(result$res), length(unique(pop$COUNTYFIPS)))
  expect_true(all(!is.na(result$res$est)))
  expect_true(all(!is.na(result$res$mse)))
  expect_contains(names(result$res), "oos_flag")
  expect_equal(sort(result$res$COUNTYFIPS[result$res$oos_flag]), sort(oos_doms))
  expect_equal(sum(result$res$oos_flag), length(oos_doms))

})

test_that("predict_oos = FALSE drops oos domains from the result", {

  set.seed(42)
  suppressWarnings(suppressMessages(
    result <- saeczi(samp_oos,
                     pop,
                     lin_formula = DRYBIO_AG_TPA_live_ADJ ~ tcc16 + elev,
                     domain_level = "COUNTYFIPS",
                     mse_est = TRUE,
                     B = 10L,
                     parallel = FALSE,
                     predict_oos = FALSE)
  ))

  expect_equal(nrow(result$res), length(unique(pop$COUNTYFIPS)) - length(oos_doms))
  expect_false(any(oos_doms %in% result$res$COUNTYFIPS))
  expect_true(all(!result$res$oos_flag))

})

test_that("oos_flag is present and all FALSE when there are no oos domains", {

  set.seed(42)
  suppressWarnings(
    result <- saeczi(samp,
                     pop,
                     lin_formula = DRYBIO_AG_TPA_live_ADJ ~ tcc16 + elev,
                     domain_level = "COUNTYFIPS",
                     mse_est = TRUE,
                     B = 10L,
                     parallel = FALSE)
  )

  expect_contains(names(result$res), "oos_flag")
  expect_true(all(!result$res$oos_flag))

})

test_that("mse_coefs includes the logistic model's domain random-effect variance", {

  out <- saeczi:::fit_zi(samp,
                        DRYBIO_AG_TPA_live_ADJ ~ tcc16 + (1 | COUNTYFIPS),
                        (DRYBIO_AG_TPA_live_ADJ != 0) ~ tcc16 + (1 | COUNTYFIPS),
                        "COUNTYFIPS")

  coefs <- saeczi:::mse_coefs(out$lmer, out$glmer)

  expect_true("sig2_b_hat" %in% names(coefs))
  expect_true(is.numeric(coefs$sig2_b_hat))
  expect_true(coefs$sig2_b_hat > 0)

})

test_that("samp_by_grp handles domains present in pop but absent from samp", {

  res <- saeczi:::samp_by_grp(samp_oos, pop, "COUNTYFIPS", 5)

  expect_length(res, 5)
  for (i in seq_along(res)) {
    # oos domains should contribute 0 rows to every bootstrap replicate
    expect_false(any(oos_doms %in% res[[i]]$COUNTYFIPS))
  }

})

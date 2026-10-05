data(samp)

res <- boot_rep(samp,
                "COUNTYFIPS",
                DRYBIO_AG_TPA_live_ADJ ~ tcc16 + (1 | COUNTYFIPS),
                DRYBIO_AG_TPA_live_ADJ != 0 ~ tcc16 + (1 | COUNTYFIPS))

names_nz <- samp[samp$DRYBIO_AG_TPA_live_ADJ > 0, ]$COUNTYFIPS |> unique()
names_all <- samp$COUNTYFIPS |> unique()


test_that("return object includes proper values", {
  expect_equal(names(res), c("beta_lm", "beta_glm", "u_lm", "u_glm",
                             "singular", "status", "message", "time"))
  expect_true(res$status %in% c("ok", "warning"))
  expect_true(is.numeric(res$time) && res$time >= 0)
})

test_that("a failed refit is reported and named after the logistic formula", {
  failed <- boot_rep(samp[0, ],
                     "COUNTYFIPS",
                     DRYBIO_AG_TPA_live_ADJ ~ tcc16 + elev + (1 | COUNTYFIPS),
                     DRYBIO_AG_TPA_live_ADJ != 0 ~ tcc16 + (1 | COUNTYFIPS))
  expect_equal(failed$status, "error")
  expect_true(nchar(failed$message) > 0)
  expect_named(failed$beta_lm, c("(Intercept)", "tcc16", "elev"))
  expect_named(failed$beta_glm, c("(Intercept)", "tcc16"))
  expect_true(all(is.na(failed$beta_glm)))
})


test_that("u_lm and u_glm are properly named vectors", {
  expect_named(res$u_lm, names_nz)
  expect_named(res$u_glm, names_all)
})

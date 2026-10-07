data(samp)
data(pop)

res <- samp_by_grp(samp, pop, "COUNTYFIPS", 10)

test_that("plots per group are correct", {
  out <- vector(mode = "logical", length = 10)
  chk <- as.data.frame(table(samp$COUNTYFIPS))
  for (i in 1:10) {
    cmp <- as.data.frame(table(res[[i]]$COUNTYFIPS))
    mtch <- all.equal(chk$Freq, cmp$Freq)
    out[i] <- mtch
  }
  expect_true(all(out))
})

test_that("each domain is sampled from its own rows when C and locale sort orders differ", {
  # "B" sorts before "a" in the C locale but after it in most other locales;
  # earlier versions mixed the two orders and sampled domains from each
  # other's rows. testthat runs tests in the C locale, so switch to a
  # typical locale to expose the difference.
  old_collate <- Sys.getlocale("LC_COLLATE")
  on.exit(Sys.setlocale("LC_COLLATE", old_collate), add = TRUE)
  locale_set <- suppressWarnings(Sys.setlocale("LC_COLLATE", "en_US.UTF-8")) != ""
  skip_if_not(locale_set, "en_US.UTF-8 locale not available")
  skip_if(identical(order(c("a", "B")), c(2L, 1L)), "locale sorts like C")

  pop_mixed <- data.frame(dom = rep(c("a", "B", "c"), times = c(5, 4, 6)),
                          x = rep(c(1, 2, 3), times = c(5, 4, 6)))
  samp_mixed <- data.frame(dom = rep(c("a", "B", "c"), times = c(1, 3, 2)))

  set.seed(1)
  res_mixed <- samp_by_grp(samp_mixed, pop_mixed, "dom", 20)

  for (i in seq_along(res_mixed)) {
    expect_equal(as.vector(table(res_mixed[[i]]$dom)[c("a", "B", "c")]), c(1, 3, 2))
    expect_equal(res_mixed[[i]]$x, c(a = 1, B = 2, c = 3)[res_mixed[[i]]$dom],
                 ignore_attr = TRUE)
  }
})

test_that("domains with no sample rows draw nothing", {
  pop_oos <- data.frame(dom = rep(c("d1", "d2", "d3"), each = 3), x = 1:9)
  samp_oos <- data.frame(dom = c("d1", "d3", "d3"))

  plan <- boot_samp_plan(samp_oos, pop_oos$dom, "dom")
  expect_equal(plan$n_samp, c(1, 0, 2))
  expect_equal(plan$n_pop, c(3, 3, 3))

  res_oos <- draw_boot_samples(plan, pop_oos, n_samp = 5)
  for (i in seq_along(res_oos)) {
    expect_false("d2" %in% res_oos[[i]]$dom)
    expect_equal(nrow(res_oos[[i]]), 3)
  }
})


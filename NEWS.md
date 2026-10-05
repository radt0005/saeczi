# saeczi (development version)

* `saeczi()` gains an `n_boot_pop` argument (default `1L`) for MSE
  estimation: the number of bootstrap populations, from 1 to `B`. The `B`
  bootstrap replicates are split as evenly as possible across `n_boot_pop`
  independently generated populations, and each replicate is compared with
  the true domain values of its own population; `B` stays the total number of
  replicates. The default reproduces earlier results exactly for the same seed
  (one population for all replicates, as in Chandra and Sud 2012). With one
  population, each domain's MSE depends on a single draw of its random effect,
  so it can vary greatly between runs however large `B` is. `n_boot_pop = B`
  draws a new population for every replicate (the standard parametric
  bootstrap), so the Monte Carlo error shrinks as `B` grows. Internally,
  `generate_boot_pop()` is split into `boot_pop_setup()` (computed once) and
  `draw_boot_pop()` (per population), and only one population is held in
  memory at a time.

* The parametric bootstrap MSE estimator now accounts for the variance
  contributed by the logistic model's domain random effect. Previously the
  bootstrap "truth" population plugged in the logistic model's estimated
  domain effects (BLUPs) directly (0 for any domain outside the sample),
  rather than drawing them from their estimated distribution the way the
  linear model's domain effect already was; this understated domain-level
  MSEs. (`mse_coefs()` / `generate_boot_pop()` in `R/utils.R`)

* `saeczi()` gains a `predict_oos` argument (default `TRUE`). When
  `pop_dat` contains domains with auxiliary data but no rows in
  `samp_dat`, `saeczi()` now automatically detects them and produces
  purely synthetic estimates for them (and MSE estimates, when
  `mse_est = TRUE`) using no domain-specific random effect, since none is
  estimable without sample data. Set `predict_oos = FALSE` to drop these
  domains from the result instead. The returned `res` data.frame gains a
  logical `oos_flag` column marking which domain estimates are synthetic.
  Fixed two bugs that had made MSE estimation error or return `NaN` for
  such domains: `samp_by_grp()` referenced a nonexistent literal `domain`
  column instead of the actual domain column, and bootstrap random-effect
  matrices (`u_lm`/`u_glm`) had no column at all (not just `NA` values)
  for domains that never appear in any bootstrap replicate.

# saeczi 0.2.0

# saeczi 0.1.1

# saeczi 0.1.0

* Initial CRAN submission.

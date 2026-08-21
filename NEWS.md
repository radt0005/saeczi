# saeczi (development version)

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

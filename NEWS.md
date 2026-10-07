# saeczi (development version)

* MSE estimation is much faster for large populations. Each bootstrap
  population now draws only its response vector instead of building a full
  copy of the population data, and bootstrap samples are drawn from domain
  row positions computed once per call instead of re-sorting the whole
  population for every sample. On a synthetic population of 4 million units
  in 250 domains (B = 10, `n_boot_pop = 10`), one call took 21 s instead of
  169 s. Results are identical for the same seed.

* Fixed: bootstrap samples could be drawn from the wrong domain when domain
  names sort differently in the C locale than in the system locale (for
  example, names that differ in letter case, such as "a" and "B"). Domains
  were ordered one way when counting and the other way when sorting rows,
  so each such domain was sampled from another domain's rows and got
  another domain's sample size. Domain names made only of digits were not
  affected.

* Bootstrap refits that fail or warn are now tracked instead of handled
  silently. Previously a refit that errored was dropped from the MSE without
  being counted, and a refit that gave a warning (e.g. a `glmer()`
  convergence warning) was used, with one console warning per refit. New
  arguments: `boot_warnings` (`"keep"`, the default, or `"fail"`) and
  `boot_failures` (`"drop"`, the default, or `"redraw"`, which replaces each
  failed replicate with a new population and sample, up to 5 attempts). One
  summary warning is given for warned refits, and another when more than 5%
  of refits fail. The defaults reproduce earlier MSE values exactly.

* Fixed: when the first bootstrap refit failed and the linear and logistic
  formulas had different predictors, every domain's MSE was `NaN`. A failed
  refit's logistic coefficients were named after the linear formula, which
  reordered the coefficient columns; they are now named correctly and matched
  to the design matrix by name.

* `res` gains two columns at the end when `mse_est = TRUE`: `n_boot_used`
  (replicates used for each domain) and `mse_se` (the Monte Carlo standard
  error of `mse`, computed from the same replicates with bootstrap
  populations as clusters; `NA` when `n_boot_pop = 1`).

* The result gains a `boot_info` element describing the bootstrap: refit
  counts by status, per-replicate status and messages, time per population
  and per phase, cost per population and per replicate, and, when
  `1 < n_boot_pop < B`, each domain's between- and within-population variance
  components with a suggested `n_boot_pop`. `print()` shows a short summary.

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

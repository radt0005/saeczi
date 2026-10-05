#' Fit a zero-inflation estimator.
#'
#' @param samp_dat A data.frame with domains, auxiliary variables, and the response variable of a sample
#' @param pop_dat A data.frame with domains and auxiliary variables of a population.
#' @param lin_formula Formula. Specification of the response and fixed effects of the linear regression model
#' @param log_formula Formula. Specification of the response and fixed effects of the logistic regression model
#' @param domain_level String. The column name in samp_dat and pop_dat that encodes the domain level
#' @param B Integer. The number of bootstraps to be used in MSE estimation.
#' @param mse_est Logical. Whether or not MSE estimation should happen.
#' @param estimand String. Whether the estimates should be 'totals' or 'means'.
#' @param parallel Logical. Should the MSE estimation be computed in parallel.
#' @param transform_fun Function. Function to be applied to the response variable prior to modeling.
#' @param inv_transform_fun Function. Inverse of transform_fun. Required if transform_fun is specified.
#' @param predict_oos Logical. If pop_dat contains domains that have no rows in samp_dat
#' (out-of-sample domains, i.e. domains with auxiliary data but no direct/sample data),
#' should estimates (and MSE estimates, if mse_est = TRUE) be produced for them anyway?
#' These are purely synthetic predictions: they use the fixed-effects portion of both
#' models only, with no domain-specific random effect, since none can be estimated
#' without sample data in the domain. Defaults to TRUE. If FALSE, such domains are
#' dropped from the returned `res` data.frame entirely.
#' @param n_boot_pop Integer. The number of bootstrap populations used for MSE
#' estimation, from 1 to `B`. The `B` bootstrap replicates are split as evenly as
#' possible across `n_boot_pop` independently generated bootstrap populations, and
#' each replicate is compared with the true domain values of its own population.
#' `B` stays the total number of replicates (model refits); for example, `B = 1000L`
#' and `n_boot_pop = 200L` give 200 populations with 5 replicates each.
#'
#' The default, `1L`, generates one population for all `B` replicates, as in the
#' published algorithm (Chandra and Sud 2012; White et al. 2024), and reproduces
#' results from earlier versions of saeczi for the same seed. With one population,
#' each domain's MSE estimate depends on a single random draw of that domain's
#' random effect, so it can vary greatly between runs however large `B` is.
#' `n_boot_pop = B` generates a new population for every replicate (the standard
#' parametric bootstrap; González-Manteiga et al. 2008), so the Monte Carlo error
#' of the MSE estimates shrinks as `B` grows. Values in between generate fewer
#' populations, which saves time when generating a population is expensive
#' compared with refitting the models. Only used when `mse_est = TRUE`.
#' @param boot_warnings String. How to treat bootstrap refits that finish with a
#' warning (for example, a convergence warning from `lme4::glmer()`): `"keep"`
#' (the default) uses them in the MSE, and `"fail"` treats them as failed
#' refits (see `boot_failures`). Either way they are counted in `boot_info`, and
#' one summary warning is given instead of one warning per refit.
#' @param boot_failures String. What to do with bootstrap refits that fail (an
#' error, or a warning when `boot_warnings = "fail"`): `"drop"` (the default)
#' leaves them out of the MSE, and `"redraw"` replaces each with a new bootstrap
#' population and sample, up to 5 attempts per replicate, so that all `B`
#' replicates are used when possible. Redrawing changes the random draws, so
#' results differ from `"drop"` whenever a refit fails. A warning is given when
#' more than 5\% of the first `B` refits fail.
#'
#' @returns
#' An object of class `zi_mod` with defined `print()` and `summary()` methods.
#' The object is structured like a list and contains the following elements:
#'
#' * call: The original function call
#'
#' * res: A data.frame containing the estimates and mse estimates. Includes a logical
#' `oos_flag` column that is TRUE for any domain whose estimate is purely synthetic,
#' i.e. domains present in pop_dat with no rows in samp_dat (see `predict_oos`), and
#' FALSE otherwise. When `mse_est = TRUE` it also has `n_boot_used`, the number
#' of bootstrap replicates used for that domain's MSE, and `mse_se`, the Monte
#' Carlo standard error of `mse`: how much `mse` would typically change if the
#' bootstrap were rerun with a different seed. `mse_se` is computed from the same
#' replicates, treating bootstrap populations as clusters, and is NA when
#' `n_boot_pop = 1`.
#'
#' * lin_mod: The modeling object used to fit the original linear model
#'
#' * log_mod: The modeling object used to fit the original logistic model
#'
#' * boot_info: Present when `mse_est = TRUE`. A list describing the bootstrap:
#' `B`, `n_boot_pop`, `n_pop_drawn` (including redraws), the `boot_warnings` and
#' `boot_failures` settings; counts of refits (`n_refits`, `n_ok`, `n_warning`,
#' `n_error`), of first-attempt failures (`n_failed_first`), of replicates left
#' out (`n_dropped`), and of singular fits (`n_singular`); `replicates`, a
#' data.frame with each replicate's population, final status, whether it was
#' used, number of attempts, refit time, and warning or error message; `time`,
#' the seconds spent per population and per phase (setup, populations,
#' sampling, refits, redraws, mse, total); `cost`, the wall-clock seconds per
#' population (drawing it and computing its true values) and per replicate
#' (sampling, refit, and prediction), excluding redraws; and, when
#' some populations have two or more replicates (`1 < n_boot_pop < B`),
#' `var_comp`, each domain's between-population (`sig2_pop`) and
#' within-population (`sig2_samp`) variance components of the squared errors
#' with the optimal replicates per population (`m_opt`), plus
#' `suggested_n_boot_pop`, the number of populations that would minimise the
#' Monte Carlo variance for the same run time (median over domains).
#'
#' @examples
#' data(pop)
#' data(samp)
#'
#' lin_formula <- DRYBIO_AG_TPA_live_ADJ ~ tcc16 + elev
#'
#' result <- saeczi(samp_dat = samp,
#'                  pop_dat = pop,
#'                  lin_formula = lin_formula,
#'                  log_formula = lin_formula,
#'                  domain_level = "COUNTYFIPS",
#'                  mse_est = FALSE)
#'
#' @export saeczi
#' @import stats
#' @importFrom rlang sym
#' @importFrom dplyr summarise group_by mutate left_join
#' @importFrom progressr progressor with_progress
#' @importFrom furrr future_map furrr_options future_map2
#' @importFrom purrr map map2 list_rbind

saeczi <- function(samp_dat,
                   pop_dat,
                   lin_formula,
                   log_formula = lin_formula,
                   domain_level,
                   B = 100L,
                   mse_est = FALSE,
                   estimand = "means",
                   parallel = FALSE,
                   transform_fun = NULL,
                   inv_transform_fun = NULL,
                   predict_oos = TRUE,
                   n_boot_pop = 1L,
                   boot_warnings = "keep",
                   boot_failures = "drop") {

  funcCall <- match.call()

  check_inherits("data.frame", samp_dat, pop_dat)
  check_inherits("formula", lin_formula, log_formula)
  check_inherits("character", domain_level, estimand, boot_warnings, boot_failures)
  check_inherits("integer", B, n_boot_pop)
  if (length(n_boot_pop) != 1 || is.na(n_boot_pop) || n_boot_pop < 1L || n_boot_pop > B) {
    stop(paste0("n_boot_pop must be a single integer from 1 to B (here 1 to ", B, ")."))
  }
  check_inherits("logical", mse_est, parallel, predict_oos)
  if (!is.null(transform_fun)) {
    if (is.null(inv_transform_fun)) {
      stop("inv_transform_fun must be specified when transform_fun is specified.")
    }
    check_inherits("function", transform_fun)
    check_inherits("function", inv_transform_fun)
    
    if (inv_transform_fun(transform_fun(14)) != 14) {
      warning("Rudimentary check on inv_transform_fun failed\nAre you sure inv_transform_fun is the inverse of transform_fun?")
    }
    
  }
  check_parallel(parallel)
  check_re(pop_dat, samp_dat, domain_level)

  if(!(estimand %in% c("means", "totals"))) {
    stop("Invalid estimand, must be either 'means' or 'totals'")
  }
  if (length(boot_warnings) != 1 || !(boot_warnings %in% c("keep", "fail"))) {
    stop("Invalid boot_warnings, must be either 'keep' or 'fail'")
  }
  if (length(boot_failures) != 1 || !(boot_failures %in% c("drop", "redraw"))) {
    stop("Invalid boot_failures, must be either 'drop' or 'redraw'")
  }

  # domains with auxiliary data in pop_dat but no rows in samp_dat: these
  # can only ever get a purely synthetic estimate (no domain-specific
  # random effect is estimable without sample data in the domain)
  oos_doms <- setdiff(unique(pop_dat[[domain_level]]), unique(samp_dat[[domain_level]]))

  if (length(oos_doms) > 0) {
    if (predict_oos) {
      message(sprintf(
        "%d domain(s) in pop_dat have no rows in samp_dat: %s.\nEstimates for these will be purely synthetic (no domain-specific random effect). See the `oos_flag` column of the result and the `predict_oos` argument.",
        length(oos_doms), paste(oos_doms, collapse = ", ")
      ))
    } else {
      message(sprintf(
        "%d domain(s) in pop_dat have no rows in samp_dat and will be dropped because predict_oos = FALSE: %s",
        length(oos_doms), paste(oos_doms, collapse = ", ")
      ))
      pop_dat <- pop_dat[pop_dat[[domain_level]] %in% unique(samp_dat[[domain_level]]), ]
      oos_doms <- character(0)
    }
  }

  Y <- toString(lin_formula[[2]])

  lin_X <- unlist(str_extract_all_base(deparse(lin_formula[[3]]), "\\w+"))
  log_X <- unlist(str_extract_all_base(deparse(log_formula[[3]]), "\\w+"))
  rand_intercept <- paste0("( 1 | ", domain_level, " )")
  lin_formula <- reformulate(c(lin_X, rand_intercept), response = Y)
  log_formula <- reformulate(c(log_X, rand_intercept), response = paste0(Y, "!= 0"))

  all_preds <- unique(c(lin_X, log_X))

  original_out <- fit_zi(samp_dat,
                         lin_formula,
                         log_formula,
                         domain_level,
                         transform_fun)

  mod1 <- original_out$lmer
  mod2 <- original_out$glmer

  .data <- pop_dat[, c(all_preds, domain_level)]

  original_pred <- collect_preds(mod1, mod2, estimand, .data, domain_level, inv_transform_fun)

  if (mse_est) {
    
    t_start <- proc.time()[["elapsed"]]
    
    pop_setup <- boot_pop_setup(original_out,
                                pop_dat,
                                domain_level,
                                log_X,
                                all_preds)

    boot_lin_formula <- reformulate(c(lin_X, rand_intercept), "response")
    boot_log_formula <- reformulate(c(log_X, rand_intercept), "response != 0")

    # Draws one bootstrap population, computes its truth, and draws n_samp
    # bootstrap samples from it. The population is then discarded, so only
    # one population is held in memory at a time.
    # Population time (draw and truth) and sampling time are recorded
    # separately, because sampling cost grows with the number of replicates,
    # not the number of populations.
    pop_times <- numeric(0)
    samp_time <- 0
    boot_truth <- list()
    draw_pop_and_samples <- function(n_samp) {
      t0 <- proc.time()[["elapsed"]]
      boot_pop_data <- draw_boot_pop(pop_setup)
      boot_truth[[length(boot_truth) + 1]] <<- compute_boot_truth(boot_pop_data,
                                                                  domain_level,
                                                                  estimand,
                                                                  inv_transform_fun)
      t1 <- proc.time()[["elapsed"]]
      samps <- samp_by_grp(samp_dat, boot_pop_data, domain_level, n_samp)
      pop_times <<- c(pop_times, t1 - t0)
      samp_time <<- samp_time + proc.time()[["elapsed"]] - t1
      samps
    }

    # With n_boot_pop = 1 the random draws happen in the same order as in
    # earlier versions.
    pop_id <- boot_pop_index(B, n_boot_pop)
    boot_samp_ls <- lapply(seq_len(n_boot_pop),
                           \(k) draw_pop_and_samples(sum(pop_id == k)))
    boot_samp_ls <- unlist(boot_samp_ls, recursive = FALSE)
    t_pops <- proc.time()[["elapsed"]]
    pop_time_first <- sum(pop_times)
    samp_time_first <- samp_time

    res <- fit_boot_reps(boot_samp_ls,
                         parallel,
                         domain_level,
                         boot_lin_formula,
                         boot_log_formula)

    is_failed <- \(r) r$status == "error" ||
      (boot_warnings == "fail" && r$status == "warning")
    statuses <- vapply(res, \(r) r$status, character(1))
    failed <- vapply(res, is_failed, logical(1))
    n_failed_first <- sum(failed)
    attempts <- rep(1L, B)

    # Optionally redraw failed replicates: each gets a new population and a
    # new sample (a new population, so its truth stays paired with it), up
    # to max_redraw attempts per replicate.
    max_redraw <- 5L
    redraw_round <- 0L
    t_redraw <- 0
    while (boot_failures == "redraw" && any(failed) && redraw_round < max_redraw) {
      redraw_round <- redraw_round + 1L
      t0 <- proc.time()[["elapsed"]]
      idx <- which(failed)
      new_samps <- lapply(idx, \(i) {
        pop_id[i] <<- length(boot_truth) + 1L
        draw_pop_and_samples(1L)[[1]]
      })
      new_res <- fit_boot_reps(new_samps,
                               parallel,
                               domain_level,
                               boot_lin_formula,
                               boot_log_formula)
      res[idx] <- new_res
      statuses <- c(statuses, vapply(new_res, \(r) r$status, character(1)))
      attempts[idx] <- attempts[idx] + 1L
      failed[idx] <- vapply(new_res, is_failed, logical(1))
      t_redraw <- t_redraw + proc.time()[["elapsed"]] - t0
    }
    t_fits <- proc.time()[["elapsed"]]

    params <- collect_boot_params(res)

    mse_out <- generate_mse(.data = pop_setup$pop_x,
                            truth = boot_truth,
                            pop_id = pop_id,
                            domain_level = domain_level,
                            beta_lm_mat = params$beta_lm_mat,
                            beta_glm_mat = params$beta_glm_mat,
                            u_lm = params$u_lm,
                            u_glm = params$u_glm,
                            lin_X = lin_X,
                            log_X = log_X,
                            estimand = estimand,
                            inv = inv_transform_fun,
                            failed = failed)

    mse_summ <- summarise_boot_mse(mse_out$sq_err, pop_id, n_boot_pop)
    t_end <- proc.time()[["elapsed"]]

    mse_df <- data.frame(mse_out$domain, mse_summ$by_domain)
    names(mse_df)[1] <- domain_level

    # Cost per population (draw and truth) and per replicate (sampling,
    # refit, and prediction), as wall-clock time, so a parallel run's refit
    # cost reflects its workers. Redraws are left out of both.
    n_fits <- length(statuses)
    t_refits <- t_fits - t_pops - t_redraw
    c_pop <- pop_time_first / n_boot_pop
    c_fit <- (samp_time_first + t_refits + (t_end - t_fits)) / B
    suggestion <- suggest_n_boot_pop(mse_summ$var_comp, c_pop, c_fit, B)

    var_comp <- mse_summ$var_comp
    if (!is.null(var_comp)) {
      var_comp <- data.frame(mse_out$domain, var_comp, m_opt = suggestion$m_opt)
      names(var_comp)[1] <- domain_level
    }

    boot_info <- list(
      B = B,
      n_boot_pop = n_boot_pop,
      n_pop_drawn = length(boot_truth),
      boot_warnings = boot_warnings,
      boot_failures = boot_failures,
      n_refits = n_fits,
      n_ok = sum(statuses == "ok"),
      n_warning = sum(statuses == "warning"),
      n_error = sum(statuses == "error"),
      n_failed_first = n_failed_first,
      n_dropped = sum(failed),
      n_singular = sum(vapply(res, \(r) isTRUE(r$singular), logical(1))),
      replicates = data.frame(
        replicate = seq_len(B),
        population = pop_id,
        status = vapply(res, \(r) r$status, character(1)),
        used = colSums(!is.na(mse_out$sq_err)) > 0,
        attempts = attempts,
        time = vapply(res, \(r) r$time, numeric(1)),
        message = vapply(res, \(r) r$message, character(1))
      ),
      time = list(
        population = pop_times,
        phase = c(setup = t_pops - t_start - pop_time_first - samp_time_first,
                  populations = pop_time_first,
                  sampling = samp_time_first,
                  refits = t_refits,
                  redraws = t_redraw,
                  mse = t_end - t_fits,
                  total = t_end - t_start)
      ),
      cost = c(population = c_pop, replicate = c_fit),
      var_comp = var_comp,
      suggested_n_boot_pop = suggestion$n_boot_pop
    )

    if (n_failed_first / B > 0.05) {
      warning(sprintf(
        "%d of %d bootstrap refits failed (%.0f%%)%s. See `boot_info$replicates`.",
        n_failed_first, B, 100 * n_failed_first / B,
        if (boot_failures == "redraw") {
          sprintf("; %d still failed after up to %d redraws and were left out of the MSE", sum(failed), max_redraw)
        } else {
          " and were left out of the MSE"
        }
      ))
    }

    if (boot_info$n_warning > 0 && boot_warnings == "keep") {
      warning(sprintf(
        "%d of %d bootstrap refits gave warnings (e.g., convergence) and were kept. See `boot_info$replicates`; use boot_warnings = \"fail\" to leave them out.",
        boot_info$n_warning, n_fits
      ))
    }

    final_df <- mse_df |>
      left_join(original_pred, by = domain_level)

  } else {

    final_df <- original_pred

  }

  oos_flag_df <- data.frame(unique(pop_dat[[domain_level]]))
  names(oos_flag_df) <- domain_level
  oos_flag_df$oos_flag <- oos_flag_df[[domain_level]] %in% oos_doms

  final_df <- final_df |>
    left_join(oos_flag_df, by = domain_level)

  # keep the original column order and add the new MSE columns at the end
  if (mse_est) {
    new_cols <- c("n_boot_used", "mse_se")
    final_df <- final_df[ , c(setdiff(names(final_df), new_cols), new_cols)]
  }

  out <- list(
    call = funcCall,
    res = final_df,
    lin_mod = original_out$lmer,
    log_mod = original_out$glmer
  )

  if (mse_est) {
    out$boot_info <- boot_info
  }

  structure(out, class = "zi_mod")

}

#' @export
print.zi_mod <- function(x, ...) {
  cat("\nCall:\n")
  cat(deparse(x$call))
  cat("\n\n")

  cat("Linear Model: \n")
  cat("- Fixed effects: \n")
  print(summary(x$lin_mod)$coefficients[ ,1])
  cat("\n")
  cat("- Random effects: \n")
  print(summary(x$lin_mod)$varcor)
  cat("\n")

  cat("Logistic Model: \n")
  cat("- Fixed effects: \n")
  print(summary(x$log_mod)$coefficients[ ,1])
  cat("\n")
  cat("- Random effects: \n")
  print(summary(x$log_mod)$varcor)
  cat("\n")

  if (!is.null(x$boot_info)) {
    bi <- x$boot_info
    cat("Bootstrap MSE: \n")
    cat(sprintf("- %d replicates from %d population(s); %d population(s) drawn in total\n",
                bi$B, bi$n_boot_pop, bi$n_pop_drawn))
    cat(sprintf("- Refits: %d ok, %d with warnings (%s), %d errors; %d replicate(s) left out\n",
                bi$n_ok, bi$n_warning, bi$boot_warnings, bi$n_error, bi$n_dropped))
    cat(sprintf("- Time: %.1f s total; %.3f s per population, %.3f s per replicate\n",
                bi$time$phase[["total"]], bi$cost[["population"]], bi$cost[["replicate"]]))
    if (!is.na(bi$suggested_n_boot_pop)) {
      cat(sprintf("- Suggested n_boot_pop for B = %d: %d\n", bi$B, bi$suggested_n_boot_pop))
    }
    cat("\n")
  }
}

#' @export
summary.zi_mod <- function(object, ...) {
  out <- list(
    lin_mod = summary(object$lin_mod),
    log_mod = summary(object$log_mod)
  )

  class(out) <- "summary.zinf_bayes"
  out
}

#' @export
print.summary.zi_mod <- function(x, ...) {
  print(x$lin_mod)
  cat("\n")
  print(x$log_mod)
}

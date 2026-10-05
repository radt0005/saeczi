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
#' FALSE otherwise.
#'
#' * lin_mod: The modeling object used to fit the original linear model
#'
#' * log_mod: The modeling object used to fit the original logistic model
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
                   n_boot_pop = 1L) {

  funcCall <- match.call()

  check_inherits("data.frame", samp_dat, pop_dat)
  check_inherits("formula", lin_formula, log_formula)
  check_inherits("character", domain_level, estimand)
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
    
    pop_setup <- boot_pop_setup(original_out,
                                pop_dat,
                                domain_level,
                                log_X,
                                all_preds)

    boot_lin_formula <- reformulate(c(lin_X, rand_intercept), "response")
    boot_log_formula <- reformulate(c(log_X, rand_intercept), "response != 0")

    # Each of the n_boot_pop bootstrap populations is generated, used for its
    # truth and its share of the B bootstrap samples, and then discarded, so
    # only one population is held in memory at a time. With n_boot_pop = 1 the
    # random draws happen in the same order as in earlier versions.
    pop_id <- boot_pop_index(B, n_boot_pop)
    boot_truth <- vector("list", length = n_boot_pop)
    boot_samp_ls <- vector("list", length = n_boot_pop)

    for (k in seq_len(n_boot_pop)) {
      boot_pop_data <- draw_boot_pop(pop_setup)
      boot_truth[[k]] <- compute_boot_truth(boot_pop_data,
                                            domain_level,
                                            estimand,
                                            inv_transform_fun)
      boot_samp_ls[[k]] <- samp_by_grp(samp_dat, boot_pop_data, domain_level, sum(pop_id == k))
    }

    rm(boot_pop_data)
    boot_samp_ls <- unlist(boot_samp_ls, recursive = FALSE)

    if (parallel) {
      with_progress({
        boot_res <- boot_rep_par(x = 1:B,
                                 boot_lst = boot_samp_ls,
                                 domain_level,
                                 boot_lin_formula,
                                 boot_log_formula,
                                 pop_setup$pop_x,
                                 boot_truth,
                                 pop_id,
                                 estimand,
                                 lin_X,
                                 log_X,
                                 inv_transform_fun)
        })

    } else {

      res <-
        map(.x = boot_samp_ls,
            .f = \(.x) {
              boot_rep(boot_samp = .x,
                       domain_level,
                       boot_lin_formula,
                       boot_log_formula)
            },
            .progress = list(
              type = "iterator",
              clear = TRUE
            ))

      beta_lm_mat <- res |>
        map(.f = ~ as.data.frame(t(.x$beta_lm))) |>
        list_rbind() |> 
        as.matrix()

      beta_glm_mat <- res |>
        map(.f = ~ as.data.frame(t(.x$beta_glm))) |>
        list_rbind() |> 
        as.matrix()

      u_lm <- res |>
        map(.f = ~ as.data.frame(t(.x$u_lm))) |>
        list_rbind() |> 
        as.matrix()

      u_glm <- res |>
        map(.f = ~ as.data.frame(t(.x$u_glm))) |>
        list_rbind() |> 
        as.matrix()

      # see comments in boot_rep_par() (utils.R) for why both are zeroed
      u_lm[is.na(u_lm)] <- 0
      u_glm[is.na(u_glm)] <- 0

      preds_full <- generate_mse(.data = pop_setup$pop_x,
                                 truth = boot_truth,
                                 pop_id = pop_id,
                                 domain_level = domain_level,
                                 beta_lm_mat = beta_lm_mat,
                                 beta_glm_mat = beta_glm_mat,
                                 u_lm = u_lm,
                                 u_glm = u_glm,
                                 lin_X = lin_X,
                                 log_X = log_X,
                                 estimand = estimand,
                                 inv = inv_transform_fun)
    

      boot_res <- preds_full

    }

    mse_df <- setNames(boot_res,
                       c(domain_level, "mse"))

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

  out <- list(
    call = funcCall,
    res = final_df,
    lin_mod = original_out$lmer,
    log_mod = original_out$glmer
  )

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

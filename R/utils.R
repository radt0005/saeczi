#' Generates B-many bootstrap samples
#' 
#' Ensures that each group in a bootstrap sample has the same number of rows as 
#' that group did in the original sample data. 
#' 
#' @param samp The sample data
#' @param pop The bootstrap population data
#' @param dom_nm Character string of the domain identifier name as it appears in samp and pop
#' @param B Integer. The number of bootstrap samples to produce
#' 
#' @return A list of length B where each item is a unique bootstrap sample
#' @noRd
#' 
samp_by_grp <- function(samp, pop, dom_nm, B) {
  
  plan <- boot_samp_plan(samp, pop[[dom_nm]], dom_nm)
  
  draw_boot_samples(plan, pop, n_samp = B)
  
}

#' Plan the stratified bootstrap sampling
#' 
#' Finds, once, the population rows of each domain and the number of rows to
#' draw from each (the domain's sample size; 0 for domains with no sample
#' rows). Domains are kept in the order `dplyr::count()` gives, so the random
#' draws happen in the same order as in earlier versions.
#' 
#' @param samp The sample data
#' @param pop_dom The domain of each population row
#' @param dom_nm Character string of the domain identifier name
#' 
#' @return A list with `rows` (the population row numbers of each domain),
#' `n_pop` (domain population sizes), and `n_samp` (rows to draw per domain)
#' @noRd
boot_samp_plan <- function(samp, pop_dom, dom_nm) {
  
  num_plots <- dplyr::count(samp, !!rlang::sym(dom_nm))
  
  pop_doms <- setNames(data.frame(pop_dom), dom_nm)
  pop_counts <- dplyr::count(pop_doms, !!rlang::sym(dom_nm))
  doms <- pop_counts[[dom_nm]]
  
  # domains with no rows in samp (out-of-sample domains) draw 0 rows
  n_samp <- num_plots$n[match(doms, num_plots[[dom_nm]])]
  n_samp[is.na(n_samp)] <- 0L
  
  rows <- split(seq_along(pop_dom), factor(pop_dom, levels = doms))
  
  list(rows = unname(rows),
       n_pop = pop_counts$n,
       n_samp = n_samp)
  
}

#' Draw stratified bootstrap samples from a population
#' 
#' Each sample draws, with replacement, `plan$n_samp` rows from each domain's
#' population rows. Only the sampled rows are copied, so the cost does not
#' grow with the population size.
#' 
#' @param plan The list returned by `boot_samp_plan()`
#' @param pop_x A data.frame of population rows (domain and predictor columns)
#' @param response Optional numeric vector of population responses, added as a
#' `response` column; NULL if `pop_x` already contains it
#' @param n_samp Integer. The number of samples to draw
#' 
#' @return A list of `n_samp` data.frames
#' @noRd
draw_boot_samples <- function(plan, pop_x, response = NULL, n_samp) {
  
  lapply(seq_len(n_samp), \(i) {
    
    ids <- unlist(lapply(seq_along(plan$rows), \(d) {
      plan$rows[[d]][sample.int(plan$n_pop[d], plan$n_samp[d], replace = TRUE)]
    }))
    
    out <- pop_x[ids, , drop = FALSE]
    if (!is.null(response)) {
      out$response <- response[ids]
    }
    out
    
  })
  
}


#' Fits the two models in a zi-estimator
#' 
#' @param samp_dat Sample data
#' @param lin_formula Formula to be used in the linear regression model
#' @param log_formula Formula to be used in the logistic regression model
#' @param domain_level Character. Domain identifier name.
#' 
#' @return A list containing the two model objects
#' @noRd
#' 
fit_zi <- function(samp_dat,
                   lin_formula,
                   log_formula,
                   domain_level,
                   transform_fun = NULL) {
  
  Y <- deparse(lin_formula[[2]])
  
  if (!is.null(transform_fun)) {
    samp_dat[[Y]] <- transform_fun(samp_dat[[Y]])
  } 
  
  # creating nonzero version of our sample data set
  nz <- samp_dat[samp_dat[[Y]] > 0, ]
  
  # fit linear mixed model on nonzero data
  lmer_nz <- suppressMessages(
    lme4::lmer(lin_formula, data = nz)
  )
  
  # Fit logistic mixed effects on ALL data
  glmer_z <- suppressMessages(
    lme4::glmer(log_formula, data = samp_dat, family = 'binomial')
  )
  
  return(list(lmer = lmer_nz, glmer = glmer_z))
  
}

#' Generates the bootstrap squared errors
#' 
#' @param .data The population data: the domain column and the predictor columns
#' @param truth A list with one data.frame per bootstrap population, each containing the true domain level values from that population
#' @param pop_id An integer vector of length B giving the bootstrap population (index into `truth`) of each replicate
#' @param domain_level Character. Domain identifier name
#' @param beta_lm_mat A matrix containing the fixed effects coefficients resulting from fitting the linear model to each bootstrap sample
#' @param beta_glm_mat A matrix containing the fixed effects coefficients resulting from fitting the logistic model to each bootstrap sample
#' @param u_lm A matrix containing the random effects values for each domain resulting from fitting the linear model to each bootstrap sample
#' @param u_glm A matrix containing the random effects values for each domain resulting from fitting the linear model to each bootstrap sample
#' @param lin_X A character vector with the names of the predictor variables used in the linear model
#' @param log_X A character vector with the names of the predictor variables used in the logistic model
#' @param estimand A string specifying whether the estimates should be 'totals' or 'means'
#' @param failed A logical vector of length B; replicates marked TRUE are excluded (their squared errors are set to NA)
#' 
#' @return A list with `domain` (the domain identifiers, in the original
#' type) and `sq_err`, a domain-by-replicate matrix of squared errors
#' (NA for excluded or unusable replicates)
#' @noRd
generate_mse <- function(.data,
                         truth,
                         pop_id,
                         domain_level,
                         beta_lm_mat,
                         beta_glm_mat,
                         u_lm,
                         u_glm,
                         lin_X,
                         log_X,
                         estimand,
                         inv,
                         failed = rep(FALSE, length(pop_id))) {
  
  boot_pop_by_dom <- split(.data, f = .data[[domain_level]])
  
  design_mat_ls <-  boot_pop_by_dom |> 
    map(.f = function(.x) {
      dmat_lm <- model.matrix(~., .x[ ,lin_X, drop = FALSE])
      dmat_glm <- model.matrix(~., .x[ ,log_X, drop = FALSE])
      return(list(design_mat_lm = dmat_lm,
                  design_mat_glm = dmat_glm))
    })
  
  n_doms <- length(design_mat_ls)
  dom_order <- names(boot_pop_by_dom)

  # Align each random-effect matrix to have exactly one column per domain
  # in dom_order, in that order (generate_preds()/preds_calc() indexes
  # u_lm/u_glm by column position j, matched positionally to
  # design_mats[[j]]). A domain can be entirely absent from u_lm (no
  # replicate ever had a positive-response row for it) or from u_glm (no
  # replicate ever had any row for it at all -- always true of
  # out-of-sample domains, which contribute 0 rows to every bootstrap
  # replicate); list_rbind() only unions the columns that did appear
  # across replicates, so such a domain has no column yet, not just NA
  # cells. Filling it with 0 gives a purely synthetic (population-average)
  # domain effect for every replicate, consistent with how NA cells within
  # an existing column are already zeroed before this function is called.
  align_to_doms <- function(u_mat, dom_order) {
    missing <- setdiff(dom_order, colnames(u_mat))
    if (length(missing) > 0) {
      to_append <- matrix(0, nrow = nrow(u_mat), ncol = length(missing),
                          dimnames = list(NULL, missing))
      u_mat <- cbind(u_mat, to_append)
    }
    u_mat[, dom_order, drop = FALSE]
  }

  u_lm <- align_to_doms(u_lm, dom_order)
  u_glm <- align_to_doms(u_glm, dom_order)

  # Align the fixed-effect columns by name to the design-matrix columns
  # (generate_preds() multiplies them positionally). A failed replicate's
  # coefficients are all NA, and list_rbind() orders columns by first
  # appearance, so without this a failed first replicate could reorder the
  # columns. A coefficient missing from every replicate (e.g. dropped as
  # rank deficient) becomes an NA column, so those replicates give NA
  # predictions and are left out, as before.
  align_to_design <- function(beta_mat, design_cols) {
    missing <- setdiff(design_cols, colnames(beta_mat))
    if (length(missing) > 0) {
      to_append <- matrix(NA_real_, nrow = nrow(beta_mat), ncol = length(missing),
                          dimnames = list(NULL, missing))
      beta_mat <- cbind(beta_mat, to_append)
    }
    beta_mat[, design_cols, drop = FALSE]
  }

  beta_lm_mat <- align_to_design(beta_lm_mat, colnames(design_mat_ls[[1]]$design_mat_lm))
  beta_glm_mat <- align_to_design(beta_glm_mat, colnames(design_mat_ls[[1]]$design_mat_glm))

  dom_res_wide <- generate_preds(beta_lm = beta_lm_mat,
                                 beta_glm = beta_glm_mat,
                                 u_lm = u_lm,
                                 u_glm = u_glm,
                                 design_mats = design_mat_ls,
                                 J = n_doms,
                                 estimand = estimand,
                                 inv = inv)
  
  truth_ordered <- truth[[1]][order(match(truth[[1]][[domain_level]], dom_order)), ]
  
  # one column of true domain values per bootstrap population, then one
  # column per replicate, so each replicate is compared with the truth of
  # the population it was sampled from
  truth_mat <- matrix(
    vapply(truth,
           \(.t) .t$domain_est[match(truth_ordered[[domain_level]], .t[[domain_level]])],
           numeric(nrow(truth_ordered))),
    nrow = nrow(truth_ordered)
  )
  
  sq_err <- (dom_res_wide - truth_mat[ , pop_id, drop = FALSE])^2
  sq_err[ , failed] <- NA
  
  list(domain = truth_ordered[[domain_level]],
       sq_err = sq_err)
  
}

#' Summarise the bootstrap squared errors by domain
#' 
#' `mse` is the mean of each domain's squared errors over the replicates
#' used. `mse_se` is its Monte Carlo standard error, treating bootstrap
#' populations as clusters: with S_k and C_k the sum and count of a domain's
#' squared errors in population k, and K the number of populations used,
#' SE^2 = K / (K - 1) * sum_k (S_k - mse * C_k)^2 / (sum_k C_k)^2. With
#' equal-size populations this equals sd(population means) / sqrt(K). It is
#' NA when fewer than two populations were requested or used.
#' 
#' When some populations have two or more replicates, the squared errors
#' are also split into between-population (`sig2_pop`) and within-population
#' (`sig2_samp`) variance components by one-way ANOVA for unbalanced groups.
#' 
#' @param sq_err A domain-by-replicate matrix of squared errors
#' @param pop_id An integer vector giving each replicate's population
#' @param n_boot_pop Integer. The number of populations requested
#' 
#' @return A list with `by_domain` (a data.frame with `mse`, `n_boot_used`
#' and `mse_se`) and `var_comp` (a data.frame with `sig2_pop` and
#' `sig2_samp`, or NULL when the components cannot be separated)
#' @noRd
summarise_boot_mse <- function(sq_err, pop_id, n_boot_pop) {
  
  used <- !is.na(sq_err)
  n_used <- rowSums(used)
  mse <- rowMeans(sq_err, na.rm = TRUE)
  
  sq_err0 <- sq_err
  sq_err0[!used] <- 0
  
  # population-by-domain sums, counts, and sums of squares
  S <- rowsum(t(sq_err0), pop_id, reorder = FALSE)
  C <- rowsum(t(used) * 1, pop_id, reorder = FALSE)
  Q <- rowsum(t(sq_err0^2), pop_id, reorder = FALSE)
  
  K_used <- colSums(C > 0)
  
  mse_se <- rep(NA_real_, length(mse))
  if (n_boot_pop >= 2) {
    resid <- S - sweep(C, 2, mse, `*`)
    se2 <- K_used / (K_used - 1) * colSums(resid^2) / n_used^2
    mse_se <- ifelse(K_used >= 2, sqrt(se2), NA_real_)
  }
  
  var_comp <- NULL
  if (n_boot_pop >= 2 && any(C >= 2)) {
    N <- colSums(C)
    ss_within <- colSums(Q) - colSums(ifelse(C > 0, S^2 / C, 0))
    ss_between <- colSums(ifelse(C > 0, S^2 / C, 0)) - N * mse^2
    ms_within <- ss_within / (N - K_used)
    ms_between <- ss_between / (K_used - 1)
    n0 <- (N - colSums(C^2) / N) / (K_used - 1)
    ok <- K_used >= 2 & N > K_used
    var_comp <- data.frame(
      sig2_pop = ifelse(ok, (ms_between - ms_within) / n0, NA_real_),
      sig2_samp = ifelse(ok, ms_within, NA_real_)
    )
  }
  
  list(by_domain = data.frame(mse = mse,
                              n_boot_used = n_used,
                              mse_se = mse_se),
       var_comp = var_comp)
  
}

#' Suggest a number of bootstrap populations
#' 
#' Uses the two-stage sampling optimum m_opt = sqrt((c_pop / c_fit) *
#' (sig2_samp / sig2_pop)) replicates per population, computed per domain
#' and summarised by the median, then n_boot_pop = B / m_opt, rounded and
#' kept within 1 to B.
#' 
#' @param var_comp The `var_comp` data.frame from `summarise_boot_mse()`
#' @param c_pop Seconds per bootstrap population
#' @param c_fit Seconds per replicate (refit and prediction)
#' @param B Integer. Total number of bootstrap replicates
#' 
#' @return A list with the per-domain `m_opt` and the suggested `n_boot_pop`
#' (NA when it cannot be computed)
#' @noRd
suggest_n_boot_pop <- function(var_comp, c_pop, c_fit, B) {
  
  if (is.null(var_comp) || !is.finite(c_pop) || !is.finite(c_fit) || c_fit <= 0) {
    return(list(m_opt = NULL, n_boot_pop = NA_integer_))
  }
  
  ratio <- var_comp$sig2_samp / var_comp$sig2_pop
  # no detectable between-population variance: fewer populations suffice
  ratio[!is.na(var_comp$sig2_pop) & var_comp$sig2_pop <= 0] <- Inf
  m_opt <- sqrt((c_pop / c_fit) * ratio)
  
  m_med <- stats::median(m_opt, na.rm = TRUE)
  if (is.na(m_med)) {
    return(list(m_opt = m_opt, n_boot_pop = NA_integer_))
  }
  
  k <- if (is.infinite(m_med)) 1 else round(B / max(m_med, 1))
  
  list(m_opt = m_opt,
       n_boot_pop = as.integer(min(max(k, 1), B)))
  
}

#' Generate the Bootstrap population data
#' 
#' Wrapper that generates a single bootstrap population. `saeczi()` calls
#' `boot_pop_setup()` once and `draw_boot_pop()` once per bootstrap
#' population instead.
#' 
#' @param original_out List containing original model objects
#' @param pop_dat The population data frame
#' @param domain_level Character. The domain column names in pop_dat
#' @param log_X Vector of characters containing logistic model predictor names
#' @param all_preds Vector of characters containing all predictor names
#' 
#' @return The population bootstrap data
#' @noRd

generate_boot_pop <- function(original_out, 
                              pop_dat,
                              domain_level,
                              log_X,
                              all_preds,
                              transform_fun) {
  
  setup <- boot_pop_setup(original_out, pop_dat, domain_level, log_X, all_preds)
  
  draw_boot_pop(setup)
  
}

#' Precompute the parts of a bootstrap population that do not change
#' between bootstrap populations
#' 
#' The design matrix, the fixed-effect parts of both models, and the
#' variance-component estimates are the same for every bootstrap population;
#' only the random draws in `draw_boot_pop()` differ.
#' 
#' @inheritParams generate_boot_pop
#' 
#' @return A list used by `draw_boot_pop()`
#' @noRd
boot_pop_setup <- function(original_out,
                           pop_dat,
                           domain_level,
                           log_X,
                           all_preds) {
  
  zi_mod_coefs <- mse_coefs(original_out$lmer, original_out$glmer)

  x_matrix <- model.matrix(
    as.formula(paste0(" ~ ", paste(all_preds, collapse = " + "))),
    data = pop_dat[ , all_preds, drop = FALSE]
  )

  # tweak for allowing new levels (also covers domains present in pop_dat
  # but absent from the original sample, i.e. out-of-sample domains)
  pop_doms <- unique(pop_dat[[domain_level]])
  all_doms <- union(pop_doms, zi_mod_coefs$domain_levels)

  dom <- pop_dat[ , domain_level, drop = TRUE]

  list(
    coefs = zi_mod_coefs,
    pop_x = pop_dat[ , c(domain_level, all_preds)],
    dom = dom,
    dom_idx = match(dom, all_doms),
    all_doms = all_doms,
    log_fixed = as.vector(x_matrix[ , c("(Intercept)", log_X)] %*% zi_mod_coefs$alpha_1),
    lin_fixed = as.vector(x_matrix[, colnames(model.matrix(original_out$lmer))] %*% zi_mod_coefs$beta_hat)
  )
  
}

#' Draw one bootstrap population
#' 
#' @param setup The list returned by `boot_pop_setup()`
#' 
#' @return The population bootstrap data: the population's domain and
#' predictor columns plus the drawn `response`
#' @noRd
draw_boot_pop <- function(setup) {
  
  data.frame(setup$pop_x, response = draw_boot_response(setup))
  
}

#' Draw the response of one bootstrap population
#' 
#' The random draws are made in a fixed order (unit errors, linear domain
#' effects, logistic domain effects, nonzero indicators), so a population
#' drawn after a given seed matches earlier versions of saeczi. Only the
#' response vector is built; the bootstrap does not need a full copy of the
#' population's predictors for each population.
#' 
#' @param setup The list returned by `boot_pop_setup()`
#' 
#' @return A numeric vector with one response per population row
#' @noRd
draw_boot_response <- function(setup) {
  
  zi_mod_coefs <- setup$coefs
  
  eps_ij <- rnorm(length(setup$dom_idx), 0, sqrt(zi_mod_coefs$sig2_eps_hat))

  area_re <- rnorm(length(setup$all_doms), 0, sqrt(zi_mod_coefs$sig2_mu_hat))

  # Draw a fresh domain random effect for the logistic model from its
  # estimated distribution N(0, sig2_b_hat), the same way u_j is drawn
  # above for the linear model, rather than plugging in the original
  # model's point-estimate BLUPs (which had no domain-level variance and
  # were fixed at 0 for any domain outside the original sample). This is
  # what lets the bootstrap "truth" population reflect the variability the
  # logistic model actually contributes to the domain-level estimates,
  # for in-sample and out-of-sample domains alike.
  b_i <- rnorm(length(setup$all_doms), 0, sqrt(zi_mod_coefs$sig2_b_hat))

  p_hat_i <- binomial()$linkinv(setup$log_fixed + b_i[setup$dom_idx])
  
  delta_i_star <- rbinom(length(p_hat_i), 1, p_hat_i)
  
  (setup$lin_fixed + area_re[setup$dom_idx] + eps_ij) * delta_i_star
  
}

#' Compute the true domain values of a bootstrap population
#' 
#' @param dom The domain of each population row
#' @param response The bootstrap population's response
#' @param domain_level Character. Domain identifier name
#' @param estimand A string specifying whether the estimates should be 'totals' or 'means'
#' @param inv_transform_fun Function or NULL. Inverse of the response transformation
#' 
#' @return A data.frame with one row per domain and a `domain_est` column
#' @noRd
compute_boot_truth <- function(dom,
                               response,
                               domain_level,
                               estimand,
                               inv_transform_fun) {
  
  if (!is.null(inv_transform_fun)) {
    response <- inv_transform_fun(response)
  }
  
  boot_pop_data <- setNames(data.frame(dom, response), c(domain_level, "response"))
  
  agg_fun <- if (estimand == "means") mean else sum
  
  boot_pop_data |>
    group_by(!!sym(domain_level)) |>
    summarise(domain_est = agg_fun(response))
  
}

#' Assign bootstrap replicates to bootstrap populations
#' 
#' Splits B replicates as evenly as possible across n_boot_pop populations;
#' when B is not divisible by n_boot_pop, the first B %% n_boot_pop
#' populations get one extra replicate.
#' 
#' @param B Integer. Total number of bootstrap replicates
#' @param n_boot_pop Integer. Number of bootstrap populations
#' 
#' @return An integer vector of length B giving each replicate's population
#' @noRd
boot_pop_index <- function(B, n_boot_pop) {
  
  sizes <- rep(B %/% n_boot_pop, n_boot_pop) +
    (seq_len(n_boot_pop) <= B %% n_boot_pop)
  
  rep(seq_len(n_boot_pop), times = sizes)
  
}

#' Bootstrap refits for the parallel option
#' 
#' @param x The vector of replicate indexes
#' @param boot_lst A list where each element contains a bootstrap sample
#' @param domain_level Character. Domain identifier name
#' @param boot_lin_formula The formula to be used for the linear model
#' @param boot_log_formula The formula to be used for the logistic model
#' 
#' @return A list with one `boot_rep()` result per bootstrap sample
#' @noRd
#' 
boot_rep_par <- function(x,
                         boot_lst,
                         domain_level,
                         boot_lin_formula,
                         boot_log_formula) {
  
  p <- progressor(steps = length(x))
  
  furrr::future_map(.x = boot_lst,
                    .f = \(.x) {
                      p()
                      boot_rep(boot_samp = .x,
                               domain_level,
                               boot_lin_formula,
                               boot_log_formula)
                    },
                    .options = furrr_options(seed = TRUE))
  
}

#' Refit both models to each bootstrap sample
#' 
#' @inheritParams boot_rep_par
#' @param parallel Logical. Whether to run the refits in parallel
#' 
#' @return A list with one `boot_rep()` result per bootstrap sample
#' @noRd
fit_boot_reps <- function(boot_lst,
                          parallel,
                          domain_level,
                          boot_lin_formula,
                          boot_log_formula) {
  
  if (parallel) {
    with_progress({
      res <- boot_rep_par(x = seq_along(boot_lst),
                          boot_lst = boot_lst,
                          domain_level,
                          boot_lin_formula,
                          boot_log_formula)
    })
  } else {
    res <-
      map(.x = boot_lst,
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
  }
  
  res
  
}

#' Combine the bootstrap refits' parameters into matrices
#' 
#' @param res A list of `boot_rep()` results
#' 
#' @return A list of four matrices with one row per replicate: fixed effects
#' (`beta_lm_mat`, `beta_glm_mat`) and domain random effects (`u_lm`, `u_glm`)
#' @noRd
collect_boot_params <- function(res) {
  
  to_mat <- \(nm) res |>
    map(.f = ~ as.data.frame(t(.x[[nm]]))) |>
    list_rbind() |>
    as.matrix()
  
  beta_lm_mat <- to_mat("beta_lm")
  beta_glm_mat <- to_mat("beta_glm")
  u_lm <- to_mat("u_lm")
  u_glm <- to_mat("u_glm")
  
  # sometimes u_lm will have fewer domains once it is filtered
  # down to positive response values
  u_lm[is.na(u_lm)] <- 0

  # domains that never appear in a given bootstrap replicate (this is
  # always true of out-of-sample domains, which contribute 0 rows to every
  # replicate) get no random-effect estimate from that replicate's glmer
  # fit; treating that as a 0 random effect gives a purely synthetic
  # (population-average) prediction for them, consistent with u_lm above
  u_glm[is.na(u_glm)] <- 0
  
  list(beta_lm_mat = beta_lm_mat,
       beta_glm_mat = beta_glm_mat,
       u_lm = u_lm,
       u_glm = u_glm)
  
}




#' Sample n rows from a data.frame
#' 
#' @param .data The data.frame to sample from
#' @param n The number of rows to sample
#' @param replace Logical. Should selected rows be added back to the data.frame for the next row selection?
#' 
#' @return A data.frame with n randomly chosen rows
#' @noRd
#' 
slice_samp <- function(.data, n, replace = TRUE) {
  .data[sample(nrow(.data), n, replace = replace),]
}


#' Extract all matches from a string
#' 
#' @param string String to extract matches from
#' @param pattern Regex pattern to use for extractions
#' 
#' @return A vector of the matches
#' @noRd
str_extract_all_base <- function(string, pattern) {
  regmatches(string, gregexpr(pattern, string))
}


#' Format model parameters for mse estimation
#' 
#' @param .fit A custom list containing model fits
#' @param ref  A list containing names to be used to fill in when a model fit failed
#' 
#' @return A list containing all of the properly formated parameters
#' @noRd
#' 
mod_param_fmt <- function(.fit, ref = NULL) {
  
  if (!is.null(.fit)) {
    .lmer <- .fit$lmer
    .glmer <- .fit$glmer
    
    beta_lm <- lme4::fixef(.lmer)
    beta_glm <- lme4::fixef(.glmer)
    
    ref_lm <- lme4::ranef(.lmer)[[1]]
    ref_glm <- lme4::ranef(.glmer)[[1]]
    
    u_lm <- setNames(
      ref_lm[ ,1],
      rownames(ref_lm)
    )
    u_glm <- setNames(
      ref_glm[ ,1],
      rownames(ref_glm)
    ) 
  } else {
    lm_terms <- ref$.lm[!grepl("\\|", ref$.lm)]
    glm_terms <- ref$.glm[!grepl("\\|", ref$.glm)]
    beta_lm <- setNames(
      rep(NA, times = length(lm_terms) + 1),
      c('(Intercept)', lm_terms)
    )
    beta_glm <- setNames(
      rep(NA, times = length(glm_terms) + 1),
      c('(Intercept)', glm_terms)
    )
    u_lm <- setNames(
      rep(NA, times = length(ref$d)),
      ref$d
    )
    u_glm <- u_lm
  }
  
  list(beta_lm = beta_lm,
       beta_glm = beta_glm,
       u_lm = u_lm,
       u_glm = u_glm)
  
}


#' Perform a single bootstrap repetition
#' 
#' Warnings from the model fits (for example, convergence warnings) are
#' recorded and muffled rather than passed on, so that `saeczi()` can count
#' them and report them once.
#' 
#' @param boot_samp data.frame, An individual bootstrap sample.
#' @param domain_level Character. Domain identifier name
#' @param boot_lin_formula The formula to be used for the linear model
#' @param boot_log_formula The formula to be used for the logistic model
#' 
#' @return A list containing the properly formated model parameters from
#' fitting the two models to the sample data, plus `status` ("ok", "warning",
#' or "error"), `message` (the warning or error messages, "" if none),
#' `singular` (whether either fit is singular; NA after an error), and `time`
#' (elapsed seconds).
#' @noRd
#' 
boot_rep <- function(boot_samp,
                     domain_level,
                     boot_lin_formula,
                     boot_log_formula) {
  
  start_time <- proc.time()[["elapsed"]]
  warn_msgs <- character(0)
  error_msg <- NULL
  
  boot_samp_fit  <- tryCatch(
    withCallingHandlers(
      {
        out <- fit_zi(boot_samp,
                      boot_lin_formula,
                      boot_log_formula,
                      domain_level)
        
        ps <- mod_param_fmt(out)
        ps$singular <- lme4::isSingular(out$lmer) || lme4::isSingular(out$glmer)
        ps
      },
      warning = function(w) {
        warn_msgs <<- c(warn_msgs, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(cond) {
      
      error_msg <<- conditionMessage(cond)
      
      doms <- unique(boot_samp[[domain_level]])
      lm_cfs <- labels(terms(boot_lin_formula))
      glm_cfs <- labels(terms(boot_log_formula))
      
      ps <- mod_param_fmt(.fit = NULL,
                          ref = list(d = doms,
                                     .lm = lm_cfs,
                                     .glm = glm_cfs))
      ps$singular <- NA
      ps
      
    }
  )
  
  if (!is.null(error_msg)) {
    boot_samp_fit$status <- "error"
    boot_samp_fit$message <- error_msg
  } else if (length(warn_msgs) > 0) {
    boot_samp_fit$status <- "warning"
    boot_samp_fit$message <- paste(unique(warn_msgs), collapse = "; ")
  } else {
    boot_samp_fit$status <- "ok"
    boot_samp_fit$message <- ""
  }
  
  boot_samp_fit$time <- proc.time()[["elapsed"]] - start_time
  
  return(boot_samp_fit)
  
}


#' Format model coefficients
#' 
#' @param lmer_model lme4::lmer model object
#' @param glmer_model lme4::glmer model object
#' 
#' @return A list with all of the necessary model parameters
#' @noRd
#' 
mse_coefs <- function(lmer_model, glmer_model) {

  # from lmer model
  beta_hat <- lmer_model@beta # linear model coefficients
  model_summary_df <- data.frame(summary(lmer_model)$varcor)

  sig2_mu_hat <- model_summary_df[1, ]$vcov
  sig2_eps_hat <- subset(model_summary_df, grp == "Residual")$vcov

  # from glmer model
  alpha_1 <- glmer_model@beta

  # random-intercept variance of the logistic model, analogous to
  # sig2_mu_hat above. glmer (binomial) fits have no residual variance
  # component, so the domain random intercept is the first (only) row.
  glm_summary_df <- data.frame(summary(glmer_model)$varcor)
  sig2_b_hat <- glm_summary_df[1, ]$vcov

  b_i <- lme4::ranef(glmer_model)[[1]][,1]
  b_domain_levels <- rownames(lme4::ranef(glmer_model)[[1]])

  return(list(
    beta_hat = beta_hat, sig2_mu_hat = sig2_mu_hat,
    sig2_eps_hat = sig2_eps_hat, alpha_1 = alpha_1, sig2_b_hat = sig2_b_hat,
    b_i = b_i, domain_levels = b_domain_levels))

}

#' A function factory that causes a function to capture messages, warnings, or errors 
#' 
#' @param .f The function 
#' 
#' @return A function that works just like .f but also returns any messages, warnings, or errors that result from using the function
#' @noRd
capture_all <- function(.f){
  
  .f <- purrr::as_mapper(.f)
  
  function(...){
    
    try_out <- suppressMessages(suppressWarnings(
      try(.f(...), silent = TRUE)
    ))
    
    res <- rlang::try_fetch(
      .f(...),
      error = function(err) rlang::abort("Failed.", parent = err),
      warning = function(warn) warn,
      message = function(message) message,
    )
    
    out <- list(
      result = NULL,
      log = NULL
    )
    
    if("error" %in% class(res)) {
      stop(res$message)
    } else if (!any(c("warning", "message") %in% class(res))){
      out$result <- try_out
      out$log <- NA
    } else {
      out$result <- try_out
      out$log <- res$message
    }
    
    return(out)
    
  }
  
}

#' Checking if a param inherits a class
#' 
#' @param what What class to check if the parameter input inherits
#' @param ... The parameter input(s) to check
#' 
#' @return Nothing if the check is passed, but an error if the check fails
#' @noRd
check_inherits <- function(what, ...) {
  opts <- list(...)
  for (i in seq_along(opts)) {
    if (!is.null(opts[[i]])) {
      if (!inherits(opts[[i]], what)) {
        stop(paste0(opts[[i]], " needs to be of class ", what))
      }
    } else {
      stop("unable to check NULL objects")
    }
  } # i
  invisible(opts)
}

#' Checking if parallel functionality is properly set up
#' 
#' @param x The parameter input to check
#' @param call The caller environment to check in
#' 
#' @return Nothing if the check is passed, but an error if the check fails
#' @noRd
check_parallel <- function(x, call = rlang::caller_env()) {
  
  if (x) {
    if (eval(inherits(future::plan(), "sequential"), envir = call)) {
      message("In order for the internal processes to be run in parallel a `future::plan()` must be specified by the user")
      message("See <https://future.futureverse.org/reference/plan.html> for reference on how to use `future::plan()`")
    }
  }
  
  if (x && future::nbrOfWorkers() == 1) {
    warning("Argument `parallel` is set to true, but only one core is being used")
  }
  
  invisible(x)
}

#' Checking random effect column
#' 
#' @param pop_dat The population dataset to check
#' @param domain_level Character. The domain level identifier.
#' 
#' @return Nothing if the check is passed, but and error if it fails
#' @noRd
check_re <- function(pop_dat, samp_dat, domain_level) {
  if (!(domain_level %in% names(pop_dat))) {
    stop(paste0("Column ", domain_level, " does not exist in pop_dat"))
  }
  if (!(domain_level %in% names(samp_dat))) {
    stop(paste0("Column ", domain_level, " does not exist in samp_dat"))
  }
  if (!inherits(pop_dat[[domain_level]], "character") || !inherits(samp_dat[[domain_level]], "character")) {
    stop(paste0("Column ", domain_level, " must be of type `character` in both pop_dat and samp_dat"))
  }
}

#' Fast aggregation
#' 
#' @noRd
agg_stat <- function(vals, nms, .f) {
  agg <- tapply(vals, nms, .f)
  out <- data.frame(nms = names(agg), vals = agg)
  out
}

#' Predict with both models and return result
#' 
#' @param mod1 Linear model object
#' @param mod2 Logistic model object
#' @param estimand Character, either "means" or "totals"
#' @param .data The data to predict on
#' @param domain_level Character name of domain variable in .data
#' 
#' @returns A data frame of results
#' @noRd

collect_preds <- function(mod1, mod2, estimand, .data, domain_level, inv_transform_fun) {
  
  lin_pred <- predict(mod1, newdata = .data, allow.new.levels = TRUE)
  log_pred <- predict(mod2, newdata = .data, type = "response", allow.new.levels = TRUE)
  
  if (!is.null(inv_transform_fun)) {
    lin_pred <- inv_transform_fun(lin_pred)
  }
  
  unit_preds <- lin_pred * log_pred

  
  out <- switch(estimand,
                "means" = agg_stat(unit_preds, .data[[domain_level]], mean),
                "totals" = agg_stat(unit_preds, .data[[domain_level]], sum),
                stop())
  
  
  names(out) <- c(domain_level, "est")
  out
  
}
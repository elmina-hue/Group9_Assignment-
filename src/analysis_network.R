# 1. Group-level multilevel VAR -------------------------------------------
# VAR: each symptom is predicted from earlier values of all model symptoms.
# This model includes participant-specific intercepts and temporal effects.
run_group_network <- function(ema, ys, pr,
                              iter = 4000, warmup = 1000,
                              chains = 4, cores = 4, seed = 2026) {
  required <- unique(c("id", "day", ys))
  if (!all(required %in% names(ema))) {
    stop("Missing data columns: ",
         paste(setdiff(required, names(ema)), collapse = ", "))
  }
  
  fit_group <- bvar(
    id_col = "id", time_col = "day", y_cols = ys,
    # id_col identifies participants; time_col specifies time.
    # y_cols identifies the outcome columns.
    x_cols = NULL, re_cols = c("Intercept"), re_temporal = TRUE,
    # x_cols refers to additional covariates,so no additional covariates are added; 
    # re_cols specifies which intercept or additional predictor effects are allowed to vary between participants,
    # so random intercepts are enabled;
    # re_temporal controls whether lagged symptom relationships can differ between participants,
    # so random temporal coefficients are enabled across participants.
    K = 1, data = ema, family = "gaussian", priors = pr,
    # K = 1 requests one lag; ema supplies all participants' data.
    # Gaussian observation model; use the priors stored in pr.
    # Verify whether lag 1 means the previous observed record or
    # a fixed time interval when measurements have gaps.
    iter = iter, warmup = warmup, chains = chains,
    cores = cores, seed = seed
    # Sampling settings: 4000 iterations and 1000 warmup for each chain
    # four chains, up to four cores
    # If iter includes warmup (a common convention), this yields
    # 3000 retained draws per chain, or 12000 in total.
  )
  
  model_summary <- summary(fit_group)
  print(model_summary)
  # Print the fitted-model summary; inspect estimates and sampling diagnostics.
  # Rhat evaluates chain agreement; ESS is effective sample size.
  # Confirm exactly which diagnostics summary() reports for this package.
  
  random_effects <- extract_random_effects(fit_group)
  print(random_effects)
  # Extract random-effect information. Verify whether this function returns
  # participant deviations, their SDs, or both before interpreting the output.
  
  invisible(list(
    fit_group = fit_group,
    summary = model_summary,
    random_effects = random_effects,
    # Extract random-effect information. Verify whether this function returns
    # participant deviations, their SDs, or both before interpreting the output
  ))
}


# 2. Individual-level single-subject VAR ----------------------------------
# Fit the VAR to Participant 1's data only.
run_individual_network <- function(ema, ys, pr, participant_id = 1,
                                   iter = 4000, warmup = 1000,
                                   chains = 4, cores = 4, seed = 2026) {
  required <- unique(c("id", "day", ys))
  if (!all(required %in% names(ema))) {
    stop("Missing data columns: ",
         paste(setdiff(required, names(ema)), collapse = ", "))
  }
  if (length(participant_id) != 1L || is.na(participant_id)) {
    stop("participant_id must be one non-missing ID.")
  }
  
  individual_data <- ema[
    !is.na(ema$id) & ema$id == participant_id,
    , drop = FALSE
  ]
  if (nrow(individual_data) < 2L) {
    stop("The selected participant has fewer than two observations.")
  }
  
  fit_individual <- bvar(
    id_col = "id", time_col = "day", y_cols = ys,
    # Same identifying, time, and outcome columns as above.
    x_cols = NULL, re_cols = NULL, re_temporal = FALSE,
    # No additional covariates or participant random effects.
    # re_cols = NULL does not necessarily remove fixed intercepts.
    K = 1, data = individual_data, family = "gaussian", priors = pr,
    # One lag, Participant 1's data, Gaussian observations,
    # and the original prior specification.
    iter = iter, warmup = warmup, chains = chains,
    cores = cores, seed = seed
  )
  
  model_summary <- summary(fit_individual)
  print(model_summary)
  
  temporal <- extract_temporal(fit_individual)
  print(temporal)
  
  bayes_factors <- bf_table(fit_individual, type = "cl")
  print(bayes_factors)
  
  invisible(list(
    fit_individual = fit_individual,
    participant_id = participant_id,
    summary = model_summary,
    # Display individual-model estimates and diagnostics.
    temporal = temporal,
    # Extract a table of temporal-edge posterior summaries, not raw draws.
    bayes_factors = bayes_factors
    # Request a table with all edges.
  ))
}

###### SENSITIVITY ANALYSIS #############

# 1. Original prior: Cauchy scale = 0.50 -----------------------------------
# Re-create the original prior; this overwrites pr with the same specification.
fit_original_model <- function(ema_p1, ys,
                               iter = 4000, warmup = 1000,
                               chains = 4, cores = 4, seed = 2026) {
  pr <- set_priors(
    phi = prior(family = "cauchy", loc = 0, scale = 0.50),
    # Original temporal-coefficient prior, centered at zero.
    sd_u = prior(family = "normal", loc = 0, scale = 0.35)
    # Hold the random-effect SD prior constant throughout this analysis.
  )
  # Refit Participant 1 with the original prior; other settings remain unchanged
  fit <- bvar(
    id_col = "id", time_col = "day", y_cols = ys,
    x_cols = NULL, re_cols = NULL, re_temporal = FALSE,
    K = 1, data = ema_p1, family = "gaussian", priors = pr,
    iter = iter, warmup = warmup, chains = chains,
    cores = cores, seed = seed
  )
  
  return(fit)
}


# 2. Narrow prior: Cauchy scale = 0.35 -------------------------------------
# Narrower temporal prior: more concentrated around zero than the original.
fit_narrow_model <- function(ema_p1, ys,
                             iter = 4000, warmup = 1000,
                             chains = 4, cores = 4, seed = 2026) {
  pr_narrow <- set_priors(
    phi = prior(family = "cauchy", loc = 0, scale = 0.35),
    sd_u = prior(family = "normal", loc = 0, scale = 0.35)
  )
  # Refit Participant 1 with a narrow temporal prior; other settings remain unchanged
  fit <- bvar(
    id_col = "id", time_col = "day", y_cols = ys,
    x_cols = NULL, re_cols = NULL, re_temporal = FALSE,
    K = 1, data = ema_p1, family = "gaussian", priors = pr_narrow,
    iter = iter, warmup = warmup, chains = chains,
    cores = cores, seed = seed
  )
  
  return(fit)
}


# 3. Wide prior: Cauchy scale = 0.65 ---------------------------------------
# Wider temporal prior: allows more prior mass at larger absolute coefficients
fit_wide_model <- function(ema_p1, ys,
                           iter = 4000, warmup = 1000,
                           chains = 4, cores = 4, seed = 2026) {
  pr_wide <- set_priors(
    phi = prior(family = "cauchy", loc = 0, scale = 0.65),
    sd_u = prior(family = "normal", loc = 0, scale = 0.35)
  )
  # Refit Participant 1 with the wider temporal prior; other settings remain unchanged.
  fit <- bvar(
    id_col = "id", time_col = "day", y_cols = ys,
    x_cols = NULL, re_cols = NULL, re_temporal = FALSE,
    K = 1, data = ema_p1, family = "gaussian", priors = pr_wide,
    iter = iter, warmup = warmup, chains = chains,
    cores = cores, seed = seed
  )
  
  return(fit)
}

### Assume we are interested in the temporal coefficient between sleep problems and feeling tired, since it might consitute a potential treatment target ###
### To assess if this is a treatment target, we must extract the temporal coefficient between these two nodes ###
# Extract sleep -> tired edge across all three models

#extract lag1sleep -> tired under the original model with scale = 0.5
extract_original_temporal <- function(fit_p1_original,
                                      predictor = NULL, outcome = NULL) {
  coefficients <- extract_temporal(fit_p1_original)
  
  if (is.null(predictor) != is.null(outcome)) {
    stop("Supply both predictor and outcome, or leave both NULL.")
  }
  if (!is.null(predictor)) {
    if (length(predictor) != 1L || length(outcome) != 1L ||
        is.na(predictor) || is.na(outcome)) {
      stop("predictor and outcome must each be one non-missing name.")
    }
    if (!all(c("predictor", "outcome") %in% names(coefficients))) {
      stop("extract_temporal() must return predictor and outcome columns.")
    }
    keep <- !is.na(coefficients$predictor) & !is.na(coefficients$outcome) &
      coefficients$predictor == paste0("lag1_", predictor) &
      coefficients$outcome == outcome
    coefficients <- coefficients[keep, , drop = FALSE]
    if (nrow(coefficients) == 0L) stop("No matching lag-1 temporal edge found.")
  }
  
  return(coefficients)
}

#extract lag1sleep -> tired under the narrow model with scale = 0.35
extract_narrow_temporal <- function(fit_p1_narrow,
                                    predictor = NULL, outcome = NULL) {
  coefficients <- extract_temporal(fit_p1_narrow)
  
  if (is.null(predictor) != is.null(outcome)) {
    stop("Supply both predictor and outcome, or leave both NULL.")
  }
  if (!is.null(predictor)) {
    if (length(predictor) != 1L || length(outcome) != 1L ||
        is.na(predictor) || is.na(outcome)) {
      stop("predictor and outcome must each be one non-missing name.")
    }
    if (!all(c("predictor", "outcome") %in% names(coefficients))) {
      stop("extract_temporal() must return predictor and outcome columns.")
    }
    keep <- !is.na(coefficients$predictor) & !is.na(coefficients$outcome) &
      coefficients$predictor == paste0("lag1_", predictor) &
      coefficients$outcome == outcome
    coefficients <- coefficients[keep, , drop = FALSE]
    if (nrow(coefficients) == 0L) stop("No matching lag-1 temporal edge found.")
  }
  
  return(coefficients)
}

#extract lag1sleep -> tired under the wide model with scale = 0.65
extract_wide_temporal <- function(fit_p1_wide,
                                  predictor = NULL, outcome = NULL) {
  coefficients <- extract_temporal(fit_p1_wide)
  
  if (is.null(predictor) != is.null(outcome)) {
    stop("Supply both predictor and outcome, or leave both NULL.")
  }
  if (!is.null(predictor)) {
    if (length(predictor) != 1L || length(outcome) != 1L ||
        is.na(predictor) || is.na(outcome)) {
      stop("predictor and outcome must each be one non-missing name.")
    }
    if (!all(c("predictor", "outcome") %in% names(coefficients))) {
      stop("extract_temporal() must return predictor and outcome columns.")
    }
    keep <- !is.na(coefficients$predictor) & !is.na(coefficients$outcome) &
      coefficients$predictor == paste0("lag1_", predictor) &
      coefficients$outcome == outcome
    coefficients <- coefficients[keep, , drop = FALSE]
    if (nrow(coefficients) == 0L) stop("No matching lag-1 temporal edge found.")
  }
  
  return(coefficients)
}

# Return posterior summaries for the three prior specifications.
# Combine the edge summaries into one table.
compare_posteriors <- function(edge_narrow, edge_original, edge_wide) {
  sensitivity_table <- rbind(
    data.frame(prior = "narrow (0.35)", edge_narrow),
    data.frame(prior = "original (0.50)", edge_original),
    data.frame(prior = "wide (0.65)", edge_wide)
  )
  # Keep all rows but only these five columns, in the specified order.
  sensitivity_table <- sensitivity_table[
    , c("prior", "mean", "median", "ci_lower", "ci_upper"),
    drop = FALSE
  ]
  rownames(sensitivity_table) <- NULL
  
  return(sensitivity_table)
}

# Return Bayes factors for sleep -> tired across the three models.
# Ho represents the hypothesis that the edge equals 0
# H1 represents the hypothesis that the edge differs from 0
# The package or script providing bf_table() must already be loaded.
prior_labels <- c("narrow (0.35)", "original (0.50)", "wide (0.65)")
compare_bayes_factors <- function(fit_p1_narrow, fit_p1_original, fit_p1_wide,
                                  variable = "sleep", outcome = "tired") {
  fits <- list(fit_p1_narrow, fit_p1_original, fit_p1_wide)
  prior_labels <- c("narrow (0.35)", "original (0.50)", "wide (0.65)")
  
  bf_values <- vapply(seq_along(fits), function(i) {
    bf <- bf_table(fits[[i]], type = "cl", variable = variable)
    
    if (!all(c("outcome", "BF10") %in% names(bf))) {
      stop("bf_table() must return columns named 'outcome' and 'BF10'.")
    }
    
    selected <- which(bf$outcome == outcome)
    if (length(selected) != 1L) {
      stop(sprintf(
        "Expected exactly one row for outcome '%s' under prior '%s'; found %d.",
        outcome, prior_labels[i], length(selected)
      ))
    }
    
    return(as.numeric(bf$BF10[selected]))
  }, numeric(1))
  # BF10 conventionally denotes evidence for H1 relative to H0.
  return(data.frame(prior = prior_labels, BF10 = bf_values))
}

# 1. Extract posterior draws for P1 (also usable for the group model).
get_draws_p1 <- function(fit_p1, parameter = "phi") {
  return(extract_draws(fit_p1, parameter = parameter))
}

# 2. Extract the sleep -> tired edge (also usable for group draws).
get_p1_sleep_tired <- function(draws_p1, edge = "phi[3,2]") {
  return(draws_p1[, edge])
}

# 3. Compute draw-wise differences: P1 minus group.
get_delta <- function(p1_sleep_tired, group_sleep_tired) {
  if (length(p1_sleep_tired) != length(group_sleep_tired)) {
    stop("P1 and group draws must have the same length; R recycling is not allowed.")
  }
  if (length(p1_sleep_tired) == 0L) {
    stop("Posterior draw vectors must not be empty.")
  }
  
  return(p1_sleep_tired - group_sleep_tired)
}

# 4. Summarize the posterior difference using the original 5%/95% bounds.
get_posterior_difference <- function(delta_sleep_tired) {
  if (!is.numeric(delta_sleep_tired) || length(delta_sleep_tired) == 0L ||
      any(!is.finite(delta_sleep_tired))) {
    stop("delta_sleep_tired must contain nonempty, finite numeric draws.")
  }
  
  return(c(
    median_difference = median(delta_sleep_tired),
    mean = mean(delta_sleep_tired),
    lower = unname(quantile(delta_sleep_tired, 0.05)),
    upper = unname(quantile(delta_sleep_tired, 0.95)),
    probability_p1_smaller = mean(delta_sleep_tired < 0),
    probability_p1_larger = mean(delta_sleep_tired > 0)
  ))
}


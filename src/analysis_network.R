# 1. Group-level multilevel VAR -------------------------------------------
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
    x_cols = NULL, re_cols = c("Intercept"), re_temporal = TRUE,
    K = 1, data = ema, family = "gaussian", priors = pr,
    iter = iter, warmup = warmup, chains = chains,
    cores = cores, seed = seed
  )
  
  model_summary <- summary(fit_group)
  print(model_summary)
  
  random_effects <- extract_random_effects(fit_group)
  print(random_effects)
  
  bayes_factors <- bf_table(fit_group, type = "temporal")
  print(bayes_factors)
  
  invisible(list(
    fit_group = fit_group,
    summary = model_summary,
    random_effects = random_effects,
    bayes_factors = bayes_factors
  ))
}


# 2. Individual-level single-subject VAR ----------------------------------
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
    x_cols = NULL, re_cols = NULL, re_temporal = FALSE,
    K = 1, data = individual_data, family = "gaussian", priors = pr,
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
    temporal = temporal,
    bayes_factors = bayes_factors
  ))
}

###### SENSITIVITY ANALYSIS #############

# 1. Original prior: Cauchy scale = 0.50 -----------------------------------
fit_original_model <- function(ema_p1, ys,
                               iter = 4000, warmup = 1000,
                               chains = 4, cores = 4, seed = 2026) {
  pr <- set_priors(
    phi = prior(family = "cauchy", loc = 0, scale = 0.50),
    sd_u = prior(family = "normal", loc = 0, scale = 0.35)
  )
  
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
fit_narrow_model <- function(ema_p1, ys,
                             iter = 4000, warmup = 1000,
                             chains = 4, cores = 4, seed = 2026) {
  pr_narrow <- set_priors(
    phi = prior(family = "cauchy", loc = 0, scale = 0.35),
    sd_u = prior(family = "normal", loc = 0, scale = 0.35)
  )
  
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
fit_wide_model <- function(ema_p1, ys,
                           iter = 4000, warmup = 1000,
                           chains = 4, cores = 4, seed = 2026) {
  pr_wide <- set_priors(
    phi = prior(family = "cauchy", loc = 0, scale = 0.65),
    sd_u = prior(family = "normal", loc = 0, scale = 0.35)
  )
  
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
compare_posteriors <- function(edge_narrow, edge_original, edge_wide) {
  sensitivity_table <- rbind(
    data.frame(prior = "narrow (0.35)", edge_narrow),
    data.frame(prior = "original (0.50)", edge_original),
    data.frame(prior = "wide (0.65)", edge_wide)
  )
  
  sensitivity_table <- sensitivity_table[
    , c("prior", "mean", "median", "ci_lower", "ci_upper"),
    drop = FALSE
  ]
  rownames(sensitivity_table) <- NULL
  
  return(sensitivity_table)
}

# Return Bayes factors for sleep -> tired across the three models.
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


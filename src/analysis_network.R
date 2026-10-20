ys <- c("down", "tired", "sleep", "nervous", "worry", "irritable")
pr <- set_priors(phi = prior(family = "cauchy", loc = 0, scale = 0.5), sd_u = prior(family = "normal", loc = 0, scale = 0.35))

## (1) GROUP-LEVEL multilevel VAR: random intercepts + random temporal effects
fit_group <- bvar(id_col = "id", time_col = "day", y_cols = ys,
                  x_cols = NULL, re_cols = c("Intercept"), re_temporal = TRUE,
                  K = 1, data = ema, family = "gaussian", priors = pr,
                  iter = 4000, warmup = 1000, chains = 4, cores = 4, seed = 2026)
summary(fit_group) # fixed effects, random-effect SDs, Rhat, ESS
qgraph(extract_network_matrix(fit_group)) # group-level temporal network
extract_random_effects(fit_group) # between-person SD of each edge
bf_table(fit_group, type = "temporal") # Bayes factors for the edges

## (2) IDIOGRAPHIC single-subject VAR for Participant 1
ema_p1 <- subset(ema, id == 1)
fit_p1 <- bvar(id_col = "id", time_col = "day", y_cols = ys,
               x_cols = NULL, re_cols = NULL, re_temporal = FALSE,
               K = 1, data = ema_p1, family = "gaussian", priors = pr,
               iter = 4000, warmup = 1000, chains = 4, cores = 4, seed = 2026)
summary(fit_p1)
qgraph(extract_network_matrix(fit_p1))
extract_temporal(fit_p1) # posterior means + intervals per edge
bf_table(fit_p1, type = "cl", variable = "sleep")
bf_table(fit_p1, type = "cl")

###### SENSITIVITY ANALYSIS #############

# Priors
pr <- set_priors(
  phi = prior(family = "cauchy", loc = 0, scale = 0.5),
  sd_u = prior(family = "normal", loc = 0, scale = 0.35)
)

pr_narrow <- set_priors(
  phi = prior(family = "cauchy", loc = 0, scale = 0.35),
  sd_u = prior(family = "normal", loc = 0, scale = 0.35)
)

pr_wide <- set_priors(
  phi = prior(family = "cauchy", loc = 0, scale = 0.65),
  sd_u = prior(family = "normal", loc = 0, scale = 0.35)
)

# Refit P1 under all three priors
fit_p1_original <- bvar(
  id_col = "id", time_col = "day", y_cols = ys,
  x_cols = NULL, re_cols = NULL, re_temporal = FALSE,
  K = 1, data = ema_p1, family = "gaussian", priors = pr,
  iter = 4000, warmup = 1000, chains = 4, cores = 4, seed = 2026
)

fit_p1_narrow <- bvar(
  id_col = "id", time_col = "day", y_cols = ys,
  x_cols = NULL, re_cols = NULL, re_temporal = FALSE,
  K = 1, data = ema_p1, family = "gaussian", priors = pr_narrow,
  iter = 4000, warmup = 1000, chains = 4, cores = 4, seed = 2026
)

fit_p1_wide <- bvar(
  id_col = "id", time_col = "day", y_cols = ys,
  x_cols = NULL, re_cols = NULL, re_temporal = FALSE,
  K = 1, data = ema_p1, family = "gaussian", priors = pr_wide,
  iter = 4000, warmup = 1000, chains = 4, cores = 4, seed = 2026
)

### Assume we are interested in the temporal coefficient between sleep problems and feeling tired, since it might consitute a potential treatment target ###
### To assess if this is a treatment target, we must extract the temporal coefficient between these two nodes ###
# Extract sleep -> tired edge across all three models
edge_original <- extract_temporal(fit_p1_original)[
  extract_temporal(fit_p1_original)$predictor == "lag1_sleep" &
    extract_temporal(fit_p1_original)$outcome == "tired", ]

edge_narrow <- extract_temporal(fit_p1_narrow)[
  extract_temporal(fit_p1_narrow)$predictor == "lag1_sleep" &
    extract_temporal(fit_p1_narrow)$outcome == "tired", ]

# Compare posteriors in a table
sensitivity_table <- rbind(
  data.frame(prior = "narrow (0.35)", edge_narrow),
  data.frame(prior = "original (0.50)", edge_original),
  data.frame(prior = "wide (0.65)", edge_wide)
)
sensitivity_table
print(sensitivity_table[, c("prior", "mean", "median", "q5", "q95")])

# BFs for sleep -> tired across all three models
bf_original <- bf_table(fit_p1_original, type = "cl", variable = "sleep")
bf_narrow   <- bf_table(fit_p1_narrow,   type = "cl", variable = "sleep")
bf_wide     <- bf_table(fit_p1_wide,     type = "cl", variable = "sleep")

bf_comparison <- data.frame(
  prior = c("narrow (0.35)", "original (0.50)", "wide (0.65)"),
  BF10  = c(
    bf_narrow$BF10[bf_narrow$outcome == "tired"],
    bf_original$BF10[bf_original$outcome == "tired"],
    bf_wide$BF10[bf_wide$outcome == "tired"]
  )
)
print(bf_comparison)
edge_wide <- extract_temporal(fit_p1_wide)[
  extract_temporal(fit_p1_wide)$predictor == "lag1_sleep" &
    extract_temporal(fit_p1_wide)$outcome == "tired", ]
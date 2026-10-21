ys <- c("down", "tired", "sleep", "anhedonia", "worthlessness","concentration", "appetite", "psychomotor", "suicidal")
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

edge_wide <- extract_temporal(fit_p1_wide)[
  extract_temporal(fit_p1_wide)$predictor == "lag1_sleep" &
    extract_temporal(fit_p1_wide)$outcome == "tired", ]

# Compare posteriors in a table
sensitivity_table <- rbind(
  data.frame(prior = "narrow (0.35)", edge_narrow),
  data.frame(prior = "original (0.50)", edge_original),
  data.frame(prior = "wide (0.65)", edge_wide)
)
sensitivity_table = sensitivity_table[, c("prior", "mean", "median", "ci_lower", "ci_upper")]
print(sensitivity_table)

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

# Overlay posterior plot using base R
means  <- c(edge_narrow$mean, edge_original$mean, edge_wide$mean)
q5s    <- c(edge_narrow$ci_lower,   edge_original$ci_lower,   edge_wide$ci_lower)
q95s   <- c(edge_narrow$ci_upper,  edge_original$ci_upper,  edge_wide$ci_upper)
labels <- c("Narrow (0.35)", "Original (0.50)", "Wide (0.65)")
colors <- c("steelblue", "black", "firebrick")

plot(1:3, means, ylim = c(min(q5s) - 0.05, max(q95s) + 0.05),
     xlim = c(0.5, 3.5), xaxt = "n", pch = 19, col = colors, cex = 1.5,
     ylab = "Posterior estimate", xlab = "Prior scale",
     main = "Sensitivity analysis: sleep -> tired (P1)")
axis(1, at = 1:3, labels = labels)
abline(h = 0, lty = 2, col = "gray")
for (i in 1:3) {
  lines(c(i, i), c(q5s[i], q95s[i]), col = colors[i], lwd = 2)
}
legend("topright", legend = labels, col = colors, lwd = 2, pch = 19)



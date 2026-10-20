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


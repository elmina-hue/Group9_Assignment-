# ============================================================
# Simulate a multilevel daily symptom dataset based on PHQ-9 domains
# ============================================================
# This is a daily EMA-style adaptation of the nine PHQ-9 symptom
# domains, NOT the standard PHQ-9 administration, which asks about
# symptoms over the previous two weeks.
#
# All coefficient values and cross-lagged relationships below are
# illustrative simulation assumptions, not empirical estimates or
# established clinical relationships.
#
# Run this script from the repository root. It creates:
#   data/ema_simulated.csv
#   data/rating_frequencies.csv

set.seed(2026)  # Makes the simulation reproducible.

N <- 25    # Number of participants
Tt <- 100  # Number of days per participant
p <- 9     # Number of symptom domains

# Keep the first three variables in the same order as the colleague's
# analysis code, which refers to sleep and tired by their indices.
symptoms <- c(
  "down", "tired", "sleep",
  "anhedonia", "worthlessness", "concentration",
  "appetite", "psychomotor", "suicidal"
)

# ============================================================
# 1. Define the group-level temporal network
# ============================================================
# Phi_group is the group-level VAR(1) coefficient matrix.
#
# Rows represent the symptom being predicted on day t.
# Columns represent the symptom measured on the previous day (t - 1).
# Thus, Phi_group["down", "sleep"] is the coefficient for the
# previous day's sleep domain predicting the next day's low mood,
# conditional on the other symptoms in the model.
#
# Diagonal entries are autoregressive (carry-over) coefficients.
# Off-diagonal entries are cross-lagged relationships.
#
# The chosen relationships and coefficient values are illustrative.

Phi_group <- matrix(
  0,
  nrow = p,
  ncol = p,
  dimnames = list(symptoms, symptoms)
)

# Autoregressive carry-over coefficients, in symptom order.
diag(Phi_group) <- c(
  0.30, 0.30, 0.35,
  0.30, 0.35, 0.30,
  0.30, 0.30, 0.30
)

# Illustrative cross-lagged relationships.
Phi_group["down", "sleep"] <- 0.25
Phi_group["down", "tired"] <- 0.15
Phi_group["anhedonia", "worthlessness"] <- 0.20
Phi_group["worthlessness", "anhedonia"] <- 0.15
Phi_group["concentration", "worthlessness"] <- 0.15

# Additional illustrative links so all nine symptom domains are
# connected in the overall network; not every pair has an effect.
Phi_group["tired", "appetite"] <- 0.15
Phi_group["tired", "psychomotor"] <- 0.15
Phi_group["suicidal", "worthlessness"] <- 0.20
Phi_group["down", "suicidal"] <- 0.15

# ============================================================
# 2. Define between-person variability
# ============================================================
# Tau gives the standard deviation of participant-specific
# deviations from each group-level coefficient.
# These are illustrative simulation settings.

Tau <- matrix(
  0.08,
  nrow = p,
  ncol = p,
  dimnames = list(symptoms, symptoms)
)

# Allow greater between-person variability for selected coefficients.
Tau["down", "sleep"] <- 0.15
Tau["anhedonia", "worthlessness"] <- 0.15

# ============================================================
# 3. Check temporal stability
# ============================================================
# A VAR process is stable when the largest absolute eigenvalue of
# its coefficient matrix is below 1. We use a stricter threshold
# of 0.95 to avoid matrices close to the stability boundary.

stationary <- function(M) {
  max(Mod(eigen(M, only.values = TRUE)$values)) < 0.95
}

# ============================================================
# 4. Set baseline and innovation variability
# ============================================================
# The latent values are shifted by mu, rounded, then clipped to
# the PHQ-9 item response range of 0 to 3.
# This is a simplified ordinal-data simulation, not a threshold model.

mu <- 1.5
innov_sd <- 0.70

# Store each participant's data separately before combining them.
long <- vector("list", N)

# ============================================================
# 5. Generate each participant's network and daily observations
# ============================================================
# Each participant receives a coefficient matrix Phi_i by adding
# random deviations to Phi_group.
#
# Participant 1 is deliberately made different from the group:
#   - sleep -> low mood is set to zero;
#   - worthlessness -> anhedonia is strengthened;
#   - anhedonia -> concentration is added/strengthened;
#   - worthlessness -> suicidal thoughts is strengthened.
#
# These are designed examples of individual differences, not
# estimates from real participant data.
#
# Each matrix is redrawn until it passes the stability check.

for (i in seq_len(N)) {

  repeat {

    # Draw one deviation for every matrix element, using Tau as
    # the standard deviation for the corresponding coefficient.
    U <- matrix(
      rnorm(p * p, mean = 0, sd = as.vector(Tau)),
      nrow = p,
      ncol = p
    )

    # Add participant-specific deviations to the group coefficients.
    Phi_i <- Phi_group + U

    # Apply the planned individual differences to Participant 1.
    if (i == 1) {
      Phi_i["down", "sleep"] <- 0.00
      Phi_i["anhedonia", "worthlessness"] <- 0.45
      Phi_i["concentration", "anhedonia"] <- 0.30
      Phi_i["suicidal", "worthlessness"] <- 0.35
    }

    # Accept this matrix only if it passes the stability check.
    if (stationary(Phi_i)) break
  }

  # ----------------------------------------------------------
  # Simulate daily latent symptom values
  # ----------------------------------------------------------
  # Y contains continuous latent values before conversion to ratings.
  # The first day is initialised with random noise. Each following
  # day depends on the previous day's symptoms plus new random noise.

  Y <- matrix(0, nrow = Tt, ncol = p)
  Y[1, ] <- rnorm(p, mean = 0, sd = innov_sd)

  for (t in 2:Tt) {
    Y[t, ] <- Phi_i %*% Y[t - 1, ] +
      rnorm(p, mean = 0, sd = innov_sd)
  }

  # ----------------------------------------------------------
  # Convert latent values to PHQ-9-style item ratings (0 to 3)
  # ----------------------------------------------------------
  # 1. Add the baseline mean.
  # 2. Round to the nearest whole number.
  # 3. Clip values below 0 to 0 and above 3 to 3.
  #
  # This is a simplified procedure and does not explicitly model
  # ordinal response thresholds.

  obs <- round(mu + Y)
  obs[obs < 0] <- 0
  obs[obs > 3] <- 3

  # Save participant ID, day, and the nine symptom-domain ratings.
  long[[i]] <- data.frame(
    id = i,
    day = 1:Tt,
    obs
  )
}

# ============================================================
# 6. Combine all participants into one dataset
# ============================================================
# ema has one row per participant-day, with columns for participant
# ID, day, and the nine symptom-domain ratings.

ema <- do.call(rbind, long)
colnames(ema) <- c("id", "day", symptoms)

# Basic checks: expected dimensions, variable names, and rating range.
stopifnot(nrow(ema) == N * Tt)
stopifnot(all(c("id", "day", symptoms) %in% names(ema)))
stopifnot(all(as.matrix(ema[, symptoms]) %in% 0:3))

# ============================================================
# 7. Calculate the frequency of each rating
# ============================================================
# Count ratings 0, 1, 2, and 3 for every symptom. Including all
# levels ensures that ratings with zero occurrences still appear.
# Proportions are calculated separately for each symptom.

rating_frequencies <- do.call(
  rbind,
  lapply(symptoms, function(s) {

    counts <- table(
      factor(ema[[s]], levels = 0:3)
    )

    data.frame(
      symptom = s,
      rating = 0:3,
      count = as.integer(counts),
      proportion = as.integer(counts) / nrow(ema)
    )
  })
)

print(rating_frequencies)

# Display only the minimum (0) and maximum (3) ratings.
boundary_frequencies <- rating_frequencies[
  rating_frequencies$rating %in% c(0, 3),
]

print(boundary_frequencies)

# Check that each symptom has exactly N * Tt ratings in total.
ratings_per_symptom <- aggregate(
  count ~ symptom,
  data = rating_frequencies,
  FUN = sum
)
stopifnot(all(ratings_per_symptom$count == N * Tt))

# ============================================================
# 8. Save the generated data and frequency table
# ============================================================
# These paths assume the script is run from the repository root.
# The folder is created if it does not already exist.

dir.create("data", showWarnings = FALSE, recursive = TRUE)

write.csv(
  ema,
  "data/ema_simulated.csv",
  row.names = FALSE
)

write.csv(
  rating_frequencies,
  "data/rating_frequencies.csv",
  row.names = FALSE
)


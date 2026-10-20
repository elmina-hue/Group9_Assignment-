
# ============================================================
# Simulate a multilevel EMA dataset
# ============================================================
# Generate daily symptom ratings for 25 participants over 100 days.
# The simulated data demonstrate how temporal symptom networks can
# differ between a group and an individual. They are not real
# participant or clinical data.
#
# The coefficients are illustrative choices for this simulation.
# They should not be interpreted as established empirical effects
# unless supported by relevant research.

set.seed(2026)  # Makes the simulation reproducible.

N <- 25    # Number of participants
Tt <- 100  # Number of days per participant
p <- 6     # Number of symptoms

symptoms <- c(
  "down", "tired", "sleep",
  "nervous", "worry", "irritable"
)

# ============================================================
# 1. Define the group-level temporal network
# ============================================================
# Phi_group is the group-level VAR(1) coefficient matrix.
#
# Rows represent symptoms being predicted on day t.
# Columns represent symptoms measured on the previous day (t - 1).
#
# For example, Phi_group["down", "sleep"] = 0.25 means that
# previous-day sleep contributes positively to the prediction of
# next-day low mood, conditional on the other symptoms in the model.
#
# Diagonal entries represent autoregressive (carry-over) effects:
# how strongly each symptom predicts itself on the following day.
#
# Off-diagonal entries represent cross-lagged relationships between
# different symptoms.
#
# The selected symptom pairs and coefficient values below are
# illustrative. They are not claimed to be empirical estimates.

Phi_group <- matrix(
  0,
  nrow = p,
  ncol = p,
  dimnames = list(symptoms, symptoms)
)

# Autoregressive carry-over coefficients, in symptom order.
diag(Phi_group) <- c(
  0.30, 0.30, 0.35,
  0.30, 0.35, 0.30
)

# Illustrative cross-lagged relationships.
Phi_group["down", "sleep"] <- 0.25
Phi_group["down", "tired"] <- 0.15
Phi_group["nervous", "worry"] <- 0.20
Phi_group["worry", "nervous"] <- 0.15
Phi_group["irritable", "worry"] <- 0.15

# ============================================================
# 2. Define between-person variability
# ============================================================
# Tau specifies the standard deviation of the random deviation
# added to each group-level coefficient for each participant.
#
# Larger values allow more variation between participants for
# the corresponding temporal relationship.
# These values are illustrative simulation settings.

Tau <- matrix(
  0.08,
  nrow = p,
  ncol = p,
  dimnames = list(symptoms, symptoms)
)

# Allow greater between-person variability for selected coefficients.
Tau["down", "sleep"] <- 0.15
Tau["nervous", "worry"] <- 0.15

# ============================================================
# 3. Check temporal stability
# ============================================================
# A VAR process is stable when the largest absolute eigenvalue
# of its coefficient matrix is below 1.
#
# Here we use a stricter threshold of 0.95 to avoid matrices that
# are too close to the stability boundary.
#
# The function returns TRUE if the matrix passes the check and
# FALSE otherwise. Unstable matrices are redrawn in the loop below.

stationary <- function(M) {
  max(Mod(eigen(M, only.values = TRUE)$values)) < 0.95
}

# ============================================================
# 4. Set the baseline and random innovation variability
# ============================================================
# mu is the baseline level around which latent symptom values vary.
# innov_sd controls the standard deviation of the random innovations
# added at each time point.

mu <- 2.5
innov_sd <- 0.70

# Store each participant's data separately before combining them.
long <- vector("list", N)

# ============================================================
# 5. Generate each participant's network and daily observations
# ============================================================
# Each participant receives their own coefficient matrix, Phi_i.
# It is created by adding random deviations to the group matrix.
# This allows the participants' temporal networks to differ.
#
# Participant 1 is deliberately made different from the group:
#   - sleep -> low mood is set to zero;
#   - worry -> nervousness is strengthened;
#   - nervousness -> irritability is added.
#
# These changes create an example of an individual network that
# departs from the group network. They are design choices for the
# demonstration, not estimates from real participant data.
#
# The repeat loop keeps generating a participant-specific matrix
# until it passes the stability check.

for (i in seq_len(N)) {
  
  repeat {
    
    # Generate random coefficient deviations using Tau as the
    # standard deviation for each corresponding matrix element.
    U <- matrix(
      rnorm(p * p, mean = 0, sd = as.vector(Tau)),
      nrow = p,
      ncol = p
    )
    
    # Combine the group-level coefficients with individual deviations.
    Phi_i <- Phi_group + U
    
    # Apply the planned individual differences to Participant 1.
    if (i == 1) {
      
      Phi_i["down", "sleep"] <- 0.00
      Phi_i["nervous", "worry"] <- 0.45
      Phi_i["irritable", "nervous"] <- 0.30
    }
    
    # Accept this matrix only if it passes the stability check.
    if (stationary(Phi_i)) break
  }
  
  # ----------------------------------------------------------
  # Simulate daily latent symptom values
  # ----------------------------------------------------------
  # Y stores continuous values before conversion to 1–4 ratings.
  # The first day is initialised using random noise.
  # Each subsequent day depends on the previous day's values
  # through the participant's coefficient matrix, plus new noise.
  
  Y <- matrix(0, nrow = Tt, ncol = p)
  
  Y[1, ] <- rnorm(p, mean = 0, sd = innov_sd)
  
  for (t in 2:Tt) {
    Y[t, ] <- Phi_i %*% Y[t - 1, ] +
      rnorm(p, mean = 0, sd = innov_sd)
  }
  
  # ----------------------------------------------------------
  # Convert latent values to ratings between 1 and 4
  # ----------------------------------------------------------
  # Keep the current rating-generation method:
  #   1. Add the baseline mean.
  #   2. Round to the nearest whole number.
  #   3. Clip values below 1 to 1 and above 4 to 4.
  #
  # This is a simplified procedure and does not explicitly model
  # ordinal response thresholds.
  
  obs <- round(mu + Y)
  obs[obs < 1] <- 1
  obs[obs > 4] <- 4
  
  # Save this participant's ID, day, and symptom ratings.
  long[[i]] <- data.frame(
    id = i,
    day = 1:Tt,
    obs
  )
}

# ============================================================
# 6. Combine all participants into one dataset
# ============================================================
# ema contains one row per participant-day, with columns for
# participant ID, day, and the six symptom ratings.

ema <- do.call(rbind, long)
colnames(ema) <- c("id", "day", symptoms)

# ============================================================
# 7. Calculate the frequency of each rating
# ============================================================
# Count how often each possible rating (1, 2, 3, and 4) occurs
# for each symptom across all participant-day observations.
#
# Including all four levels ensures that a rating with zero
# occurrences is still shown in the table.
#
# Proportions are calculated separately for each symptom.

rating_frequencies <- do.call(
  rbind,
  lapply(symptoms, function(s) {
    
    counts <- table(
      factor(ema[[s]], levels = 1:4)
    )
    
    data.frame(
      symptom = s,
      rating = 1:4,
      count = as.integer(counts),
      proportion = as.integer(counts) / nrow(ema)
    )
  })
)

# Display the complete rating-frequency table in R.
print(rating_frequencies)

# Display only ratings at the lower and upper boundaries (1 and 4).
boundary_frequencies <- rating_frequencies[
  rating_frequencies$rating %in% c(1, 4),
]

print(boundary_frequencies)

# ============================================================
# 8. Save the generated data and frequency table
# ============================================================
# The repository already contains a data/ folder.
# Run this script from the repository root so these relative paths
# point to that folder.
#
# write.csv() creates these CSV files if they do not exist yet
# and overwrites them if they already exist.

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



## ---- Simulate the multilevel EMA dataset (run as-is) -------------------
set.seed(2026)
N <- 25 # participants
Tt <- 100 # days per participant
p <- 6 # symptoms
symptoms <- c("down", "tired", "sleep", "nervous", "worry", "irritable")
# Group-level temporal (VAR-1) coefficient matrix.
# Phi[j, k] = effect of symptom k on day t on symptom j on day t+1.
Phi_group <- matrix(0, p, p, dimnames = list(symptoms, symptoms))
diag(Phi_group) <- c(0.30, 0.30, 0.35, 0.30, 0.35, 0.30) # carry-over
Phi_group["down", "sleep"] <- 0.25
Phi_group["down", "tired"] <- 0.15
Phi_group["nervous", "worry"] <- 0.20
Phi_group["worry", "nervous"] <- 0.15
Phi_group["irritable", "worry"] <- 0.15
# Between-person SD of each coefficient (how much people differ).
Tau <- matrix(0.08, p, p, dimnames = list(symptoms, symptoms))
Tau["down", "sleep"] <- 0.15
Tau["nervous", "worry"] <- 0.15
stationary <- function(M) max(Mod(eigen(M, only.values = TRUE)$values)) < 0.95
mu <- 2.5; innov_sd <- 0.70
long <- vector("list", N)
for (i in seq_len(N)) {
  repeat {
    U <- matrix(rnorm(p * p, 0, as.vector(Tau)), p, p)
    Phi_i <- Phi_group + U
    if (i == 1) { # Participant 1 is idiographic
      Phi_i["down", "sleep"] <- 0.00
      Phi_i["nervous", "worry"] <- 0.45
      Phi_i["irritable", "nervous"] <- 0.30
    }
    if (stationary(Phi_i)) break
  }
  Y <- matrix(0, Tt, p); Y[1, ] <- rnorm(p, 0, innov_sd)
  for (t in 2:Tt) Y[t, ] <- Phi_i %*% Y[t - 1, ] + rnorm(p, 0, innov_sd)
  obs <- round(mu + Y); obs[obs < 1] <- 1; obs[obs > 4] <- 4
  long[[i]] <- data.frame(id = i, day = 1:Tt, obs)
}
ema <- do.call(rbind, long)
colnames(ema) <- c("id", "day", symptoms)

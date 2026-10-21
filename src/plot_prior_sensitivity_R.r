# plot_prior_sensitivity.R
# Source this file to define plot_prior_sensitivity().
# Then use plot_prior_sensitivity(edge_narrow, edge_original, edge_wide)
# Sourcing does not draw a plot or require edge summary objects to exist.
# Uses only packages distributed with R; no additional installation is needed.
# Each input must contain numeric scalar fields:
# mean, ci_lower, ci_upper.
# The interval level is determined by the supplied bounds, not this function.
# For example, supply interval_label = "90% credible interval" if appropriate.

plot_prior_sensitivity <- function(
    edge_narrow,
    edge_original,
    edge_wide,
    title = "How sensitive is the estimate to the prior?",
    subtitle = "Participant 1 | Sleep to tired | Temporal coefficient",
    interval_label = "Credible interval",
    digits = 3
) {
  edges <- list(edge_narrow, edge_original, edge_wide)
  required <- c("mean", "ci_lower", "ci_upper")
  
  get_value <- function(x, field) {
    value <- x[[field]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
      stop("Each summary must contain one finite numeric value for: ", field)
    }
    as.numeric(value)
  }
  for (x in edges) {
    if (!all(required %in% names(x))) {
      stop("Each edge summary must contain mean, ci_lower, and ci_upper.")
    }
  }
  
  means <- vapply(edges, get_value, numeric(1), field = "mean")
  lower <- vapply(edges, get_value, numeric(1), field = "ci_lower")
  upper <- vapply(edges, get_value, numeric(1), field = "ci_upper")
  if (any(lower > upper)) stop("ci_lower must not exceed ci_upper.")
  
  labels <- c("Narrow", "Original", "Wide")
  scales <- c("Prior scale: 0.35", "Prior scale: 0.50", "Prior scale: 0.65")
  colors <- c("#148A87", "#526BB3", "#C25773")
  y <- c(3, 2, 1)
  ink <- "#20334A"
  muted <- "#65768A"
  
  # Include zero so all three estimates share a visible reference point.
  bounds <- range(c(lower, upper, means, 0))
  span <- diff(bounds)
  if (span == 0) span <- max(abs(bounds), 1) * 0.2
  xlim <- bounds + c(-0.12, 0.12) * span
  ticks <- pretty(xlim, n = 5)
  ticks <- ticks[ticks >= xlim[1] & ticks <= xlim[2]]
  
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(
    bg = "white", fg = ink, family = "sans",
    mar = c(5.5, 9.5, 6.5, 2),
    oma = c(2.5, 0, 0, 0),
    las = 1, xaxs = "i", yaxs = "i"
  )
  
  graphics::plot(
    NA, xlim = xlim, ylim = c(0.5, 3.5),
    axes = FALSE, xlab = "", ylab = "", bty = "n"
  )
  
  # Subtle horizontal bands, with emphasis on the original prior.
  for (i in seq_along(y)) {
    graphics::rect(
      xlim[1], y[i] - 0.43, xlim[2], y[i] + 0.43,
      col = if (i == 2L) "#EDF1FA" else "#F7F9FC", border = NA
    )
  }
  
  for (tick in ticks) {
    graphics::segments(tick, 0.57, tick, 3.43, col = "#E2E8F0", lwd = 1)
  }
  graphics::segments(0, 0.57, 0, 3.43, col = "#8795A7", lty = 2, lwd = 1.4)
  
  for (i in seq_along(y)) {
    # A soft halo and crisp central line show the supplied interval.
    graphics::segments(
      lower[i], y[i], upper[i], y[i],
      col = grDevices::adjustcolor(colors[i], alpha.f = 0.15),
      lwd = 14, lend = "round"
    )
    graphics::segments(
      lower[i], y[i], upper[i], y[i],
      col = colors[i], lwd = 3, lend = "round"
    )
    graphics::segments(
      c(lower[i], upper[i]), y[i] - 0.07,
      c(lower[i], upper[i]), y[i] + 0.07,
      col = colors[i], lwd = 2
    )
    graphics::points(means[i], y[i], pch = 21, bg = colors[i], col = "white", cex = 2, lwd = 2)
    
    # Numbers are offset below each interval to avoid covering the estimate.
    estimate_text <- sprintf(
      paste0("Mean %.", digits, "f   |   [%.", digits, "f, %.", digits, "f]"),
      means[i], lower[i], upper[i]
    )
    graphics::text(
      mean(xlim), y[i] - 0.27, labels = estimate_text,
      cex = 0.78, col = muted
    )
    
    # Labels in the left margin, aligned to their respective rows.
    graphics::mtext(
      labels[i], side = 2, at = y[i] + 0.09,
      line = 0.9, adj = 1, las = 1, font = 2,
      cex = 0.95, col = colors[i]
    )
    graphics::mtext(
      scales[i], side = 2, at = y[i] - 0.14,
      line = 0.9, adj = 1, las = 1,
      cex = 0.76, col = muted
    )
  }
  
  graphics::axis(
    1, at = ticks, col = NA, col.ticks = "#BCC7D5",
    col.axis = muted, cex.axis = 0.85, tck = -0.015
  )
  graphics::mtext(
    "Temporal coefficient", side = 1, line = 2.6,
    cex = 0.95, col = ink
  )
  graphics::mtext(
    title, side = 3, line = 4, adj = 0,
    cex = 1.15, font = 2, col = ink
  )
  graphics::mtext(
    subtitle, side = 3, line = 2.5, adj = 0,
    cex = 0.86, col = muted
  )
  graphics::mtext(
    paste0("Dots: posterior means   |   Lines: ", interval_label,
           "   |   Dashed line: zero"),
    side = 1, outer = TRUE, line = 0.8,
    cex = 0.76, col = muted
  )
  
  invisible(data.frame(
    prior = labels, scale = c(0.35, 0.50, 0.65),
    mean = means, ci_lower = lower, ci_upper = upper
  ))
}

# No top-level plotting call: call the function from your analysis script.

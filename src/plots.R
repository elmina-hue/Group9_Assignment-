# plots.R
# Source this file, then call plot_temporal_network().
# Requires the igraph package. Sourcing does not draw a plot.

plot_temporal_network <- function(
  network,
  level = c("group", "individual"),
  participant = NULL,
  model = "VAR",
  lag = 1,
  node_names = NULL,
  labels = NULL,
  symptom_colors = c(
    down = "#DCE8F5",
    tired = "#FAE6D3",
    sleep = "#D8EEE8",
    nervous = "#E7DFF4",
    worry = "#F4DDE5",
    irritable = "#E8EDCF",
    anhedonia = "#E7DFF4",
    worthlessness = "#F4DDE5",
    concentration = "#E8EDCF",
    appetite = "#FAE6D3",
    psychomotor = "#D8EEE8",
    suicidal = "#F2D9DD"
  ),
  title = NULL,
  subtitle = NULL,
  extractor = NULL,
  transpose = FALSE,
  threshold = 0.10,
  show_loops = TRUE,
  show_edge_labels = FALSE
) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Install igraph first: install.packages('igraph')")
  }

  level <- match.arg(level)

  # For fitted objects, supply the appropriate extraction function.
  # Otherwise, network should already be a coefficient matrix.
  if (!is.null(extractor)) {
    if (!is.function(extractor)) stop("extractor must be a function.")
    network <- extractor(network)
  }
  network <- as.matrix(network)

  if (!is.numeric(network) || nrow(network) != ncol(network) ||
      nrow(network) < 1L) {
    stop("network must be a non-empty, square numeric matrix.")
  }
  if (any(!is.finite(network))) {
    stop("network must contain only finite values (no NA, NaN, or Inf).")
  }
  if (!is.numeric(threshold) || length(threshold) != 1L ||
      !is.finite(threshold) || threshold < 0) {
    stop("threshold must be one finite, non-negative number.")
  }
  if (!is.logical(show_loops) || length(show_loops) != 1L || is.na(show_loops)) {
    stop("show_loops must be TRUE or FALSE.")
  }
  if (!is.logical(show_edge_labels) || length(show_edge_labels) != 1L ||
      is.na(show_edge_labels)) {
    stop("show_edge_labels must be TRUE or FALSE.")
  }

  # Rows are lagged predictors and columns are outcomes in bvarnet's matrix.
  # Transpose only when a different extractor returns the reverse convention.
  if (transpose) network <- t(network)

  if (is.null(node_names)) node_names <- colnames(network)
  if (is.null(node_names)) node_names <- rownames(network)
  if (is.null(node_names)) {
    node_names <- paste("Node", seq_len(ncol(network)))
  }
  if (length(node_names) != ncol(network)) {
    stop("node_names must have one entry per node, in matrix order.")
  }
  if (is.null(labels)) labels <- node_names
  if (length(labels) != ncol(network)) {
    stop("labels must have one entry per node.")
  }
  labels <- as.character(labels)
  # Break long project symptom names at readable syllable boundaries.
  syllable_breaks <- c(
    concentration = "concen-\ntration",
    worthlessness = "worth-\nlessness",
    psychomotor = "psycho-\nmotor",
    anhedonia = "anhe-\ndonia"
  )
  split_labels <- unname(syllable_breaks[tolower(labels)])
  needs_break <- !is.na(split_labels) & !grepl("\\n", labels)
  labels[needs_break] <- split_labels[needs_break]

  node_colors <- unname(symptom_colors[node_names])
  node_colors[is.na(node_colors)] <- "#E2E8F0"

  if (is.null(title)) {
    title <- if (level == "group") {
      "Group-level temporal symptom network"
    } else {
      "Individual-level temporal symptom network"
    }
  }
  if (is.null(subtitle)) {
    scope <- if (level == "group") {
      "Group level"
    } else if (!is.null(participant)) {
      paste("Participant", participant)
    } else {
      "Individual level"
    }
    subtitle <- paste(scope, model, paste("Lag", lag), sep = " | ")
  }

  # Hide coefficients smaller than the selected threshold; optionally hide loops.
  shown <- network
  shown[abs(shown) < threshold] <- 0
  if (!show_loops) diag(shown) <- 0

  graph <- igraph::graph_from_adjacency_matrix(
    shown, mode = "directed", weighted = TRUE, diag = TRUE
  )
  weights <- igraph::E(graph)$weight
  edge_base_colors <- ifelse(weights >= 0, "#0D9488", "#E11D48")
  # Keep weaker edges lighter, as in the earlier qgraph version with fade = TRUE.
  max_weight <- max(abs(network))
  edge_alpha <- if (max_weight > 0) {
    pmax(0.08, abs(weights) / max_weight)
  } else {
    rep(1, length(weights))
  }
  edge_rgb <- grDevices::col2rgb(edge_base_colors) / 255
  edge_colors <- grDevices::rgb(
    edge_rgb[1, ], edge_rgb[2, ], edge_rgb[3, ], alpha = edge_alpha
  )
  edge_labels <- if (show_edge_labels) sprintf("%.2f", weights) else NA

  # Match the app's fixed circular layout and edge styling.
  n_nodes <- ncol(network)
  theta <- pi / 2 - (seq_len(n_nodes) - 1) * 2 * pi / n_nodes
  # Scale the layout with the limits to preserve its size, with extra room
  # around the loops and a slight downward shift.
  layout <- 1.29 * cbind(cos(theta), sin(theta))
  layout[, 2] <- layout[, 2] - 0.4

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  # Leave room below the standalone plot for its explanatory notes.
  graphics::par(
    bg = "white", fg = "#20334A", family = "sans",
    oma = c(4.5, 0, 0, 0),
    xpd = TRUE
  )
  plot(
    graph,
    layout = layout,
    rescale = FALSE,
    xlim = c(-1.55, 1.55),
    ylim = c(-1.55, 1.55),
    main = title,
    vertex.label = labels,
    vertex.size = 54,
    vertex.color = node_colors,
    vertex.frame.color = "#7D8FA2",
    vertex.label.color = "#20334A",
    vertex.label.family = "sans",
    vertex.label.cex = 0.8,
    edge.color = edge_colors,
    edge.width = pmin(4, 0.8 + abs(weights) * 4),
    edge.arrow.size = 0.55,
    edge.curved = 0.15,
    edge.label = edge_labels,
    edge.label.color = "#65768A",
    edge.label.cex = 0.8,
    margin = 0.1
  )

  graphics::mtext(
    subtitle, side = 3, line = 0.5,
    cex = 0.9, col = "#65768A"
  )

  graphics::mtext(
    "Teal: positive coefficient   |   Rose: negative coefficient",
    side = 1, outer = TRUE, line = 0.7,
    cex = 0.8, col = "#65768A"
  )
  graphics::mtext(
    "Thicker edges: larger absolute weights | Loops: autoregressive effects",
    side = 1, outer = TRUE, line = 2.0,
    cex = 0.75, col = "#65768A"
  )
  graphics::mtext(
    "Edge appearance does not indicate statistical evidence or uncertainty.",
    side = 1, outer = TRUE, line = 3.3,
    cex = 0.7, col = "#65768A"
  )

  invisible(graph)
}

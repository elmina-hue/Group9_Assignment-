# plots.R
# Source this file, then call plot_temporal_network().
# Requires the qgraph package. Sourcing does not draw a plot.

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
    irritable = "#E8EDCF"
  ),
  title = NULL,
  subtitle = NULL,
  extractor = NULL,
  transpose = FALSE,
  maximum = NULL,
  ...   #extra qgraph options; these can override the plotting defaults
) {
  if (!requireNamespace("qgraph", quietly = TRUE)) {
    stop("Install qgraph first: install.packages('qgraph')")
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

  # qgraph interprets entry [i, j] as an edge from node i to node j.
  # bvarnet's extractor returns rows as lagged predictors and columns as
  # outcomes, so use transpose = TRUE for its output.
  if (transpose) network <- t(network)

  if (is.null(node_names)) node_names <- colnames(network)
  if (is.null(node_names)) node_names <- rownames(network)
  if (is.null(node_names)) {
    node_names <- paste("Node", seq_len(ncol(network)))
  }
  if (length(node_names) != ncol(network)) {
    stop("node_names must have one entry per node, in matrix order.")
  }

  if (is.null(labels)) {
    labels <- tools::toTitleCase(gsub("_", " ", node_names))
  }
  if (length(labels) != ncol(network)) {
    stop("labels must have one entry per node.")
  }

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

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(
    bg = "white", fg = "#20334A", family = "sans",
    oma = c(6, 0, 5, 0)
  )

  plot_args <- list(
    input = network,
    layout = "circle",
    layoutScale = c(0.90, 0.95),
    labels = labels,
    shape = "circle",
    color = node_colors,
    borders = TRUE,
    border.color = "#7D8FA2",
    border.width = 1.5,
    vsize = 10,
    label.cex = 1.15,
    label.scale = FALSE,
    label.color = "#20334A",
    directed = TRUE,
    # Draw diagonal coefficients as self-loops (autoregressive effects).
    diag = TRUE,
    posCol = "#148A87",
    negCol = "#C25773",
    fade = TRUE,
    minimum = 0,
    esize = 7,
    asize = 4,
    curve = 0.2,
    edge.labels = FALSE,
    mar = c(4, 3, 4, 3)
  )

  # Set the same maximum in each plot to compare edge strengths on one scale.
  if (!is.null(maximum)) plot_args$maximum <- maximum

  # Named extra qgraph arguments override the defaults above, for example
  # directed = FALSE or layout = "spring".
  extra_args <- list(...)
  if (length(extra_args)) {
    if (is.null(names(extra_args)) ||
        any(is.na(names(extra_args)) | names(extra_args) == "")) {
      stop("All extra qgraph arguments in ... must be named.")
    }
    plot_args[names(extra_args)] <- extra_args
  }

  graph <- do.call(qgraph::qgraph, plot_args)

  graphics::mtext(
    title, side = 3, outer = TRUE, line = 2.5,
    cex = 1.15, font = 2, col = "#20334A"
  )
  graphics::mtext(
    subtitle, side = 3, outer = TRUE, line = 0.9,
    cex = 0.9, col = "#65768A"
  )
  graphics::mtext(
    "Teal: positive coefficient   |   Rose: negative coefficient",
    side = 1, outer = TRUE, line = 1.5,
    cex = 0.8, col = "#65768A"
  )
  graphics::mtext(
    "Thicker edges: larger absolute weights | Loops: autoregressive effects",
    side = 1, outer = TRUE, line = 2.9,
    cex = 0.75, col = "#65768A"
  )
  graphics::mtext(
    "Edge appearance does not indicate statistical evidence or uncertainty.",
    side = 1, outer = TRUE, line = 4.3,
    cex = 0.7, col = "#65768A"
  )
  invisible(graph)
}

# Description of app for user:

# This Shiny app compares a group-level symptom network with a participant-level
# network. It uses the shared fitting functions in analysis_network.R and fits
# the group model at startup so it can use posterior draws to calculate intervals
# for participant-minus-group edge differences. The group posterior mean network
# is also saved to data/group_network.csv. Each arrow represents a lag-one
# temporal association: a symptom on one day predicting another symptom the next
# day. The example data are simulated; these associations do not prove causation.

# Brief description of script:

# Two-page Shiny app for participant data input and network comparison.
# Run from the project root with:
#   shiny::runApp("Group9_Assignment--main/src/symptom_network_app.R")
# The group model uses the `ema` data frame created by simulate_data.R. Its
# posterior mean matrix is written to src/data/group_network.csv for reference;
# the app refits the model at startup because edge-difference intervals require
# posterior draws, which are not stored in that CSV.

required_packages <- c("shiny", "bvarnet", "igraph", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                                logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Install required packages first: ", paste(missing_packages, collapse = ", "))

symptoms <- c("down", "tired", "sleep", "anhedonia", "worthlessness",
              "concentration", "appetite", "psychomotor", "suicidal")

# Locate this app's supporting scripts and data relative to the app file.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
source_files <- vapply(sys.frames(), function(frame) if (!is.null(frame$ofile)) as.character(frame$ofile)[1] else "", character(1))
source_files <- source_files[nzchar(source_files)]
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else if (length(source_files)) tail(source_files, 1) else "symptom_network_app.R"
app_dir <- dirname(normalizePath(script_path, mustWork = FALSE))
simulation_env <- new.env(parent = globalenv())
source(file.path(app_dir, "simulate_data.R"), local = simulation_env)
source(file.path(app_dir, "plots.R"), local = TRUE)
analysis_env <- new.env(parent = asNamespace("bvarnet"))
source(file.path(app_dir, "analysis_network.R"), local = analysis_env)
group_data <- simulation_env$ema

# Path for the saved group posterior mean network matrix.
group_network_path <- file.path(app_dir, "data", "group_network.csv")

# Validate uploaded participant data: one row per consecutive day, with all
# nine symptom ratings recorded as integers from 0 to 3.
validate_symptom_data <- function(data) {
  required <- c("day", symptoms)
  missing <- setdiff(required, names(data))
  if (length(missing)) stop("Missing CSV columns: ", paste(missing, collapse = ", "))
  if (!all(vapply(data[symptoms], is.numeric, logical(1)))) stop("Symptom columns must be numeric ratings.")
  if (!is.numeric(data$day) || any(!is.finite(data$day)) || anyNA(data[c("day", symptoms)])) stop("Day and symptom columns cannot contain missing or non-finite values.")
  if (anyDuplicated(data$day)) stop("Each day must appear only once.")
  data <- data[order(data$day), , drop = FALSE]
  if (nrow(data) < 30) stop("At least 30 consecutive daily observations are required.")
  if (any(diff(data$day) != 1)) stop("Days must be consecutive.")
  ratings <- as.matrix(data[symptoms])
  if (any(ratings < 0 | ratings > 3)) stop("PHQ-9 ratings must be between 0 and 3.")
  if (any(ratings != round(ratings))) stop("PHQ-9 ratings must be whole numbers from 0 to 3.")
  data
}

# The shared analysis functions call bvarnet::bvar(). If its Stan models have
# not been compiled yet, compile them once and retry that shared analysis call.
run_with_stan_setup <- function(analysis_call) {
  tryCatch(
    analysis_call(),
    error = function(e) {
      if (!grepl("Stan models are not set up", conditionMessage(e), fixed = TRUE)) {
        stop(e)
      }
      bvarnet::bvarnet_setup_models(method = "compile", ask = FALSE)
      analysis_call()
    }
  )
}

# The sourced analysis functions resolve bvarnet functions in their enclosing
# environment. Use the same priors as the project's analysis specification.
network_priors <- bvarnet::set_priors(
  phi = bvarnet::prior(family = "cauchy", loc = 0, scale = 0.50),
  sd_u = bvarnet::prior(family = "normal", loc = 0, scale = 0.35)
)

# Fit the group model through the shared analysis function so the app and report
# use the same specification. Retain its fit for posterior edge comparisons.
prepare_group_reference <- function() {
  analysis <- run_with_stan_setup(function() {
    analysis_env$run_group_network(
      ema = as.data.frame(group_data), ys = symptoms, pr = network_priors,
      iter = 4000, warmup = 1000, chains = 4, cores = 1, seed = 2026
    )
  })
  fit <- analysis$fit_group
  net <- bvarnet::extract_network_matrix(fit, stat = "mean")
  net <- net[symptoms, symptoms, drop = FALSE]
  dir.create(dirname(group_network_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(net, group_network_path, row.names = TRUE)
  list(network = net, fit = fit)
}

# Fit the selected participant with the shared individual-level analysis
# function. Return both the posterior mean network and the fitted model.
fit_person_network <- function(d) {
  d <- validate_symptom_data(d)
  participant_data <- transform(d, id = 1)
  analysis <- run_with_stan_setup(function() {
    analysis_env$run_individual_network(
      ema = participant_data, ys = symptoms, pr = network_priors,
      participant_id = 1, iter = 4000, warmup = 1000, chains = 4,
      cores = 1, seed = 2026
    )
  })
  fit <- analysis$fit_individual
  m <- bvarnet::extract_network_matrix(fit, stat = "mean")
  list(network = m[symptoms, symptoms, drop = FALSE], fit = fit)
}
# bvarnet network matrices use rows for lagged predictors and columns for outcomes.

# Generate example data using one simulated participant from the project data.
make_example_data <- function(person_id, n_days) {
  d <- group_data[group_data$id == person_id, c("day", symptoms), drop = FALSE]
  utils::head(d, n_days)
}
#returns one participant's simulated rows with a day column and nine symptom columns

# Build the group reference once when the app starts. The mean matrix is used
# for plotting; posterior draws are used to summarize edge differences.
group_analysis <- prepare_group_reference()
group_network <- group_analysis$network
group_draws <- analysis_env$get_draws_p1(group_analysis$fit, parameter = "phi")

ui <- shiny::fluidPage(
  #sets up the page and title
  shiny::tags$head(shiny::tags$style(shiny::HTML("body{background:#f5f7fb;color:#0f172a}.container-fluid{max-width:1250px}.well{background:white}"))),
  shiny::titlePanel("Participant and group symptom networks"),
  shiny::p("Explore how symptoms predict one another from one day to the next. Temporal associations do not establish causation."),
  shiny::tabsetPanel(id = "pages",
    #defines page one, participant data:
    shiny::tabPanel("1. Participant data", value = "input",
      shiny::radioButtons("source", "Choose participant data", inline = TRUE,
        choices = c("Use a simulated participant" = "example", "Upload a participant CSV" = "upload")),
      shiny::conditionalPanel("input.source == 'example'",
        shiny::selectInput("person_id", "Simulated participant", choices = seq_len(simulation_env$N), selected = 1),
        shiny::numericInput("days", "Daily observations", 100, min = 30, max = simulation_env$Tt),
        shiny::actionButton("generate", "Generate and analyse")),
      shiny::conditionalPanel("input.source == 'upload'",
        shiny::fileInput("file", "Participant CSV", accept = ".csv"),
        shiny::helpText("Required columns: day plus the nine PHQ-9 symptom domains (down, tired, sleep, anhedonia, worthlessness, concentration, appetite, psychomotor, suicidal). One row per consecutive day; ratings 0–3."),
        shiny::actionButton("analyse", "Analyse uploaded data")),
      shiny::verbatimTextOutput("status"), shiny::h4("Daily ratings (first 10 rows)"), shiny::tableOutput("preview")),
    #defines page two, comparing nerworks:
    shiny::tabPanel("2. Compare networks", value = "networks",
      shiny::textOutput("counts"),
      shiny::fluidRow(shiny::column(6, shiny::plotOutput("group_plot", height = "500px")),
                      shiny::column(6, shiny::plotOutput("person_plot", height = "500px"))),
      shiny::p("Teal: positive | Rose: negative | Arrow: yesterday’s predictor → today’s outcome | Self-loop: carry-over"),
      shiny::sliderInput("threshold", "Minimum displayed absolute coefficient", min = 0, max = .5, value = 0, step = .01),
      shiny::checkboxInput("loops", "Show self-loops", TRUE), shiny::checkboxInput("labels", "Show edge values", TRUE),
      shiny::h4("Edge estimates"), shiny::tableOutput("edges")))
)

# Connect the existing interface to the analysis. Model errors are returned to
# the status message so they do not stop the Shiny session.
server <- function(input, output, session) {
  result <- shiny::eventReactive(input$generate + input$analyse, {  #eventReactive runs the analysis
    tryCatch({  #tryCatch captures errors so the app can show them as status message
      d <- if (identical(input$source, "example")) make_example_data(as.integer(input$person_id), input$days) else {
        shiny::req(input$file); readr::read_csv(input$file$datapath, show_col_types = FALSE)
      }
      d <- validate_symptom_data(d)
      person <- fit_person_network(d)
      list(data = d, individual = person$network, fit = person$fit, error = NULL)
    }, error = function(e) list(data = NULL, individual = NULL, error = conditionMessage(e)))
  }, ignoreInit = TRUE)
  #filling in the interface:
  output$status <- shiny::renderText({ x <- result(); if (is.null(x$error)) "Model fitted using the project’s Bayesian VAR approach." else paste("Could not analyse data:", x$error) })
  output$preview <- shiny::renderTable({ x <- result(); shiny::req(is.null(x$error)); utils::head(x$data, 10) })
  output$counts <- shiny::renderText({ x <- result(); shiny::req(is.null(x$error)); sprintf("%d daily observations | %d consecutive transitions", nrow(x$data), nrow(x$data) - 1) })
  output$group_plot <- shiny::renderPlot({
    plot_temporal_network(
      group_network, level = "group", title = "Group network", subtitle = "",
      node_names = symptoms, threshold = input$threshold,
      show_loops = input$loops, show_edge_labels = input$labels,
      show_notes = FALSE
    )
  })
  output$person_plot <- shiny::renderPlot({
    x <- result(); shiny::req(is.null(x$error))
    plot_temporal_network(
      x$individual, level = "individual", title = "Participant network", subtitle = "",
      node_names = symptoms, threshold = input$threshold,
      show_loops = input$loops, show_edge_labels = input$labels,
      show_notes = FALSE
    )
  })
  output$edges <- shiny::renderTable({  #one row per directed symptom pair
    x <- result(); shiny::req(is.null(x$error))
    person_draws <- analysis_env$get_draws_p1(x$fit, parameter = "phi")
    n_draws <- min(nrow(group_draws), nrow(person_draws))
    if (n_draws < 1L) stop("The fitted models did not return posterior draws.")
    # The models retain different numbers of draws, so take equally sized sets
    # before calculating independent draw-wise participant-minus-group deltas.
    group_idx <- unique(round(seq(1, nrow(group_draws), length.out = n_draws)))
    person_idx <- unique(round(seq(1, nrow(person_draws), length.out = n_draws)))
    n_draws <- min(length(group_idx), length(person_idx))
    group_idx <- group_idx[seq_len(n_draws)]
    person_idx <- person_idx[seq_len(n_draws)]
    rows <- vector("list", length(symptoms)^2)
    k <- 0L
    for (i in seq_along(symptoms)) for (j in seq_along(symptoms)) {
      k <- k + 1L
      edge <- sprintf("phi[%d,%d]", i, j)
      group_edge <- analysis_env$get_p1_sleep_tired(group_draws[group_idx, , drop = FALSE], edge = edge)
      person_edge <- analysis_env$get_p1_sleep_tired(person_draws[person_idx, , drop = FALSE], edge = edge)
      delta <- analysis_env$get_delta(person_edge, group_edge)
      difference <- analysis_env$get_posterior_difference(delta)
      rows[[k]] <- data.frame(
        `Previous-day symptom` = symptoms[i], Outcome = symptoms[j],
        `Group estimate` = group_network[i, j],
        `Participant estimate` = x$individual[i, j],
        Difference = x$individual[i, j] - group_network[i, j],
        `90% CI lower` = difference[["lower"]],
        `90% CI upper` = difference[["upper"]],
        check.names = FALSE
      )
    }
    z <- do.call(rbind, rows)
    z[order(-abs(z$Difference)), ]
  }, digits = 3)
}

# Starting the app:
app <- shiny::shinyApp(ui, server)
if (sys.nframe() == 0L || interactive()) shiny::runApp(app)

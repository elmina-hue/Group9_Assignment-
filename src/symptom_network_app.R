# Participant and group symptom network app
#
# Use this app to explore how yesterday's symptom ratings are associated with
# today's ratings. Choose a simulated participant or upload a CSV, preview the
# daily ratings, and then fit that participant's model. The app compares the
# participant network with a group reference network and shows posterior
# intervals for the difference between their edge estimates. Bold intervals
# do not include zero. These temporal associations describe prediction, not
# causation; the built-in participant data are simulated for demonstration.
#
# App details:
# - Uploaded data need one row per consecutive day, a `day` column, and all
#   nine symptom columns shown in the upload instructions. Ratings must be
#   whole numbers from 0 to 3, with at least 30 days of observations.
# - The participant model is fitted each time the user selects Analyse. This
#   can take a few minutes; the daily ratings can be previewed beforehand.
# - The group network is loaded from the saved posterior cache when available.
#   If the cache is missing, the app fits the group model once and saves the
#   network and posterior draws for later launches.
#
# Run from the repository root with:
#   shiny::runApp("Group9_Assignment--main/src/symptom_network_app.R")

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
group_posterior_cache_path <- file.path(app_dir, "data", "group_posterior_cache.rds")

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

# Load the compact posterior cache when available. Otherwise fit through the
# shared analysis function, then cache the mean matrix and phi draws needed here.
prepare_group_reference <- function() {
  if (file.exists(group_posterior_cache_path)) {
    cached <- readRDS(group_posterior_cache_path)
    if (!is.list(cached) || !all(c("network", "draws") %in% names(cached))) {
      stop("Group posterior cache must contain 'network' and 'draws'.")
    }
    net <- as.matrix(cached$network)
    draws <- as.matrix(cached$draws)
    required_edges <- as.vector(outer(
      seq_along(symptoms), seq_along(symptoms),
      function(i, j) sprintf("phi[%d,%d]", i, j)
    ))
    if (!all(symptoms %in% rownames(net)) || !all(symptoms %in% colnames(net))) {
      stop("Cached group network must label all nine symptoms on rows and columns.")
    }
    if (!all(required_edges %in% colnames(draws)) || nrow(draws) < 1L) {
      stop("Group posterior cache is missing temporal-coefficient draws.")
    }
    net <- net[symptoms, symptoms, drop = FALSE]
    draws <- draws[, required_edges, drop = FALSE]
    return(list(network = net, draws = draws, fit = NULL))
  }

  analysis <- run_with_stan_setup(function() {
    analysis_env$run_group_network(
      ema = as.data.frame(group_data), ys = symptoms, pr = network_priors,
      iter = 4000, warmup = 1000, chains = 4, cores = 1, seed = 2026
    )
  })
  fit <- analysis$fit_group
  net <- bvarnet::extract_network_matrix(fit, stat = "mean")
  net <- net[symptoms, symptoms, drop = FALSE]
  draws <- as.matrix(analysis_env$get_draws_p1(fit, parameter = "phi"))
  dir.create(dirname(group_network_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(network = net, draws = draws), group_posterior_cache_path)
  utils::write.csv(net, group_network_path, row.names = TRUE)
  list(network = net, draws = draws, fit = fit)
}

# Fit a fresh individual-level model for the selected participant.
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
  net <- m[symptoms, symptoms, drop = FALSE]
  draws <- as.matrix(analysis_env$get_draws_p1(fit, parameter = "phi"))
  list(network = net, draws = draws)
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
group_draws <- group_analysis$draws

ui <- shiny::fluidPage(
  # Page title and short explanation of the network comparison for app users.
  shiny::tags$head(shiny::tags$style(shiny::HTML("body{background:#f5f7fb;color:#0f172a}.container-fluid{max-width:1250px}.well{background:white}"))),
  shiny::titlePanel("Participant and group symptom networks"),
  shiny::p("Compare a participant's symptom network with the group network. Arrows show how one day's symptom ratings are associated with the next day's ratings; they do not establish cause and effect."),
  shiny::tabsetPanel(id = "pages",
    # Page 1: select, upload, and preview the participant's daily ratings.
    shiny::tabPanel("1. Participant data", value = "input",
      shiny::radioButtons("source", "Choose participant data", inline = TRUE,
        choices = c("Use a simulated participant" = "example", "Upload a participant CSV" = "upload")),
      shiny::conditionalPanel("input.source == 'example'",
        shiny::selectInput("person_id", "Simulated participant", choices = seq_len(simulation_env$N), selected = 1),
        shiny::numericInput("days", "Daily observations", 100, min = 30, max = simulation_env$Tt),
        shiny::helpText("Generate the selected participant's ratings first. Review the preview, then choose Analyse participant to fit the model."),
        shiny::actionButton("generate", "Generate participant data"),
        shiny::actionButton("analyse_example", "Analyse participant")),
      shiny::conditionalPanel("input.source == 'upload'",
        shiny::fileInput("file", "Participant CSV", accept = ".csv"),
        shiny::helpText("Required columns: day plus the nine PHQ-9 symptom domains (down, tired, sleep, anhedonia, worthlessness, concentration, appetite, psychomotor, suicidal). One row per consecutive day; ratings 0–3."),
        shiny::helpText("Preview the ratings before selecting Analyse uploaded data. Fitting the participant model may take a few minutes."),
        shiny::actionButton("preview_upload", "Preview uploaded data"),
        shiny::actionButton("analyse_upload", "Analyse uploaded data")),
      shiny::verbatimTextOutput("status"), shiny::h4("Daily ratings (first 10 rows)"), shiny::tableOutput("preview")),
    # Page 2: compare network plots and inspect the estimated edge differences.
    shiny::tabPanel("2. Compare networks", value = "networks",
      shiny::textOutput("counts"),
      shiny::fluidRow(shiny::column(6, shiny::plotOutput("group_plot", height = "500px")),
                      shiny::column(6, shiny::plotOutput("person_plot", height = "500px"))),
      shiny::p("Teal: positive | Rose: negative | Arrow: yesterday’s predictor → today’s outcome | Self-loop: carry-over"),
      shiny::sliderInput("threshold", "Minimum displayed absolute coefficient", min = 0, max = .5, value = 0, step = .01),
      shiny::checkboxInput("loops", "Show self-loops", TRUE), shiny::checkboxInput("labels", "Show edge values", TRUE),
      shiny::h4("Edge estimates"), shiny::uiOutput("edges"),
      shiny::p("Bold 90% credible intervals do not include zero, indicating a difference between the participant and group edge estimates.")))
)

# Connect the existing interface to the analysis. Model errors are returned to
# the status message so they do not stop the Shiny session.
server <- function(input, output, session) {
  participant_data <- shiny::reactiveVal(NULL)
  model_result <- shiny::reactiveVal(NULL)
  status_message <- shiny::reactiveVal("Generate or preview participant data to begin.")

  shiny::observeEvent(input$source, {
    participant_data(NULL)
    model_result(NULL)
    status_message("Generate or preview participant data to begin.")
  }, ignoreInit = TRUE)

  shiny::observeEvent(input$generate, {
    tryCatch({
      d <- make_example_data(as.integer(input$person_id), input$days)
      participant_data(validate_symptom_data(d))
      model_result(NULL)
      status_message("Ratings are ready. Click Analyse participant to fit the model.")
    }, error = function(e) {
      participant_data(NULL)
      status_message(paste("Could not generate data:", conditionMessage(e)))
    })
  })

  shiny::observeEvent(input$preview_upload, {
    tryCatch({
      shiny::req(input$file)
      d <- readr::read_csv(input$file$datapath, show_col_types = FALSE)
      participant_data(validate_symptom_data(d))
      model_result(NULL)
      status_message("Ratings are ready. Click Analyse uploaded data to fit the model.")
    }, error = function(e) {
      participant_data(NULL)
      status_message(paste("Could not load data:", conditionMessage(e)))
    })
  })

  shiny::observeEvent(input$analyse_example + input$analyse_upload, {
    d <- participant_data()
    if (is.null(d)) {
      status_message("Generate or preview participant data before analysis.")
      return()
    }
    tryCatch({
      status_message("Fitting individual model…")
      person <- shiny::withProgress(
        message = "Fitting individual model…",
        detail = "Sampling the participant-level VAR. This may take a few minutes.",
        value = 0,
        {
          fit_person_network(d)
        }
      )
      model_result(list(data = d, individual = person$network,
                        draws = person$draws, error = NULL))
      status_message("Analysis complete. Showing the network comparison.")
      shiny::updateTabsetPanel(session, "pages", selected = "networks")
    }, error = function(e) {
      model_result(list(error = conditionMessage(e)))
      status_message(paste("Could not analyse data:", conditionMessage(e)))
    })
  })

  output$status <- shiny::renderText(status_message())
  output$preview <- shiny::renderTable({
    d <- participant_data()
    shiny::req(!is.null(d))
    utils::head(d, 10)
  })
  output$counts <- shiny::renderText({
    d <- participant_data()
    shiny::req(!is.null(d))
    sprintf("%d daily observations | %d consecutive transitions", nrow(d), nrow(d) - 1)
  })
  output$group_plot <- shiny::renderPlot({
    plot_temporal_network(
      group_network, level = "group", title = "Group network", subtitle = "",
      node_names = symptoms, threshold = input$threshold,
      show_loops = input$loops, show_edge_labels = input$labels,
      show_notes = FALSE
    )
  })
  output$person_plot <- shiny::renderPlot({
    x <- model_result(); shiny::req(!is.null(x), is.null(x$error))
    plot_temporal_network(
      x$individual, level = "individual", title = "Participant network", subtitle = "",
      node_names = symptoms, threshold = input$threshold,
      show_loops = input$loops, show_edge_labels = input$labels,
      show_notes = FALSE
    )
  })
  output$edges <- shiny::renderUI({  #one row per directed symptom pair
    x <- model_result(); shiny::req(!is.null(x), is.null(x$error))
    person_draws <- x$draws
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
    z <- z[order(-abs(z$Difference)), , drop = FALSE]
    excludes_zero <- z[["90% CI lower"]] > 0 | z[["90% CI upper"]] < 0
    interval_columns <- match(c("90% CI lower", "90% CI upper"), names(z))
    body_rows <- lapply(seq_len(nrow(z)), function(i) {
      cells <- lapply(seq_along(z), function(j) {
        value <- z[[j]][i]
        formatted <- if (is.numeric(z[[j]])) sprintf("%.3f", value) else as.character(value)
        if (j %in% interval_columns && excludes_zero[i]) {
          shiny::tags$td(shiny::tags$strong(formatted))
        } else {
          shiny::tags$td(formatted)
        }
      })
      do.call(shiny::tags$tr, cells)
    })
    shiny::tags$table(
      class = "table table-striped",
      shiny::tags$thead(shiny::tags$tr(lapply(names(z), shiny::tags$th))),
      do.call(shiny::tags$tbody, body_rows)
    )
  })
}

# Starting the app:
app <- shiny::shinyApp(ui, server)
if (sys.nframe() == 0L || interactive()) shiny::runApp(app)

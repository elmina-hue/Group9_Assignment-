# Description of app for user:

# The app loads a saved group network or fits a multilevel Bayesian VAR to the 
# group data if no saved network exists. It separately fits a Bayesian VAR to 
# one participant’s daily symptom ratings. Both produce lag-one 
# temporal networks, where an arrow means one symptom predicts another on the 
# next day. The app displays those networks side by side and provides edge 
# estimates for comparison. The example data are simulated, and the network 
# associations should not be interpreted as proof of causation.

# Brief description of script:

# Combined two-page Shiny app: participant data input and side-by-side networks.
# Run from the project root with: shiny::runApp("symptom_network_app.R")
# The group reference is fitted from data/group_symptoms.csv (or generated from
# src/simulate_data.R) and cached in data/group_network.csv.

required_packages <- c("shiny", "bvarnet", "igraph", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                                logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Install required packages first: ", paste(missing_packages, collapse = ", "))

symptoms <- c("down", "tired", "sleep", "nervous", "worry", "irritable")

#find the script's own location and locate project files
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
source_files <- vapply(sys.frames(), function(frame) if (!is.null(frame$ofile)) as.character(frame$ofile)[1] else "", character(1))
source_files <- source_files[nzchar(source_files)]
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else if (length(source_files)) tail(source_files, 1) else "symptom_network_app.R"
project_root <- dirname(normalizePath(script_path, mustWork = FALSE))
if (basename(project_root) != "group_project") project_root <- getwd()
source(file.path(project_root, "Group9_Assignment--main", "src", "simulate_data.R"), local = TRUE)

#set expected file locations
group_network_path <- file.path(project_root, "data", "group_network.csv") # group model's estimated network
group_data_path <- file.path(project_root, "data", "group_symptoms.csv") # observations for group model

#check participant data; requires day column and symptom columns
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
  if (any(as.matrix(data[symptoms]) < 1 | as.matrix(data[symptoms]) > 4)) stop("Ratings must be between 1 and 4.")
  data
}

# Preparing and fitting the group model:
#reads saved network and checks that it contains all symptoms as row and column names
read_network <- function(path) {
  x <- as.matrix(utils::read.csv(path, row.names = 1, check.names = FALSE))
  if (!all(symptoms %in% rownames(x)) || !all(symptoms %in% colnames(x))) stop("Group network must label all six symptoms on rows and columns.")
  x[symptoms, symptoms, drop = FALSE]
}

#tries to fit a model
fit_bvarnet <- function(...) tryCatch(bvarnet::bvar(...), error = function(e) {
  if (!grepl("Stan models are not set up", conditionMessage(e), fixed = TRUE)) stop(e)
  bvarnet::bvarnet_setup_models(method = "compile", ask = FALSE)
  bvarnet::bvar(...)
})

#gets the group network when the app starts; if the group network and symptom csv are both missing, it generates group data - we can adjust this
prepare_group_reference <- function() {
  if (file.exists(group_network_path)) return(read_network(group_network_path))
  d <- if (file.exists(group_data_path)) readr::read_csv(group_data_path, show_col_types = FALSE) else simulate_group_data()
  required <- c("id", "day", symptoms)
  if (!all(required %in% names(d))) stop("Group data needs columns: ", paste(required, collapse = ", "))
  #sets Cauchy prior for temporal coefficients & normal prior for residual standard deviations
  priors <- bvarnet::set_priors(phi = bvarnet::prior(family = "cauchy", loc = 0, scale = .5),
                                sd_u = bvarnet::prior(family = "normal", loc = 0, scale = .35))
  #fit lag-one Gaussian model; runs 4 chains with 2000 iterations per chain
  fit <- fit_bvarnet(id_col = "id", time_col = "day", y_cols = symptoms, re_cols = "Intercept",
                     re_temporal = TRUE, K = 1, data = as.data.frame(d), family = "gaussian",
                     priors = priors, iter = 2000, warmup = 500, chains = 4, cores = 1, seed = 2026)
  #extracts the posterior mean network and saves it
  net <- bvarnet::extract_network_matrix(fit, stat = "mean")
  dir.create(dirname(group_network_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(net, group_network_path, row.names = TRUE)
  net
}

# Fitting the individual participant model:
#checks participant data, fits a lag-one VAR, extracts posterior mean coefficiants
fit_person_network <- function(d) {
  d <- validate_symptom_data(d)
  fit <- fit_bvarnet(id_col = ".person", time_col = "day", y_cols = symptoms, re_cols = NULL,
                     re_temporal = FALSE, K = 1, data = transform(d, .person = 1), family = "gaussian",
                     iter = 1200, warmup = 400, chains = 2, cores = 1, seed = 2026)
  m <- bvarnet::extract_network_matrix(fit, stat = "mean")
  m[symptoms, symptoms, drop = FALSE]
}
# network matrix: rows->symptoms being predicted today; columns->symptoms used as a predictor yesterday.

#generates example data using the group network as the data-generating pattern
make_example_data <- function(group_matrix, n_days, seed) {
  set.seed(seed)
  phi <- group_matrix[symptoms, symptoms, drop = FALSE]
  if (max(Mod(eigen(phi)$values)) >= 1) stop("Group network is not stationary; example data cannot be generated safely.")
  y <- matrix(0, n_days, length(symptoms), dimnames = list(NULL, symptoms))
  #generates each days symptom values from the previous day's values plus random noise
  for (i in 2:n_days) y[i, ] <- phi %*% y[i - 1, ] + stats::rnorm(length(symptoms), 0, .7)
  ratings <- matrix(pmin(4, pmax(1, round(2.5 + y))), nrow = n_days,
                    ncol = length(symptoms), dimnames = list(NULL, symptoms))
  cbind(data.frame(day = seq_len(n_days)), as.data.frame(ratings, check.names = FALSE))
}
#results in data frame with day column and 6 symptom columns

# Drawing the networks:
#turns coefficient matrix into network diagram
draw_network <- function(A, title, threshold, loops, labels) {
  shown <- A  #A->coefficient matrix
  shown[abs(shown) < threshold] <- 0  #threshold = minimal abs coefficient to display
  if (!loops) diag(shown) <- 0  #loops - whether to show self-connections
  g <- igraph::graph_from_adjacency_matrix(t(shown), mode = "directed", weighted = TRUE, diag = TRUE) #converts matrix into graph
  theta <- pi / 2 - (0:5) * 2 * pi / 6  #calculates 6 evenly spaced points around a circle, one for each symptom
  weights <- igraph::E(g)$weight  #gets coefficient attached to each edge in the graph
  graphics::par(mar = c(1, 1, 3, 1))
  #draws the graph:
  plot(g, layout = cbind(cos(theta), sin(theta)), main = title, vertex.label = symptoms,
       vertex.size = 34, vertex.color = "#f1f5f9", vertex.frame.color = "#334155",
       vertex.label.color = "#0f172a", edge.color = ifelse(weights >= 0, "#0d9488", "#e11d48"),
       edge.width = pmin(7, 1 + abs(weights) * 8), edge.arrow.size = .55, edge.curved = .15,
       edge.label = if (labels) sprintf("%.2f", weights) else NA,
       edge.label.color = "#334155", edge.label.cex = .8, margin = .25)
}

# Building the two pages of the app:
#loads/fits the group network when the app starts
group_network <- prepare_group_reference()

ui <- shiny::fluidPage(
  #sets up the page and title
  shiny::tags$head(shiny::tags$style(shiny::HTML("body{background:#f5f7fb;color:#0f172a}.container-fluid{max-width:1250px}.well{background:white}"))),
  shiny::titlePanel("Participant and group symptom networks"),
  shiny::p("Explore how symptoms predict one another from one day to the next. Temporal associations do not establish causation."),
  shiny::tabsetPanel(id = "pages",
    #defines page one, participant data:
    shiny::tabPanel("1. Participant data", value = "input",
      shiny::radioButtons("source", "Choose participant data", inline = TRUE,
        choices = c("Generate an example participant" = "example", "Upload a participant CSV" = "upload")),
      shiny::conditionalPanel("input.source == 'example'",
        shiny::numericInput("days", "Daily observations", 60, min = 30, max = 365),
        shiny::numericInput("seed", "Random seed", 1, min = 1, max = 999999),
        shiny::actionButton("generate", "Generate and analyse")),
      shiny::conditionalPanel("input.source == 'upload'",
        shiny::fileInput("file", "Participant CSV", accept = ".csv"),
        shiny::helpText("Required columns: day, down, tired, sleep, nervous, worry, irritable. One row per consecutive day; ratings 1–4."),
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

# Connecting the interface to the analysis:
#tells the app what to do when someone interacts with it
server <- function(input, output, session) {
  result <- shiny::eventReactive(input$generate + input$analyse, {  #eventReactive runs the analysis
    tryCatch({  #tryCatch captures errors so the app can show them as status message
      d <- if (identical(input$source, "example")) make_example_data(group_network, input$days, input$seed) else {
        shiny::req(input$file); readr::read_csv(input$file$datapath, show_col_types = FALSE)
      }
      d <- validate_symptom_data(d)
      list(data = d, individual = fit_person_network(d), error = NULL)
    }, error = function(e) list(data = NULL, individual = NULL, error = conditionMessage(e)))
  }, ignoreInit = TRUE)
  #filling in the interface:
  output$status <- shiny::renderText({ x <- result(); if (is.null(x$error)) "Model fitted using the project’s Bayesian VAR approach." else paste("Could not analyse data:", x$error) })
  output$preview <- shiny::renderTable({ x <- result(); shiny::req(is.null(x$error)); utils::head(x$data, 10) })
  output$counts <- shiny::renderText({ x <- result(); shiny::req(is.null(x$error)); sprintf("%d daily observations | %d consecutive transitions", nrow(x$data), nrow(x$data) - 1) })
  output$group_plot <- shiny::renderPlot({ draw_network(group_network, "Group network", input$threshold, input$loops, input$labels) })
  output$person_plot <- shiny::renderPlot({ x <- result(); shiny::req(is.null(x$error)); draw_network(x$individual, "Participant network", input$threshold, input$loops, input$labels) })
  output$edges <- shiny::renderTable({  #makes 36-row table, one row for each directed pair of 6 symptoms
    x <- result(); shiny::req(is.null(x$error))
    g <- as.data.frame(as.table(group_network), stringsAsFactors = FALSE); names(g) <- c("Outcome", "Previous-day symptom", "Group estimate")
    p <- as.data.frame(as.table(x$individual), stringsAsFactors = FALSE); names(p) <- c("Outcome", "Previous-day symptom", "Participant estimate")
    z <- merge(g, p, by = c("Outcome", "Previous-day symptom")); z$Difference <- z$`Participant estimate` - z$`Group estimate`
    z[order(-abs(z$Difference)), ]
  }, digits = 3)
}

# Starting the app:
app <- shiny::shinyApp(ui, server)
if (sys.nframe() == 0L || interactive()) shiny::runApp(app)

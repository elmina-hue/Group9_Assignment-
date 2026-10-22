# PIPS Group Project: Temporal Symptom Networks

This repository contains our group project for the **Programming in Psychological Science (PIPS)** course in the Research Master's in Psychology at the University of Amsterdam.

The project investigates how psychological symptoms predict one another over time and whether temporal relationships estimated across a group adequately describe the symptom dynamics of an individual.

Using simulated ecological momentary assessment (EMA) data and Bayesian vector autoregressive (VAR) models, we estimate and compare group-level and individual-level temporal symptom networks. The project includes data simulation, statistical analysis, network visualisations, sensitivity analysis, and an interactive Shiny application.

> **Note:** All data are simulated for methodological demonstration. They do not represent observations collected from real participants and should not be interpreted as clinical evidence.

## Research Question

To what extent does a group-level temporal symptom network adequately describe the symptom dynamics of an individual, and which temporal relationships may differ between the individual and the group?

## Theoretical Background

The project follows the network approach to psychopathology, which conceptualises psychological disorders as systems of interconnected symptoms. In a symptom network, symptoms are represented as nodes and relationships between symptoms are represented as edges (Borsboom & Cramer, 2013; Robinaugh et al., 2020).

Ecological momentary assessment (EMA) involves repeatedly measuring people's experiences in everyday life (Myin-Germeys et al., 2018). A first-order vector autoregressive model, or VAR(1), estimates temporal relationships between symptoms measured at consecutive time points.

The temporal network contains two types of relationships:

- **Autoregressive edges:** Represent the persistence of a symptom from one measurement to the next.
- **Cross-lagged edges:** Represent predictive relationships between different symptoms across consecutive measurements.

We distinguish between group-level symptom dynamics, which describe temporal relationships across participants, and individual-level symptom dynamics, which describe relationships estimated for a particular participant. These patterns may differ, illustrating why group-level findings may not fully represent every individual.This distinction has societal and clinical relevance because treatment approaches informed by group-level findings may not adequately reflect the symptom dynamics and needs of every individual. Understanding individual differences may therefore help highlight the limitations of generalising group-level findings to individual treatment needs and the potential value of more personalised approaches to mental health care.

Temporal associations represent predictive relationships and should not automatically be interpreted as causal effects.

## Project Objectives

The project aims to:

- Generate simulated longitudinal data representing nine PHQ-9 symptom domains.
- Estimate group-level temporal relationships using a Bayesian multilevel VAR(1) model.
- Estimate an individual-level symptom network for a selected participant.
- Compare group-level and individual-level temporal coefficients.
- Visualise symptom networks and differences between selected relationships.
- Examine the sensitivity of selected estimates to alternative prior specifications.
- Develop an interactive Shiny application for exploring group and individual symptom networks.
- Discuss the methodological assumptions and limitations of using simulated data.

## Methodology

### 1. Data Simulation

The simulation generates daily symptom ratings for **25 simulated participants across 100 days per participant**.

The nine symptom domains are based on the Patient Health Questionnaire-9 (PHQ-9; Kroenke et al., 2001), a validated measure of depression severity. The original questionnaire assesses nine symptoms over the previous two weeks, with response options scored from 0 to 3.

- `down` — Depressed mood
- `tired` — Fatigue or low energy
- `sleep` — Sleep problems
- `anhedonia` — Loss of interest or pleasure
- `worthlessness` — Feelings of worthlessness or guilt
- `concentration` — Difficulty concentrating
- `appetite` — Appetite changes
- `psychomotor` — Psychomotor slowing or restlessness
- `suicidal` — Thoughts of death or self-harm

The simulation specifies an illustrative group-level temporal coefficient matrix and introduces participant-specific variation around the group-level coefficients. Participant 1 has predefined coefficient differences to demonstrate how an individual temporal network may depart from the group pattern.

Ratings range from **0 to 3**, corresponding to the PHQ-9 item response options. The simulation uses a simplified Gaussian process followed by rounding and clipping to this range. It represents a daily adaptation of PHQ-9 symptom domains rather than the standard PHQ-9 administration, which assesses symptoms over the previous two weeks.

The coefficient values and selected cross-lagged relationships are illustrative simulation assumptions, not empirical estimates of clinical relationships.

The simulation script generates the dataset and calculates rating frequencies for quality checks.

### 2. Statistical Analysis

The statistical analysis uses the `bvarnet` package to estimate Bayesian VAR(1) models.

Two complementary analyses are conducted:

- **Group-level analysis:** A Bayesian multilevel VAR estimates temporal relationships across participants while accounting for between-person variability.
- **Individual-level analysis:** A separate Bayesian VAR estimates temporal relationships using observations from Participant 1.

The analysis extracts temporal network estimates, examines between-person variation in coefficients, and calculates Bayes factors for selected relationships where appropriate. These temporal estimates represent lagged relationships, in which a symptom at the previous time point ((t-1)) predicts a symptom at the current time point ((t)). For example, the analysis examines whether sleep problems on one day predict tiredness on the following day.

A sensitivity analysis examines the lagged relationship between sleep problems at day (t-1) and tiredness at day (t) for Participant 1. The individual-level Bayesian VAR model is refitted using three prior specifications: narrow (0.35), original (0.50), and wide (0.65).

The estimated coefficient and corresponding Bayes factor are compared across the three models. The results are summarised in tables, and a plot displays the posterior means and intervals for the selected coefficient under each prior specification. This analysis assesses whether the estimated relationship is sensitive to the choice of prior.

### 3. Visualisation

Network plots represent the estimated temporal relationships:

- **Nodes** represent symptom domains.
- **Directed edges** represent lagged predictive relationships.
- **Self-loops** represent autoregressive effects.

Group-level and individual-level networks are compared to identify similarities and differences in estimated temporal relationships. Additional plots and tables summarise selected coefficient estimates, uncertainty intervals, and sensitivity to alternative prior specifications.

### 4. Interactive Shiny Application

The Shiny application provides an interactive interface for exploring individual and group temporal symptom networks.

The application is designed to allow users to:

- Generate example participant data or upload a participant dataset in CSV format.
- View group-level and individual-level temporal networks side by side.
- Adjust network display settings, including edge thresholds, self-loops, and edge labels.
- Inspect daily ratings and the number of observations analysed.
- Compare group and participant edge estimates and their differences.

The application uses the project's Bayesian VAR modelling approach to estimate the participant network and compares it with the group reference network.

Uploaded data must contain the required symptom variables and consecutive daily observations, with whole-number ratings from 0 to 3.

## Repository Structure

The repository is organised as follows:

- `data/` — Contains the simulated dataset and rating-frequency table generated by the simulation script.
- `src/simulate_data.R` — Generates the simulated longitudinal symptom data.
- `src/analysis_network.R` — Fits Bayesian VAR models and compares group-level and individual-level temporal relationships.
- `src/plots.R` — Creates symptom network visualisations and comparison plots.
- `main.Rmd` — Main project deliverable, combining the methodological explanation, executable code, analysis, visualisations, results, and discussion.
- `README.md` — Provides the project overview, dependencies, and instructions for running the project.
- `src/symptom_network_app.R` — Shiny application for exploring and comparing group-level and individual-level temporal networks.

The CSV files are generated when the simulation script is run and may not be present in a fresh checkout until that step is completed.

## Requirements and Dependencies

The project requires **R** and the following packages:

- `bvarnet` — Bayesian vector autoregressive modelling.
- `qgraph` — Network visualisation used in the analysis.
- `shiny` — Interactive web application framework.
- `igraph` — Network construction and plotting in the Shiny application.
- `readr` — Reading CSV data.
- `cmdstanr` — R interface for running Bayesian models using CmdStan.

Install the required packages in R:

```r
install.packages(c(
  "bvarnet",
  "qgraph",
  "shiny",
  "igraph",
  "readr",
  "cmdstanr"
))
```

## Getting Started

### 1. Obtain the Repository

Clone or download the repository from GitHub and open the project folder in RStudio.

### 2. Run the Main Deliverable

Open main.Rmd in RStudio and knit the document. The main file automatically checks whether the required R packages are installed and loaded, installing and/or loading them as needed. It also checks and handles the necessary CmdStan setup before running the analysis.

The main file is intended to run the project workflow from beginning to end, including generating the simulated dataset, conducting the statistical analysis, producing visualisations, and presenting the results and discussion.

The simulation script, src/simulate_data.R, generates daily symptom ratings and saves the dataset and rating-frequency table in the data/ directory. It should run before the analysis steps that depend on its output.

### 3. Launch the Shiny Application

From the repository root, run:

```r
shiny::runApp("src/symptom_network_app.R")
```

The application opens an interactive interface for generating or uploading participant data and comparing the participant's temporal network with the group-level network.

## Data and Interpretation

The data are simulated and designed to demonstrate temporal network modelling. They are not real participant data and do not establish clinical relationships.

Although the ratings are represented on the PHQ-9 response scale, the simulation uses a simplified continuous Gaussian process followed by rounding and clipping. The Bayesian analysis also uses a Gaussian model as an approximation to the ordinal ratings. These assumptions should be considered when interpreting the results.

Temporal edges represent estimated predictive relationships between consecutive days, not proof of causation. The group-level network may not fully capture an individual's symptom dynamics, but conclusions from this project remain limited by the illustrative simulation design, model assumptions, and amount of simulated data.

## References

Borsboom, D., & Cramer, A. O. J. (2013). Network analysis: An integrative approach to the structure of psychopathology. *Annual Review of Clinical Psychology, 9*, 91–121. https://doi.org/10.1146/annurev-clinpsy-050212-185608

Kroenke, K., Spitzer, R. L., & Williams, J. B. W. (2001). The PHQ-9: Validity of a brief depression severity measure. *Journal of General Internal Medicine, 16*(9), 606–613. https://doi.org/10.1046/j.1525-1497.2001.016009606.x

Myin-Germeys, I., Kasanova, Z., Vaessen, T., Vachon, H., Kirtley, O., Viechtbauer, W., & Reininghaus, U. (2018). Experience sampling methodology in mental health research: New insights and technical developments. *World Psychiatry, 17*(2), 123–132. https://doi.org/10.1002/wps.20513

Robinaugh, D. J., Hoekstra, R. H. A., Toner, E. R., & Borsboom, D. (2020). The network approach to psychopathology: A review of the literature 2008–2018 and an agenda for future research. *Psychological Medicine, 50*(3), 353–366. https://doi.org/10.1017/S0033291719003404




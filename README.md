# Group9_Assignment-
# PIPS Group Project

## Project Overview

Psychological symptoms, such as depressed mood, fatigue, sleep problems, and nervousness, can change over time and may predict one another across consecutive days. The network approach to psychopathology conceptualises psychological disorders as systems of interconnected symptoms, represented as nodes and edges (Borsboom & Cramer, 2013; Robinaugh et al., 2020).

This project investigates temporal relationships between symptoms using ecological momentary assessment (EMA), which involves repeatedly measuring people's experiences in everyday life (Myin-Germeys et al., 2018). We use a first-order vector autoregressive model (VAR(1)) to examine whether symptoms on one day predict symptom levels on the following day.

The estimated temporal networks contain two types of relationships: autoregressive edges, which represent the persistence of a symptom over time, and cross-lagged edges, which represent predictive relationships between different symptoms across consecutive days.

An important distinction is between group-level and individual-level symptom dynamics. A group-level network describes average temporal relationships across participants, whereas an individual-level network represents the relationships estimated for a particular participant. These patterns may differ, meaning that relationships observed across a group do not necessarily describe every individual's symptom dynamics. This issue is relevant to ergodicity and the generalisation of group-level findings to individuals (Molenaar, 2004; Fisher et al., 2018).

## Research Question

To what extent does the group-level temporal symptom network adequately describe the symptom dynamics of an individual participant, and which relationships may differ between the individual and the group?

## Objectives

- Review scientific literature relevant to the selected psychological topic.
- Select an appropriate questionnaire or psychological measure.
- Generate simulated psychological data using assumptions informed by the literature.
- Conduct an appropriate inferential statistical analysis.
- Create informative visualisations to communicate the results.
- Interpret the findings and discuss the limitations of simulated data.

## Repository Structure

- `README.md`: Project overview, requirements, and instructions.
- `main.Rmd`: Main report combining scientific explanations, code execution, results, visualisations, and discussion.
- `src/simulate_data.R`: Code for generating the simulated dataset.
- `src/analysis.R`: Code for the statistical analysis.
- `src/plots.R`: Code for creating visualisations.
- `data/`: Folder for generated datasets and related documentation.
- `references.md`: Scientific papers and methodological references used in the project.

## Requirements

- R (version to be specified).
- RStudio (recommended).
- Required R packages: To be specified once the analysis and visualisation methods have been selected.

## Getting Started

1. Clone or download this repository.
2. Install R and, if preferred, RStudio.
3. Install the required R packages listed in this README.
4. Open `main.Rmd` in RStudio.
5. Follow the report from top to bottom to run the simulation and statistical analysis and generate the visualisations.

Detailed instructions will be updated once the workflow has been implemented and tested.

## Data

The project will use simulated psychological data. The simulation procedure, variables, and assumptions will be documented in `src/simulate_data.R` and explained in the main report.

Simulated observations will not be presented as data collected from real participants.


## References

Scientific papers, questionnaire sources, and methodological references will be documented in `references.md` and cited in the main report.


# Group9_Assignment-
# PIPS Group Project

1. Project Overview

This repository contains our group project for the Programming in Psychological Science (PIPS) course in the Research Master's in Psychology at the University of Amsterdam. Psychological symptoms, such as depressed mood, fatigue, sleep problems, and nervousness, can change over time and may predict one another across consecutive days. However, the relationships observed across a group of people may not fully represent the symptom dynamics of a particular individual.

This project investigates temporal relationships between psychological symptoms using simulated ecological momentary assessment (EMA) data and Bayesian vector autoregressive (VAR) models. We compare group-level and individual-level symptom networks to explore whether the average pattern of relationships adequately describes an individual participant.

The project also includes an interactive Shiny application that will allow users to generate example participant data or upload a CSV file and explore differences between individual-level and group-level symptom networks.

All simulated data are generated for methodological demonstration and do not represent observations collected from real participants.

2. Research Question

To what extent does the group-level temporal symptom network adequately describe the symptom dynamics of an individual participant, and which relationships may differ between the individual and the group?

3. Theoretical Background

The project follows the network approach to psychopathology, which conceptualises psychological disorders as systems of interconnected symptoms. In a symptom network, symptoms are represented as nodes and relationships between symptoms are represented as edges (Borsboom & Cramer, 2013; Robinaugh et al., 2020).

We use ecological momentary assessment (EMA), which involves repeatedly measuring people's experiences in everyday life, to investigate how symptoms relate to one another over time (Myin-Germeys et al., 2018). A first-order vector autoregressive model (VAR(1)) estimates whether symptom levels on one day predict symptom levels on the following day.

The resulting temporal network contains two types of relationships:

Autoregressive edges: Represent the persistence of a symptom from one measurement to the next.

Cross-lagged edges: Represent predictive relationships between different symptoms across consecutive measurements.

We distinguish between group-level symptom dynamics, which describe average temporal relationships across participants, and individual-level symptom dynamics, which describe relationships estimated for a particular participant. These patterns may differ, highlighting the limitations of assuming that group-level findings apply equally to every individual.

4. Project Objectives

The project aims to:

Generate simulated longitudinal data containing daily ratings of six psychological symptoms.

Estimate group-level temporal relationships using a Bayesian multilevel VAR(1) model.

Estimate an individual-level symptom network and compare it with the group-level network.

Visualise symptom networks and compare selected temporal relationships.

Explore uncertainty in estimated relationships and differences between the group and the individual.

Develop a Shiny app for generating example data or uploading participant data and exploring symptom networks.

Discuss the implications and limitations of using simulated data to demonstrate individual differences.

5. Methodological Approach

5.1 Data Simulation

The initial simulation generates daily ratings for six symptoms across 25 simulated participants, with 100 observations per participant.

The symptom variables are:

down: Depressed mood

tired: Tiredness

sleep: Sleep problems

nervous: Nervousness

worry: Worry

irritable: Irritability

The simulation specifies a group-level pattern of temporal relationships and allows individual participants to differ from this pattern. These assumptions are used to demonstrate how individual-level networks may resemble or differ from the group-level network.

5.2 Statistical Analysis

A Bayesian VAR(1) model is used to estimate temporal relationships between symptoms across consecutive days.

Two analyses are planned:

Group-level analysis: A Bayesian multilevel VAR model estimates average temporal relationships across participants while accounting for between-person variation.

Individual-level analysis: A separate VAR model estimates temporal relationships using observations from an individual participant.

The analyses will compare selected temporal coefficients and their uncertainty to investigate similarities and differences between the individual-level and group-level networks.

Bayesian estimation provides posterior distributions for model parameters and allows comparisons between hypotheses using Bayes factors where appropriate.

5.3 Visualisation


5.4 Interactive Shiny Application

6. Repository Structure

The repository is organised as follows:

README.md: Project overview, methodology, requirements, and instructions.

main.Rmd: Main project deliverable, combining explanations, code execution, analysis, visualisations, results, and discussion.

references.md: Scientific and methodological references.

src/simulate_data.R: Code for generating simulated longitudinal symptom data.

src/analysis.R: Code for fitting the Bayesian VAR models and comparing temporal relationships.

src/plots.R: Code for creating symptom network visualisations and comparison plots.

data/: Directory for the saved example dataset and related documentation.

app.R: Shiny application for generating or uploading participant data and exploring symptom networks.

7. Requirements and Dependencies

  
8. Getting Started
   
8.1 Obtaining the Repository
   
8.2 Installing Dependencies
   
8.3 Running the Main Deliverable
   
8.4 Launching the Shiny App

9. Data
  
   10. Interpretation
  
   11. References
   12. 


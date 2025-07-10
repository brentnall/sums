### PROSTATE - Sensitivity (Stage III-IV) Estimates
setwd("../../../code")
source("Allcancer_functions.R")

# Load data ---------------------------------------------------------------
setwd("../data")
mortality_data <- read_excel("SUMS_Results_Master_final.xlsx", sheet="Master_Mortality_Results")
surrogate_data <- read_excel("SUMS_Results_Master_final.xlsx", sheet="Master_Surrogate_Results")

mortality_data <- tibble::as_tibble(mortality_data)
surrogate_data <- tibble::as_tibble(surrogate_data)


# PRIMARY ANALYSIS - Estimate Selection -----------------------------------

# Mortality  ---------------------------------------------------------------



mortality_data <-
  filter_mortality_data(mortality_data, cancer = "Prostate")

mortality_data_1 <- mortality_data %>%
  #rate ratio > relative risk/HR, for the studies that report PY
  filter(type_indicator == "Rate ratio") %>%
  calculate_logRR() %>%
  replace_with_reported()

mortality_data_2 <- mortality_data %>%
  #relative risk > adjusted HR (Norrkoping); only relative risk available (ERSPC Pilot)
  filter(type_indicator == "Relative risk") %>%
  # already have rate ratio for UK CAP
  filter(trial_acronym != "UK CAP") %>%
  calculate_logrel() %>%
  replace_with_reported()

mortality_data <- rbind(mortality_data_1, mortality_data_2)

# Surrogate 1. ------------------------------------------------------
surrogate_data_1 <-
  filter_surrogate_data(surrogate_data,
                        "Prostate",
                        "1./4. Absolute incidence of late / early stage cancer")

#PLCO
PLCO <- surrogate_data_1 %>%
  filter(trial_acronym == "PLCO (Prostate)") %>%
  #  10y FU
  filter(period_after_rnd == "10") %>%
  filter(stage_category == "III" | stage_category == "IV") %>%
  summarise(
    screening_test = first(screening_test),
    outcome_description = first(outcome_description),
    numerator_screening = sum(as.numeric(numerator_screening)),
    numerator_comparator = sum(as.numeric(numerator_comparator)),
    denominator_screening = mean(as.numeric(denominator_screening)),
    denominator_comparator = mean(as.numeric(denominator_comparator)),
    type_indicator = first(type_indicator),
    reported_definition_late = first(reported_definition_late),
    reported_definition_early = first(reported_definition_early)
  ) %>%
  mutate(trial_acronym = "PLCO (Prostate)", stage_category = "Stage III-IV")

prostate_surr1 <- PLCO %>%
  calculate_logrel()

finaldata_primary_1 <- final_table(mortality_data, prostate_surr1)


# Surrogate 3. -------------------------------------------------------------
surrogate_data_3 <-
  filter_surrogate_data(surrogate_data,
                        "Prostate",
                        "3. % of target cancer diagnosed at late stage")

#PLCO
PLCO <- surrogate_data_3 %>%
  filter(trial_acronym == "PLCO (Prostate)") %>%
  #  10y FU
  filter(period_after_rnd == "10") %>%
  filter(stage_category != "8-10") %>%
  summarise(
    screening_test = first(screening_test),
    numerator_screening = sum(as.numeric(numerator_screening)),
    numerator_comparator = sum(as.numeric(numerator_comparator)),
    denominator_screening = mean(as.numeric(denominator_screening)),
    denominator_comparator = mean(as.numeric(denominator_comparator)),
    comparator_other_screening = first(comparator_other_screening),
    comparator_noscreening = first(comparator_noscreening),
    reported_definition_late = first(reported_definition_late),
    reported_definition_early = first(reported_definition_early)
  ) %>%
  mutate(
    trial_acronym = "PLCO (Prostate)",
    stage_category = "Stage III-IV",
    outcome_description = "% prostate cancers diagnosed at stage III-IV; through 10 years of follow-up"
  )


prostate_surr3 <- PLCO %>%
  calculate_logrel()

finaldata_primary_3 <- final_table(mortality_data, prostate_surr3)

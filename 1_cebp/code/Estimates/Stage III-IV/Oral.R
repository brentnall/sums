### ORAL - Sensitivity (Stage III-IV) Estimates
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
  filter_mortality_data(mortality_data, cancer = "Oral")

mortality_data <- mortality_data %>%
  # rate ratio > relative risk
  filter(type_indicator == "Rate ratio") %>%
  calculate_logRR() %>%
  replace_with_reported()


# Surrogate 1. ------------------------------------------------------
surrogate_data_1 <-
  filter_surrogate_data(surrogate_data,
                        "Oral",
                        "1./4. Absolute incidence of late / early stage cancer")

oral_surr1 <- surrogate_data_1 %>%
  # follow-up time = mortality outcome follow-up
  filter(end_fu_date == "31/12/2004") %>%
  # reported definition of late stage: III + IV
  filter(stage_category == "Stage III-IV") %>%
  #rate ratio > relative risk
  filter(type_indicator == "Rate ratio") %>%
  select(
    trial_acronym,
    screening_test,
    outcome_description,
    numerator_screening,
    numerator_comparator,
    denominator_screening,
    denominator_comparator,
    reported_definition_late,
    reported_definition_early,
    stage_category
  )


oral_surr1 <- calculate_logRR(oral_surr1)

finaldata_primary_1 <- final_table(mortality_data, oral_surr1)


# Surrogate 3. -------------------------------------------------------------
surrogate_data_3 <-
  filter_surrogate_data(surrogate_data,
                        "Oral",
                        "3. % of target cancer diagnosed at late stage")
N <- min_N_trials(surrogate_data_3)
N ### 1 trial

oral_surr3 <- surrogate_data_3 %>%
  filter(end_fu_date == "31/12/2004") %>%
  select(
    trial_acronym,
    screening_test,
    outcome_description,
    numerator_screening,
    numerator_comparator,
    denominator_screening,
    denominator_comparator,
    reported_definition_late,
    reported_definition_early,
    stage_category,
    comparator_other_screening,
    comparator_noscreening
  )

oral_surr3 <- calculate_logrel(oral_surr3)

finaldata_primary_3 <- final_table(mortality_data, oral_surr3)

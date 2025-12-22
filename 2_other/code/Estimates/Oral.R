### ORAL 

setwd("../../code")
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

# Surrogate 4. -------------------------------------------------------------
surrogate_data_4 <-
  filter_surrogate_data(surrogate_data,
                        "Oral",
                        "1./4. Absolute incidence of late / early stage cancer")

oral_surr4 <- surrogate_data_4 %>%
  # follow-up time = mortality outcome follow-up
  filter(end_fu_date == "31/12/2004") %>%
  # reported definition of early stage: I + II
  filter(stage_category == "Stage I-II") %>%
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


oral_surr4 <- calculate_logRR(oral_surr4)

finaldata_primary_4 <- final_table(mortality_data, oral_surr4)


# Surrogate 5. ------------------------------------------------------------
surrogate_data_5 <-
  filter_surrogate_data(surrogate_data,
                        "Oral",
                        "5. % of target cancers that are screen-detected")

oral_surr5 <- surrogate_data_5 %>%
  # follow-up time = mortality outcome follow-up
  filter(end_fu_date == "31/12/2004") %>%
  select(
    trial_acronym,
    screening_test,
    comparator_other_screening,
    outcome_description,
    numerator_screening,
    denominator_screening
  ) %>%
  calculate_prop()

finaldata_primary_5 <- final_table(mortality_data, oral_surr5)


# Surrogate 6. ------------------------------------------------------------
surrogate_data_6 <- filter_surrogate_data(surrogate_data, "Oral","6. % of high-grade target cancers that are screen-detected")
N <- min_N_trials(surrogate_data_6)
N ### 0 trials


# Surrogate 7. ------------------------------------------------------------
surrogate_data_7 <-
  filter_surrogate_data(surrogate_data, "Oral", "7. Diagnostic yield at screening")

oral_surr7 <- surrogate_data_7 %>%
  # follow-up time = mortality outcome follow-up
  filter(end_fu_date == "31/12/2004") %>%
  summarise(
    trial_acronym = first(trial_acronym),
    screening_test = first(screening_test),
    outcome_description = last(outcome_description),
    numerator_screening = sum(as.numeric(numerator_screening)),
    denominator_screening = mean(as.numeric(denominator_screening))
  )

oral_surr7 <- calculate_prop(oral_surr7)

finaldata_primary_7 <- final_table(mortality_data, oral_surr7)

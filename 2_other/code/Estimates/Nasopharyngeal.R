### NASOPHARYNGEAL
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

mortality_data <- filter_mortality_data(mortality_data,cancer="Nasopharyngeal")

mortality_data <- mortality_data %>%
  #rate ratio > relative risk
  filter(type_indicator == "Rate ratio")

mortality_data <- mortality_data %>%
  calculate_logRR() %>%
  replace_with_reported()

# Surrogate 4. -------------------------------------------------------------
surrogate_data_4 <- filter_surrogate_data(surrogate_data, "Nasopharyngeal","1./4. Absolute incidence of late / early stage cancer")

nasoph_surr4 <- surrogate_data_4 %>%
  # late stage definition: stage I+II (reported)
  filter(stage_category == "Stage I or II") %>%
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

#only relative risk reported
nasoph_surr4 <- calculate_logrel(nasoph_surr4)

finaldata_primary_4 <- final_table(mortality_data, nasoph_surr4)

# Surrogate 5. ------------------------------------------------------------
surrogate_data_5 <-
  filter_surrogate_data(surrogate_data,
                        "Nasopharyngeal",
                        "5. % of target cancers that are screen-detected")

nasoph_surr5 <- surrogate_data_5 %>%
  select(
    trial_acronym,
    screening_test,
    comparator_other_screening,
    outcome_description,
    numerator_screening,
    denominator_screening
  ) %>%
  calculate_prop()

finaldata_primary_5 <- final_table(mortality_data, nasoph_surr5)


# Surrogate 6. ------------------------------------------------------------
surrogate_data_6 <- filter_surrogate_data(surrogate_data, "Nasopharyngeal","6. % of high-grade target cancers that are screen-detected")
N <- min_N_trials(surrogate_data_6)
N ### 0 trials

# Surrogate 7. ------------------------------------------------------------
surrogate_data_7 <-
  filter_surrogate_data(surrogate_data,
                        "Nasopharyngeal",
                        "7. Diagnostic yield at screening")

nasoph_surr7 <- surrogate_data_7 %>%
  select(
    trial_acronym,
    screening_test,
    comparator_other_screening,
    outcome_description,
    numerator_screening,
    denominator_screening
  )

nasoph_surr7 <- calculate_prop(nasoph_surr7)

finaldata_primary_7 <- final_table(mortality_data, nasoph_surr7)

### LIVER
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
  filter_mortality_data(mortality_data, cancer = "Liver")

mortality_data <- mortality_data %>%
  #rate ratio > relative risk
  filter(type_indicator == "Rate ratio") %>%
  calculate_logRR() %>%
  replace_with_reported()

# Surrogate 4. -------------------------------------------------------------
surrogate_data_4 <- filter_surrogate_data(surrogate_data, "Liver","1./4. Absolute incidence of late / early stage cancer")

Qidong <- surrogate_data_4 %>%
  filter(trial_acronym == "Qidong Liver") %>%
  # early stage definition: stage I (reported)
  filter(stage_category == "Stage I") %>%
  # rate ratio > relative risk
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
Qidong <- calculate_logRR(Qidong)

#Shanghai removed because (0 cancers detected in control arm)
Shanghai <- surrogate_data_4 %>%
  filter(trial_acronym == "Shanghai Liver") %>%
  # early stage definition: stage I (reported)
  filter(stage_category == "Stage I") %>%
  # 0 cancers detcted in control arm - give value 1
  mutate(numerator_comparator=1) %>%
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
Shanghai <- calculate_logrel(Shanghai)

liver_surr4 <- rbind(Qidong, Shanghai)


finaldata_primary_4 <- final_table(mortality_data, liver_surr4)

# Surrogate 5. ------------------------------------------------------------
surrogate_data_5 <-
  filter_surrogate_data(surrogate_data,
                        "Liver",
                        "5. % of target cancers that are screen-detected")

liver_surr5 <- surrogate_data_5 %>%
  #exclude Qidong because screened control arm
  filter(trial_acronym!="Qidong Liver") %>%
  select(
    trial_acronym,
    screening_test,
    comparator_other_screening,
    outcome_description,
    numerator_screening,
    denominator_screening
  ) %>%
  calculate_prop()

finaldata_primary_5 <- final_table(mortality_data, liver_surr5)

# Surrogate 6. ------------------------------------------------------------
surrogate_data_6 <- filter_surrogate_data(surrogate_data, "Liver","6. % of high-grade target cancers that are screen-detected")
N <- min_N_trials(surrogate_data_6)
N ### 0 trials

# Surrogate 7. ------------------------------------------------------------
surrogate_data_7 <-
  filter_surrogate_data(surrogate_data, "Liver", "7. Diagnostic yield at screening")

liver_surr7 <- surrogate_data_7 %>%
  #exclude Qidong because screened control arm
  filter(trial_acronym!="Qidong Liver") %>%
  select(
    trial_acronym,
    screening_test,
    comparator_other_screening,
    outcome_description,
    numerator_screening,
    denominator_screening
  ) %>%
  calculate_prop()

finaldata_primary_7 <- final_table(mortality_data, liver_surr7)

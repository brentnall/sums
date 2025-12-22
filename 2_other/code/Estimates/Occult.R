### OCCULT 

setwd("../../code")
source("Allcancer_functions.R")

# Load data ---------------------------------------------------------------
setwd("../data")
mortality_data <- read_excel("SUMS_Results_Master_final.xlsx", sheet="Master_Mortality_Results")
surrogate_data <- read_excel("SUMS_Results_Master_final.xlsx", sheet="Master_Surrogate_Results")

mortality_data <- tibble::as_tibble(mortality_data)
surrogate_data <- tibble::as_tibble(surrogate_data)


# PRIMARY ANALYSIS - Estimate Selection -----------------------------------

# Mortality---------------------------------------------------------------
mortality_data <-
  filter_mortality_data(mortality_data, cancer = "Occult cancer")

# only relative risks available
mortality_data <- mortality_data %>%
  calculate_logrel() %>%
  replace_with_reported()


# Surrogate 4. -------------------------------------------------------------
surrogate_data_4 <- filter_surrogate_data(surrogate_data, "Occult cancer","1./4. Absolute incidence of late / early stage cancer")


#SOMIT
SOMIT <- surrogate_data_4 %>%
  filter(trial_acronym == "SOMIT") %>%
  #early stage definition: T1 or T2; N0 and M0 (reported)
  filter(stage_category == "T1 or T2; N0 and M0") %>%
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


#D'Aquapendente
DAquapendente <- surrogate_data_4 %>%
  filter(trial_acronym == "D'Aquapendente") %>%
  #early stage definition: stage I-III (advanced = stage IV reported)
  filter(stage_category != "Stage IV") %>%
  summarise(
    screening_test = first(screening_test),
    outcome_description = first(outcome_description),
    numerator_screening = sum(as.numeric(numerator_screening)),
    numerator_comparator = sum(as.numeric(numerator_comparator)),
    denominator_screening = mean(as.numeric(denominator_screening)),
    denominator_comparator = mean(as.numeric(denominator_comparator)),
    reported_definition_late = first(reported_definition_late),
    reported_definition_early = first(reported_definition_early)
  ) %>%
  mutate(trial_acronym = "D'Aquapendente", stage_category = "Stage I-III")

# MVTEP
MVTEP <- surrogate_data_4 %>%
  filter(trial_acronym == "MVTEP") %>%
  #early stage definition: T1-2N0M0 (reported)
  filter(stage_category == "Early (T1-2N0M0)") %>%
  filter(str_detect(outcome_description, "Absolute incidence")) %>%
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

# SOME
SOME <- surrogate_data_4 %>%
  filter(trial_acronym == "SOME") %>%
  #early stage definition: T1-2, N0, M0 (reported)
  filter(stage_category == "T1-2, N0, M0") %>%
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

occult_surr4 <- rbind(SOMIT, DAquapendente, MVTEP, SOME)

occult_surr4 <- calculate_logrel(occult_surr4)

finaldata_primary_4 <- final_table(mortality_data, occult_surr4)


# Surrogate 5. ------------------------------------------------------------
surrogate_data_5 <-
  filter_surrogate_data(surrogate_data,
                        "Occult cancer",
                        "5. % of target cancers that are screen-detected")

occult_surr5 <- surrogate_data_5 %>%
  select(
    trial_acronym,
    screening_test,
    comparator_other_screening,
    outcome_description,
    numerator_screening,
    denominator_screening
  ) %>%
  calculate_prop() %>%
  mutate(comparator_other_screening="Not applicable")

finaldata_primary_5 <- final_table(mortality_data, occult_surr5)


# Surrogate 6. ------------------------------------------------------------
surrogate_data_6 <- filter_surrogate_data(surrogate_data, "Occult cancer","6. % of high-grade target cancers that are screen-detected")
N <- min_N_trials(surrogate_data_6)
N ### 0 trials


# Surrogate 7. ------------------------------------------------------------
surrogate_data_7 <-
  filter_surrogate_data(surrogate_data,
                        "Occult cancer",
                        "7. Diagnostic yield at screening")

occult_surr7 <- surrogate_data_7 %>%
  select(
    trial_acronym,
    screening_test,
    comparator_other_screening,
    outcome_description,
    numerator_screening,
    denominator_screening
  ) %>%
  mutate(comparator_other_screening= "Not applicable")

occult_surr7 <- calculate_prop(occult_surr7)

finaldata_primary_7 <- final_table(mortality_data, occult_surr7)

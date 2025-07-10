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


# Surrogate 1. ------------------------------------------------------
surrogate_data_1 <-
  filter_surrogate_data(surrogate_data,
                        "Occult cancer",
                        "1./4. Absolute incidence of late / early stage cancer")
surrogate_data_1 <- surrogate_data_1 %>%
  #SOME only reports early stage
  filter(trial_acronym != "SOME")

#SOMIT
SOMIT <- surrogate_data_1 %>%
  filter(trial_acronym == "SOMIT") %>%
  #late stage definition: N1 and/or M1
  filter(str_detect(stage_category, "N1")) %>%
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
DAquapendente <- surrogate_data_1 %>%
  filter(trial_acronym == "D'Aquapendente") %>%
  #late stage definition: Stage IV (=advanced reported)
  filter(stage_category == "Stage IV") %>%
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

# MVTEP
MVTEP <- surrogate_data_1 %>%
  filter(trial_acronym == "MVTEP") %>%
  #late stage definition: Advanced
  filter(stage_category == "Advanced") %>%
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

occult_surr1 <- rbind(SOMIT, DAquapendente, MVTEP)

occult_surr1 <- calculate_logrel(occult_surr1)

finaldata_primary_1 <- final_table(mortality_data, occult_surr1)



# Surrogate 3. -------------------------------------------------------------
surrogate_data_3 <-
  filter_surrogate_data(surrogate_data,
                        "Occult cancer",
                        "3. % of target cancer diagnosed at late stage")

occult_surr3 <- surrogate_data_3 %>%
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

# SOMIT
# calculate surrogate 3 using numbers from surrogate 1 file
SOMIT <- surrogate_data_1 %>%
  filter(trial_acronym == "SOMIT") %>%
  filter(stage_category != "T3") %>%
  #calculate all cancers diagnosed by arm and use as denominator (all stages + missing)
  mutate(denominator_screening = sum(as.numeric(numerator_screening))) %>%
  mutate(denominator_comparator = sum(as.numeric(numerator_comparator))) %>%
  filter(stage_category == "N1 and/or M1") %>%
  mutate(outcome_description = "% cancers with loco-regional or distant metastasis (N1 and/or M1); follow-up for 24 months; Extensive screening vs usual care") %>%
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

occult_surr3 <- rbind(occult_surr3, SOMIT)

occult_surr3 <- calculate_logrel(occult_surr3)

finaldata_primary_3 <- final_table(mortality_data, occult_surr3)

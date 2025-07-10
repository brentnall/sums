### OCCULT - Sensitivity (Stage III-IV) Estimates
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


#D'Aquapendente
DAquapendente <- surrogate_data_1 %>%
  filter(trial_acronym == "D'Aquapendente") %>%
  filter(stage_category == "Stage III" |
           stage_category == "Stage IV") %>%
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
  mutate(trial_acronym = "D'Aquapendente", stage_category = "Stage III-IV")

occult_surr1 <-  DAquapendente %>%
  calculate_logrel()

finaldata_primary_1 <- final_table(mortality_data, occult_surr1)


# Surrogate 3. -------------------------------------------------------------
surrogate_data_3 <-
  filter_surrogate_data(surrogate_data,
                        "Occult cancer",
                        "3. % of target cancer diagnosed at late stage")

#D'Aquapendente
# calculate surrogate 3 using numbers from surrogate 1 file
DAquapendente <- surrogate_data_1 %>%
  filter(trial_acronym == "D'Aquapendente") %>%
  #calculate all cancers diagnosed by arm and use as denominator (all stages + missing)
  mutate(denominator_screening = sum(as.numeric(numerator_screening)) +
           as.numeric(missing_screening)) %>%
  mutate(denominator_comparator = sum(as.numeric(numerator_comparator)) +
           as.numeric(missing_comparator))


DAquapendente_denominator_s <-
  DAquapendente$denominator_screening[1]
DAquapendente_denominator_c <-
  DAquapendente$denominator_comparator[1]

DAquapendente <- surrogate_data_1 %>%
  filter(trial_acronym == "D'Aquapendente") %>%
  filter(stage_category == "Stage III" |
           stage_category == "Stage IV") %>%
  summarise(
    screening_test = first(screening_test),
    numerator_screening = sum(as.numeric(numerator_screening)),
    numerator_comparator = sum(as.numeric(numerator_comparator)),
    denominator_screening = DAquapendente_denominator_s,
    denominator_comparator = DAquapendente_denominator_c,
    type_indicator = first(type_indicator),
    reported_definition_late = first(reported_definition_late),
    reported_definition_early = first(reported_definition_early)
  ) %>%
  mutate(
    trial_acronym = "DAquapendente",
    stage_category = "Stage III-IV",
    outcome_description = "% stage III-IV cancers; follow-up for 24 months; Torso CT + FOBT vs usual care (involving non-standardised screening)"
  )


occult_surr3 <- DAquapendente %>%
  calculate_logrel()

finaldata_primary_3 <- final_table(mortality_data, occult_surr3)

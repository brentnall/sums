### OVARIAN - Sensitivity (Stage III-IV) Estimates
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
  filter_mortality_data(mortality_data, cancer = "Ovarian")

mortality_data <- mortality_data %>%
  calculate_logrel() %>%
  replace_with_reported()

# Surrogate 1. ------------------------------------------------------

surrogate_data_1 <-
  filter_surrogate_data(surrogate_data,
                        "Ovarian",
                        "1./4. Absolute incidence of late / early stage cancer")


#UK Pilot Ovarian
UK_Pilot_Ovarian <- surrogate_data_1 %>%
  filter(trial_acronym == "UK Pilot Ovarian") %>%
  # reported stage category (note, not sure if those numbers were reported or the sum of ind stages was calculated during extraction)
  filter(stage_category == "Stage III-IV") %>%
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
  ) %>%
  calculate_logrel()


#UKCTOCS
UKCTOCS <- surrogate_data_1 %>%
  filter(trial_acronym == "UKCTOCS") %>%
  # use same cancer definition as mortality outcome (=ovarian cancer)
  filter(!str_detect(outcome_description, "ovarian and")) %>%
  # follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[mortality_data$trial_acronym ==
                                                     "UKCTOCS"]) %>%
  # late stage defined as stage III+IV (only individual stage categories available for this cutoff)
  filter(stage_category == "Stage III" |
           stage_category == "Stage IV") %>%
  filter(str_detect(denominator_definition, "Person-years")) %>%
  group_by(screening_test) %>%
  summarise(
    outcome_description = first(outcome_description),
    numerator_screening = sum(as.numeric(numerator_screening)),
    numerator_comparator = sum(as.numeric(numerator_comparator)),
    denominator_screening = mean(as.numeric(denominator_screening)),
    denominator_comparator = mean(as.numeric(denominator_comparator)),
    reported_definition_late = first(reported_definition_late),
    reported_definition_early = first(reported_definition_early)
  ) %>%
  mutate(trial_acronym = "UKCTOCS", stage_category = "Stage III-IV") %>%
  calculate_logRR()


#PLCO
PLCO <- surrogate_data_1 %>%
  filter(trial_acronym == "PLCO (Ovarian)") %>%
  # follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[mortality_data$trial_acronym ==
                                                     "PLCO (Ovarian)"]) %>%
  filter(period_after_rnd == mortality_data$period_after_rnd[mortality_data$trial_acronym ==
                                                               "PLCO (Ovarian)"]) %>%
  # late stage defined as stage III+IV (only individual stage categories available)
  filter(stage_category == "Stage III" |
           stage_category == "Stage IV") %>%
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
  mutate(trial_acronym = "PLCO (Ovarian)", stage_category = "Stage III-IV") %>%
  calculate_logRR()

ovarian_surr1 <- rbind(UK_Pilot_Ovarian, UKCTOCS, PLCO)

finaldata_primary_1 <- final_table(mortality_data, ovarian_surr1)



# Surrogate 3. -------------------------------------------------------------
surrogate_data_3 <-
  filter_surrogate_data(surrogate_data,
                        "Ovarian",
                        "3. % of target cancer diagnosed at late stage")


#PLCO
PLCO <- surrogate_data_3 %>%
  filter(trial_acronym == "PLCO (Ovarian)") %>%
  # follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[mortality_data$trial_acronym ==
                                                     "PLCO (Ovarian)"]) %>%
  filter(period_after_rnd == mortality_data$period_after_rnd[mortality_data$trial_acronym ==
                                                               "PLCO (Ovarian)"]) %>%
  # late stage defined as stage III+IV (only individual stage categories available)
  filter(stage_category == "Stage III" |
           stage_category == "Stage IV") %>%
  summarise(
    screening_test = first(screening_test),
    outcome_description = first(outcome_description),
    numerator_screening = sum(as.numeric(numerator_screening)),
    numerator_comparator = sum(as.numeric(numerator_comparator)),
    denominator_screening = mean(as.numeric(denominator_screening)),
    denominator_comparator = mean(as.numeric(denominator_comparator)),
    reported_definition_late = first(reported_definition_late),
    reported_definition_early = first(reported_definition_early),
    comparator_other_screening = first(comparator_other_screening),
    comparator_noscreening = first(comparator_noscreening)
  ) %>%
  mutate(trial_acronym = "PLCO (Ovarian)", stage_category = "Stage III-IV")

#UK Pilot Ovarian
# calculate surrogate 3 using numbers from surrogate 1 file
UK_Pilot_Ovarian <- surrogate_data_1 %>%
  filter(trial_acronym == "UK Pilot Ovarian") %>%
  filter(str_detect(stage_category, "-")) %>%
  #calculate all cancers diagnosed by arm and use as denominator (all stages + missing)
  mutate(denominator_screening = sum(as.numeric(numerator_screening)) +
           as.numeric(missing_screening)) %>%
  mutate(denominator_comparator = sum(as.numeric(numerator_comparator)) +
           as.numeric(missing_comparator)) %>%
  filter(stage_category == "Stage III-IV") %>%
  mutate(outcome_description = "% stage III-IV ovarian cancers; Follow-up until 1997; Screening vs no screening") %>%
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

#UKCTOCS
# calculate surrogate 3 using numbers from surrogate 1 file
UKCTOCS_MMS <- surrogate_data_1 %>%
  filter(trial_acronym == "UKCTOCS") %>%
  filter(str_detect(outcome_description, "MMS")) %>%
  # use same cancer definition as mortality outcome (=ovarian cancer)
  filter(!str_detect(outcome_description, "ovarian and")) %>%
  # follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[mortality_data$trial_acronym ==
                                                     "UKCTOCS"]) %>%
  filter(type_indicator == "Rate ratio") %>%
  #calculate all cancers diagnosed by arm and use as denominator (all stages + missing)
  mutate(denominator_screening = sum(as.numeric(numerator_screening)) +
           as.numeric(missing_screening)) %>%
  mutate(denominator_comparator = sum(as.numeric(numerator_comparator)) +
           as.numeric(missing_comparator)) %>%
  filter(stage_category == "Stage III" |
           stage_category == "Stage IV") %>%
  summarise(
    screening_test = first(screening_test),
    numerator_screening = sum(as.numeric(numerator_screening)),
    numerator_comparator = sum(as.numeric(numerator_comparator)),
    denominator_screening = mean(as.numeric(denominator_screening)),
    denominator_comparator = mean(as.numeric(denominator_comparator)),
    reported_definition_late = first(reported_definition_late),
    reported_definition_early = first(reported_definition_early),
    comparator_other_screening = first(comparator_other_screening),
    comparator_noscreening = first(comparator_noscreening)
  ) %>%
  mutate(
    trial_acronym = "UKCTOCS",
    stage_category = "Stage III-IV",
    outcome_description = "% stage III-IV primary ovarian cancers; follow-up until 31/12/2014; MMS vs no screening"
  )

UKCTOCS_USS <- surrogate_data_1 %>%
  filter(trial_acronym == "UKCTOCS") %>%
  filter(str_detect(outcome_description, "USS")) %>%
  # use same cancer definition as mortality outcome (=ovarian cancer)
  filter(!str_detect(outcome_description, "ovarian and")) %>%
  # follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[mortality_data$trial_acronym ==
                                                     "UKCTOCS"]) %>%
  filter(type_indicator == "Rate ratio") %>%
  #calculate all cancers diagnosed by arm and use as denominator (all stages + missing)
  mutate(denominator_screening = sum(as.numeric(numerator_screening)) +
           as.numeric(missing_screening)) %>%
  mutate(denominator_comparator = sum(as.numeric(numerator_comparator)) +
           as.numeric(missing_comparator)) %>%
  filter(stage_category == "Stage III" |
           stage_category == "Stage IV") %>%
  summarise(
    screening_test = first(screening_test),
    numerator_screening = sum(as.numeric(numerator_screening)),
    numerator_comparator = sum(as.numeric(numerator_comparator)),
    denominator_screening = mean(as.numeric(denominator_screening)),
    denominator_comparator = mean(as.numeric(denominator_comparator)),
    reported_definition_late = first(reported_definition_late),
    reported_definition_early = first(reported_definition_early),
    comparator_other_screening = first(comparator_other_screening),
    comparator_noscreening = first(comparator_noscreening)
  ) %>%
  mutate(
    trial_acronym = "UKCTOCS",
    stage_category = "Stage III-IV",
    outcome_description = "% stage III-IV primary ovarian cancers; follow-up until 31/12/2014; USS vs no screening"
  )


ovarian_surr3 <-
  rbind(PLCO, UK_Pilot_Ovarian, UKCTOCS_MMS, UKCTOCS_USS)

ovarian_surr3 <- calculate_logrel(ovarian_surr3)

finaldata_primary_3 <- final_table(mortality_data, ovarian_surr3)


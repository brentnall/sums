setwd("../../../code")
source("Allcancer_functions.R")

# Load data ---------------------------------------------------------------
setwd("../data")
mortality_data <- read_excel("SUMS_Results_Master_final.xlsx", sheet="Master_Mortality_Results")
surrogate_data <- read_excel("SUMS_Results_Master_final.xlsx", sheet="Master_Surrogate_Results")

mortality_data <- tibble::as_tibble(mortality_data)
surrogate_data <- tibble::as_tibble(surrogate_data)

setwd("../code/Estimates")
source("Lung.R")

# PRIMARY ANALYSIS ESTIMATES ----------------------------------------------

#select trials that have a midpoint estimate (earlier than mortality FU)
mid_timepoint <- finaldata_primary_1 %>%
  filter(
    trial_acro_sm == "PLCO (Lung)" |
      trial_acro_sm == "UKLS" |
      trial_acro_sm == "Lung Screening Study (LSS)" |
      trial_acro_sm == "ITALUNG" |
      trial_acro_sm == "Czech study" |
      trial_acro_sm == "LUSI"
  ) %>%
  mutate(timepoint = "Earlier than mortality")

# SAME FU AS MORTALITY ----------------------------------------------------

#for trials we used surrogate measured at time of mortality
main_timepoint_1 <- finaldata_primary_1 %>%
  filter(
    trial_acro_sm != "PLCO (Lung)" &
      trial_acro_sm != "UKLS" &
      trial_acro_sm != "Lung Screening Study (LSS)" &
      trial_acro_sm != "ITALUNG" &
      trial_acro_sm != "Czech study" &
      trial_acro_sm != "LUSI"
  )


PLCO <-
  filter_surrogate_data(surrogate_data,
                        "Lung",
                        "1./4. Absolute incidence of late / early stage cancer") %>%
  filter(trial_acronym == "PLCO (Lung)") %>%
  #late stage definition: stage III+IV non-small lung cancer
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
  mutate(trial_acronym = "PLCO (Lung)", stage_category = "Stage III-IV (non-small)")

UKLS <-
  filter_surrogate_data(surrogate_data,
                        "Lung",
                        "1./4. Absolute incidence of late / early stage cancer") %>%
  filter(trial_acronym == "UKLS") %>%
  filter(str_detect(outcome_description, "(stage III-IV)")) %>%
  filter(str_detect(outcome_description, "31/12/201")) %>%
  filter(type_indicator == "Rate ratio") %>%
  select(
    trial_acronym,
    screening_test,
    outcome_description,
    numerator_screening,
    numerator_comparator,
    denominator_screening,
    denominator_comparator,
    type_indicator,
    reported_definition_late,
    reported_definition_early,
    stage_category
  )


# LSS # no estimate available

# ITALUNG # no estimate available

# czech_study # no estimate available

LUSI <-
  filter_surrogate_data(surrogate_data,
                        "Lung",
                        "1./4. Absolute incidence of late / early stage cancer") %>%
  filter(trial_acronym == "LUSI") %>%
  filter(subgroup == "Whole trial") %>%
  filter(str_detect(outcome_description, "30 April 2018")) %>%
  filter(type_indicator == "Relative risk") %>%
  select(
    trial_acronym,
    screening_test,
    outcome_description,
    numerator_screening,
    numerator_comparator,
    denominator_screening,
    denominator_comparator,
    type_indicator,
    reported_definition_late,
    reported_definition_early,
    stage_category
  )


lung_surr1 <- rbind(PLCO, UKLS, LUSI)

lung_surr1_rate <- lung_surr1 %>%
  filter(type_indicator == "Rate ratio") %>%
  calculate_logRR()

lung_surr1_risk <- lung_surr1 %>%
  filter(type_indicator == "Relative risk") %>%
  calculate_logrel()

lung_surr1 <- rbind(lung_surr1_rate, lung_surr1_risk)

lung_surr1 <- final_table(mortality_data, lung_surr1)

main_timepoint <- rbind(main_timepoint_1, lung_surr1) %>%
  mutate(timepoint = "Mortality")


# ONGOING SCREENING --------------------------------------------------------
DLCST <-
  filter_surrogate_data(surrogate_data,
                        "Lung",
                        "1./4. Absolute incidence of late / early stage cancer") %>%
  filter(trial_acronym == "DLCST") %>%
  #FU until 2010
  filter(end_fu_date == "40268") %>%
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
  mutate(trial_acronym = "DLCST", stage_category = "Stage III-IV (non-small)")

# Johns_Hopkins  #no estimate available

# PLCO # no estimate available

UKLS <-
  filter_surrogate_data(surrogate_data,
                        "Lung",
                        "1./4. Absolute incidence of late / early stage cancer") %>%
  filter(trial_acronym == "UKLS") %>%
  filter(str_detect(outcome_description, "1 year since study entry")) %>%
  select(
    trial_acronym,
    screening_test,
    outcome_description,
    numerator_screening,
    numerator_comparator,
    denominator_screening,
    denominator_comparator,
    type_indicator,
    reported_definition_late,
    reported_definition_early,
    stage_category
  )

# NELSON #no estimate available

# LSS # no estimate available

MSK <-
  filter_surrogate_data(surrogate_data,
                        "Lung",
                        "1./4. Absolute incidence of late / early stage cancer") %>%
  filter(trial_acronym == "Memorial Sloan-Kettering") %>%
  filter(str_detect(outcome_description, "whole screening period")) %>%
  filter(str_detect(outcome_description, "stage II-III")) %>%
  select(
    trial_acronym,
    screening_test,
    outcome_description,
    numerator_screening,
    numerator_comparator,
    denominator_screening,
    denominator_comparator,
    type_indicator,
    reported_definition_late,
    reported_definition_early,
    stage_category
  )

DANTE <-
  filter_surrogate_data(surrogate_data,
                        "Lung",
                        "1./4. Absolute incidence of late / early stage cancer") %>%
  filter(trial_acronym == "DANTE") %>%
  #screening ongoing
  filter(str_detect(outcome_description, "2008")) %>%
  filter(
    stage_category == "Stage II" |
      stage_category == "Stage IIIA" |
      stage_category == "Stage IIIB"  | stage_category == "Stage IV"
  ) %>%
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
  mutate(trial_acronym = "DANTE", stage_category = "Stage II-IV")

NLST <-
  filter_surrogate_data(surrogate_data,
                        "Lung",
                        "1./4. Absolute incidence of late / early stage cancer") %>%
  filter(trial_acronym == "NLST") %>%
  filter(subgroup == "Whole trial") %>%
  filter(
    stage_category == "Stage IIIA" |
      stage_category == "Stage IIIB"  | stage_category == "Stage IV"
  ) %>%
  filter(str_detect(outcome_description, "screening round")) %>%
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
  mutate(trial_acronym = "NLST", stage_category = "Stage III-IV")

# MILD #no estimate

# ITALUNG # no estimate available

czech_study <-
  filter_surrogate_data(surrogate_data,
                        "Lung",
                        "1./4. Absolute incidence of late / early stage cancer") %>%
  filter(trial_acronym == "Czech study") %>%
  filter(period_after_rnd == "3 years") %>%
  filter(stage_category == "Stage III") %>%
  select(
    trial_acronym,
    screening_test,
    outcome_description,
    numerator_screening,
    numerator_comparator,
    denominator_screening,
    denominator_comparator,
    type_indicator,
    reported_definition_late,
    reported_definition_early,
    stage_category
  )

LUSI <-
  filter_surrogate_data(surrogate_data,
                        "Lung",
                        "1./4. Absolute incidence of late / early stage cancer") %>%
  filter(trial_acronym == "LUSI") %>%
  filter(subgroup == "Whole trial") %>%
  filter(str_detect(outcome_description, "5 years after randomisation")) %>%
  filter(str_detect(outcome_description, "stage II-IV")) %>%
  filter(type_indicator == "Relative risk") %>%
  select(
    trial_acronym,
    screening_test,
    outcome_description,
    numerator_screening,
    numerator_comparator,
    denominator_screening,
    denominator_comparator,
    type_indicator,
    reported_definition_late,
    reported_definition_early,
    stage_category
  )

Mayo <-
  filter_surrogate_data(surrogate_data,
                        "Lung",
                        "1./4. Absolute incidence of late / early stage cancer") %>%
  filter(trial_acronym == "Mayo Lung Project") %>%
  filter(str_detect(outcome_description, "01/04/1980")) %>%
  filter(stage_category == "Stage III") %>%
  filter(type_indicator == "Relative risk") %>%
  select(
    trial_acronym,
    screening_test,
    outcome_description,
    numerator_screening,
    numerator_comparator,
    denominator_screening,
    denominator_comparator,
    type_indicator,
    reported_definition_late,
    reported_definition_early,
    stage_category
  )


lung_surr1 <-
  rbind(DLCST, UKLS, MSK, DANTE, NLST, czech_study, LUSI, Mayo) %>%
  calculate_logrel()

screening_timepoint <- final_table(mortality_data, lung_surr1) %>%
  mutate(timepoint = "Intervention phase")

# COMBINE ALL ESTIMATES ---------------------------------------------------
finaldata_subgroup_1 <-
  rbind(mid_timepoint, main_timepoint, screening_timepoint)

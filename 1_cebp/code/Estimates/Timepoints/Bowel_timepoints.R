setwd("../../../code")
source("Allcancer_functions.R")

# Load data ---------------------------------------------------------------
setwd("../data")
mortality_data <- read_excel("SUMS_Results_Master_final.xlsx", sheet="Master_Mortality_Results")
surrogate_data <- read_excel("SUMS_Results_Master_final.xlsx", sheet="Master_Surrogate_Results") %>%
  #rename for consistency
  mutate(trial_acronym=ifelse(trial_acronym=="PLCO (Bowel)", "PLCO (CRC)", trial_acronym))

mortality_data <- tibble::as_tibble(mortality_data)
surrogate_data <- tibble::as_tibble(surrogate_data)

setwd("../code/Estimates")
source("Bowel.R")

# PRIMARY ANALYSIS ESTIMATES ----------------------------------------------
#select trials that have a midpoint estimate (earlier than mortality FU)
mid_timepoint <- finaldata_primary_1 %>%
  filter(
    trial_acro_sm == "SCORE" |
      trial_acro_sm == "Finnish RHS CRC" |
      trial_acro_sm == "NORCCAP" |
      trial_acro_sm == "Telemark Polyp I"
  ) %>%
  mutate(timepoint = "Earlier than mortality")


# SAME FU AS MORTALITY ----------------------------------------------------
main_timepoint_1 <- crc_surr1 %>%
  #for all other trials we used surrogate measured at time of mortality
  filter(
    trial_acronym != "SCORE" &
      trial_acronym != "Finnish RHS CRC" &
      trial_acronym != "NORCCAP" &
      trial_acronym != "Telemark Polyp I"
  )

#from those with mid-point available, only Telemark also has mortality timepoint
main_timepoint_2 <-
  filter_surrogate_data(surrogate_data,
                        "Bowel",
                        "1./4. Absolute incidence of late / early stage cancer") %>%
  filter(trial_acronym == "Telemark Polyp I") %>%
  filter(end_fu_date == "1995") %>%
  filter(stage_category == "C or D") %>%
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
  ) %>%
  calculate_logrel()

crc_surr1 <- rbind(main_timepoint_1, main_timepoint_2)

main_timepoint <- final_table(mortality_data, crc_surr1) %>%
  mutate(timepoint = "Mortality")



# ONGOING SCREENING --------------------------------------------------------

# SCORE  # no earlier estimate available

FinnishRHS <-
  filter_surrogate_data(surrogate_data,
                        "Bowel",
                        "1./4. Absolute incidence of late / early stage cancer") %>%
  filter(trial_acronym == "Finnish RHS CRC") %>%
  #whole trial estimate
  filter(subgroup == "Whole trial") %>%
  filter(str_detect(outcome_description, "follow-up until 31/12/2007")) %>%
  filter(str_detect(stage_category, "Non-localised")) %>%
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

Funen <-
  filter_surrogate_data(surrogate_data,
                        "Bowel",
                        "1./4. Absolute incidence of late / early stage cancer") %>%
  filter(trial_acronym == "Funen") %>%
  filter(period_after_rnd == "Up to 38 months") %>%
  filter(stage_category == "C" |
           stage_category == "Distant spread (now D)") %>%
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
  mutate(trial_acronym = "Funen", stage_category = "C/Distant spread")

Gothenburg <-
  filter_surrogate_data(surrogate_data,
                        "Bowel",
                        "1./4. Absolute incidence of late / early stage cancer") %>%
  filter(trial_acronym == "Gothenburg CRC") %>%
  filter(stage_category == "D") %>%
  filter(str_detect(period_after_rnd, "2-7 years")) %>%
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

# Minnesota  #no earlier estimate available

Nottingham <-
  filter_surrogate_data(surrogate_data,
                        "Bowel",
                        "1./4. Absolute incidence of late / early stage cancer") %>%
  filter(trial_acronym == "Nottingham") %>%
  filter(period_after_rnd == "Up to 5 years") %>%
  filter(type_indicator == "Relative risk") %>%
  filter(stage_category == "C" | stage_category == "D") %>%
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
  mutate(trial_acronym = "Nottingham", stage_category = "C/D")

# PLCO #no earlier estimate available

# UKFSST #no earlier estimate available

# NORCCAP # no earlier estimate available

# NordICC #no earlier estimate available

# Telemark #no earlier estimate available

crc_surr1 <- rbind(FinnishRHS, Funen, Gothenburg, Nottingham) %>%
  calculate_logrel()

screening_timepoint <- final_table(mortality_data, crc_surr1) %>%
  mutate(timepoint = "Intervention phase")

# COMBINE ALL ESTIMATES ---------------------------------------------------
finaldata_subgroup_1 <-
  rbind(mid_timepoint, main_timepoint, screening_timepoint)


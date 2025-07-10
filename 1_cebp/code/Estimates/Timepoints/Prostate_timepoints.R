setwd("../../../code")
source("Allcancer_functions.R")

# Load data ---------------------------------------------------------------
setwd("../data")
mortality_data <- read_excel("SUMS_Results_Master_final.xlsx", sheet="Master_Mortality_Results")
surrogate_data <- read_excel("SUMS_Results_Master_final.xlsx", sheet="Master_Surrogate_Results")

mortality_data <- tibble::as_tibble(mortality_data)
surrogate_data <- tibble::as_tibble(surrogate_data)

setwd("../code/Estimates")
source("Prostate.R")

# PRIMARY ANALYSIS ESTIMATES ----------------------------------------------

#select trials that have a midpoint estimate (earlier than mortality FU)
mid_timepoint <- finaldata_primary_1 %>%
  filter(
    trial_acro_sm == "UK CAP" |
      trial_acro_sm == "PLCO (Prostate)" |
      trial_acro_sm == "Norrkoping"
  ) %>%
  mutate(timepoint = "Earlier than mortality")

# SAME FU AS MORTALITY ----------------------------------------------------
UKCAP <-
  filter_surrogate_data(surrogate_data,
                        "Prostate",
                        "1./4. Absolute incidence of late / early stage cancer") %>%
  filter(trial_acronym == "UK CAP") %>%
  filter(str_detect(outcome_description, "10-year")) %>%
  filter(stage_category == ">=8") %>%
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
    stage_category,
    period_after_rnd,
    end_fu_date
  )

PLCO <-
  filter_surrogate_data(surrogate_data,
                        "Prostate",
                        "1./4. Absolute incidence of late / early stage cancer") %>%
  filter(trial_acronym == "PLCO (Prostate)") %>%
  filter(str_detect(outcome_description, "8-10")) %>%
  filter(str_detect(outcome_description, "31 December 2009")) %>%
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
    stage_category,
    period_after_rnd,
    end_fu_date
  )

#Norrkoping  #no estimate available

ERSPC <- surrogate_data_1 %>%
  filter(trial_acronym == "ERSPC") %>%
  # same subgroup as mortality
  filter(subgroup_details == mortality_data$subgroup_details[(mortality_data$trial_acronym ==
                                                                "ERSPC")]) %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "ERSPC")]) %>%
  #late stage definition: gleason 8-10 (reported)
  filter(stage_category == "8-10") %>%
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
    stage_category,
    period_after_rnd,
    end_fu_date
  )

ERSPC_pilot <- surrogate_data_1 %>%
  filter(trial_acronym == "ERSPC Pilot 1") %>%
  # late stage definition: gleason >=7 (3+4)
  filter(stage_category == ">=7 (3+4)") %>%
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
    stage_category,
    period_after_rnd,
    end_fu_date
  )


prostate_surr1 <- rbind(UKCAP, PLCO, ERSPC, ERSPC_pilot) %>%
  calculate_logRR()

main_timepoint <- final_table(mortality_data, prostate_surr1) %>%
  mutate(timepoint = "Mortality")


# ONGOING SCREENING --------------------------------------------------------
UKCAP <-
  filter_surrogate_data(surrogate_data,
                        "Prostate",
                        "1./4. Absolute incidence of late / early stage cancer") %>%
  filter(trial_acronym == "UK CAP") %>%
  filter(stage_category == ">=8") %>%
  filter(period_after_rnd == "<3 years") %>%
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
    stage_category,
    period_after_rnd,
    end_fu_date
  )

#PLCO #no estimate available

#Norrkoping #no estimate available

# ERSPC #no estimate available

#ERSPC_pilot #no estimate available

prostate_surr1 <- UKCAP %>%
  calculate_logrel()

screening_timepoint <-
  final_table(mortality_data, prostate_surr1) %>%
  mutate(timepoint = "Intervention phase")

# COMBINE ALL ESTIMATES ---------------------------------------------------
finaldata_subgroup_1 <-
  rbind(mid_timepoint, main_timepoint, screening_timepoint)

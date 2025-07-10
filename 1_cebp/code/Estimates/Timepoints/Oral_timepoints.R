setwd("../../../code")
source("Allcancer_functions.R")

# Load data ---------------------------------------------------------------
setwd("../data")
mortality_data <- read_excel("SUMS_Results_Master_final.xlsx", sheet="Master_Mortality_Results")
surrogate_data <- read_excel("SUMS_Results_Master_final.xlsx", sheet="Master_Surrogate_Results")

mortality_data <- tibble::as_tibble(mortality_data)
surrogate_data <- tibble::as_tibble(surrogate_data)

setwd("../code/Estimates")
source("Oral.R")

# SAME FU AS MORTALITY ----------------------------------------------------

oral_surr1 <- surrogate_data_1 %>%
  # follow-up time = mortality outcome follow-up
  filter(end_fu_date == "31/12/2004") %>%
  # reported definition of late stage: III + IV
  filter(stage_category == "Stage III-IV") %>%
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


oral_surr1 <- calculate_logRR(oral_surr1)

main_timepoint <- final_table(mortality_data, oral_surr1) %>%
  mutate(timepoint = "Mortality")


# ONGOING SCREENING --------------------------------------------------------
oral_surr1 <- surrogate_data_1 %>%
  # follow-up time = mortality outcome follow-up
  filter(end_fu_date == "30/06/2002") %>%
  # reported definition of late stage: III + IV
  filter(stage_category == "Stage III-IV") %>%
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


oral_surr1 <- calculate_logRR(oral_surr1)

screening_timepoint <- final_table(mortality_data, oral_surr1) %>%
  mutate(timepoint = "Intervention phase")

# COMBINE ALL ESTIMATES ---------------------------------------------------
finaldata_subgroup_1 <- rbind(main_timepoint, screening_timepoint)

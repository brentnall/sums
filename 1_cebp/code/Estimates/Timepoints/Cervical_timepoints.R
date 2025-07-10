setwd("../../../code")
source("Allcancer_functions.R")

# Load data ---------------------------------------------------------------
setwd("../data")
mortality_data <- read_excel("SUMS_Results_Master_final.xlsx", sheet="Master_Mortality_Results")
surrogate_data <- read_excel("SUMS_Results_Master_final.xlsx", sheet="Master_Surrogate_Results")

mortality_data <- tibble::as_tibble(mortality_data)
surrogate_data <- tibble::as_tibble(surrogate_data)

setwd("../code/Estimates")
source("Cervical.R")

# PRIMARY ANALYSIS ESTIMATES ----------------------------------------------

#all are same timpoint as mortality


# SAME FU AS MORTALITY ----------------------------------------------------
main_timepoint <- finaldata_primary_1 %>%
  mutate(timepoint = "Mortality")


# ONGOING SCREENING --------------------------------------------------------
# Finnish #no earlier timepoint

# Mumbai #no earlier timepoint

Osmanabad_Cervix <-
  filter_surrogate_data(surrogate_data,
                        "Cervical",
                        "1./4. Absolute incidence of late / early stage cancer") %>%
  filter(trial_acronym == "Osmanabad Cervix") %>%
  filter(end_fu_date == "November 2003") %>%
  filter(stage_category == "Invasive") %>%
  filter(comparator_noscreening == "Usual care") %>%
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

Tamil_Nadu_Dindigul <-
  filter_surrogate_data(surrogate_data,
                        "Cervical",
                        "1./4. Absolute incidence of late / early stage cancer") %>%
  filter(trial_acronym == "Tamil Nadu/Dindigul") %>%
  #end of screening -2003
  filter(end_fu_date == "37741") %>%
  filter(stage_category == "Stage II+") %>%
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
  ) %>%
  calculate_logRR()

cervical_surr1 <- rbind(Osmanabad_Cervix, Tamil_Nadu_Dindigul)

screening_timepoint <-
  final_table(mortality_data, cervical_surr1) %>%
  mutate(timepoint = "Intervention phase")

# COMBINE ALL ESTIMATES ---------------------------------------------------
finaldata_subgroup_1 <- rbind(main_timepoint, screening_timepoint)


setwd("../../../code")
source("Allcancer_functions.R")

# Load data ---------------------------------------------------------------
setwd("../data")
mortality_data <- read_excel("SUMS_Results_Master_final.xlsx", sheet="Master_Mortality_Results")
surrogate_data <- read_excel("SUMS_Results_Master_final.xlsx", sheet="Master_Surrogate_Results")

mortality_data <- tibble::as_tibble(mortality_data)
surrogate_data <- tibble::as_tibble(surrogate_data)

setwd("../code/Estimates")
source("Ovarian.R")

# PRIMARY ANALYSIS ESTIMATES ----------------------------------------------

#all are same timepoint as mortality

# SAME FU AS MORTALITY ----------------------------------------------------
main_timepoint <- finaldata_primary_1 %>%
  mutate(timepoint = "Mortality")


# ONGOING SCREENING --------------------------------------------------------
UK_Pilot_Ovarian <-
  filter_surrogate_data(surrogate_data,
                        "Ovarian",
                        "1./4. Absolute incidence of late / early stage cancer") %>%
  filter(trial_acronym == "UK Pilot Ovarian")

UKCTOCS <-
  filter_surrogate_data(surrogate_data,
                        "Ovarian",
                        "1./4. Absolute incidence of late / early stage cancer") %>%
  filter(trial_acronym == "UKCTOCS")

PLCO <-
  filter_surrogate_data(surrogate_data,
                        "Ovarian",
                        "1./4. Absolute incidence of late / early stage cancer") %>%
  filter(trial_acronym == "PLCO (Ovarian)") %>%
  filter(period_after_rnd == "0-5") %>%
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

screening_timepoint <- final_table(mortality_data, PLCO) %>%
  mutate(timepoint = "Intervention phase")

# COMBINE ALL ESTIMATES ---------------------------------------------------
finaldata_subgroup_1 <- rbind(main_timepoint, screening_timepoint)

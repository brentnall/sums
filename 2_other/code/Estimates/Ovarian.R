### OVARIAN 
setwd("../../code")
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


# Surrogate 4. -------------------------------------------------------------
surrogate_data_4 <- filter_surrogate_data(surrogate_data, "Ovarian","1./4. Absolute incidence of late / early stage cancer")


#UK Pilot Ovarian
UK_Pilot_Ovarian <- surrogate_data_4 %>%
  filter(trial_acronym == "UK Pilot Ovarian") %>%
  # reported stage category (note, not sure if those numbers were reported or the sum of ind stages was calculated during extraction)
  filter(stage_category == "Stage I-II") %>%
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
UKCTOCS <- surrogate_data_4 %>%
  filter(trial_acronym == "UKCTOCS") %>%
  # use same cancer definition as mortality outcome (=ovarian cancer)
  filter(!str_detect(outcome_description, "ovarian and")) %>%
  # follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[mortality_data$trial_acronym ==
                                                     "UKCTOCS"]) %>%
  # early stage defined as stage I+II (only individual stage categories available for this cutoff)
  filter(stage_category == "Stage I" |
           stage_category == "Stage II") %>%
  filter(str_detect(denominator_definition, "Person-years")) %>%
  group_by(screening_test) %>%
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
  mutate(trial_acronym = "UKCTOCS", stage_category = "Stage I-II") %>%
  calculate_logRR()


#PLCO
PLCO <- surrogate_data_4 %>%
  filter(trial_acronym == "PLCO (Ovarian)") %>%
  # follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[mortality_data$trial_acronym ==
                                                     "PLCO (Ovarian)"]) %>%
  filter(period_after_rnd == mortality_data$period_after_rnd[mortality_data$trial_acronym ==
                                                               "PLCO (Ovarian)"]) %>%
  # late stage defined as stage I+II (only individual stage categories available)
  filter(stage_category == "Stage I" |
           stage_category == "Stage II") %>%
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
  mutate(trial_acronym = "PLCO (Ovarian)", stage_category = "Stage I-II") %>%
  calculate_logRR()

ovarian_surr4 <- rbind(UK_Pilot_Ovarian, UKCTOCS, PLCO)

finaldata_primary_4 <- final_table(mortality_data, ovarian_surr4)

# Surrogate 5. ------------------------------------------------------------
surrogate_data_5 <-
  filter_surrogate_data(surrogate_data,
                        "Ovarian",
                        "5. % of target cancers that are screen-detected")

#UK Pilot Ovarian
UK_Pilot_Ovarian <- surrogate_data_5 %>%
  filter(trial_acronym == "UK Pilot Ovarian") 

#UKCTOCS
UKCTOCS <- surrogate_data_5 %>%
  filter(trial_acronym == "UKCTOCS") %>%
  # use same cancer definition as mortality outcome (=ovarian cancer)
  filter(!str_detect(outcome_description, "ovarian and")) %>%
  filter(!str_detect(outcome_description, "ovarian plus")) 


# one row per screening arm
screening1 <- UKCTOCS %>%
  mutate(comparator_other_screening="Not applicable")

screening2 <- UKCTOCS %>%
  mutate(
    screening_test = comparator_other_screening,
    numerator_screening = numerator_comparator,
    denominator_screening = denominator_comparator
  ) %>%
  mutate(comparator_other_screening="Not applicable")


UKCTOCS <- rbind(screening1, screening2)

#PLCO
PLCO <- surrogate_data_5 %>%
  filter(trial_acronym == "PLCO (Ovarian)") 


ovarian_surr5 <- rbind(UK_Pilot_Ovarian, UKCTOCS, PLCO)

ovarian_surr5 <- ovarian_surr5 %>%
  select(
    trial_acronym,
    screening_test,
    comparator_other_screening,
    outcome_description,
    numerator_screening,
    denominator_screening
  ) %>%
  calculate_prop()

finaldata_primary_5 <- final_table(mortality_data, ovarian_surr5)

# Surrogate 6. ------------------------------------------------------------
surrogate_data_6 <- filter_surrogate_data(surrogate_data, "Ovarian","6. % of high-grade target cancers that are screen-detected")
N <- min_N_trials(surrogate_data_6)
N ### 2 trials

ovarian_surr6 <- surrogate_data_6 %>%
  select(trial_acronym, screening_test, numerator_screening, denominator_screening)

ovarian_surr6 <- calculate_prop(ovarian_surr6)

finaldata_primary_6 <- final_table(mortality_data, ovarian_surr6)

# Surrogate 7. ------------------------------------------------------------
surrogate_data_7 <-
  filter_surrogate_data(surrogate_data, "Ovarian", "7. Diagnostic yield at screening")

#UK Pilot
UK_Pilot_Ovarian <- surrogate_data_7 %>%
  filter(trial_acronym == "UK Pilot Ovarian") %>%
  select(
    trial_acronym,
    screening_test,
    comparator_other_screening,
    outcome_description,
    numerator_screening,
    denominator_screening
  )

#UKCTOCS
UKCTOCS <- surrogate_data_7 %>%
  filter(trial_acronym == "UKCTOCS") %>%
  # use same cancer definition as mortality outcome (=ovarian cancer)
  filter(!str_detect(outcome_description, "ovarian and")) %>%
  filter(!str_detect(outcome_description, "ovarian plus"))

# one row per screening arm
screening1 <- UKCTOCS %>%
  select(
    trial_acronym,
    screening_test,
    comparator_other_screening,
    outcome_description,
    numerator_screening,
    denominator_screening
  )%>%
  mutate(comparator_other_screening = "Not applicable")

screening2 <- UKCTOCS %>%
  select(
    trial_acronym,
    screening_test = comparator_other_screening,
    outcome_description,
    numerator_screening = numerator_comparator,
    denominator_screening = denominator_comparator) %>%
  mutate(comparator_other_screening = "Not applicable")

UKCTOCS <- rbind(screening1, screening2)


#PLCO
PLCO <- surrogate_data_7 %>%
  filter(trial_acronym == "PLCO (Ovarian)") %>%
  select(
    trial_acronym,
    screening_test,
    comparator_other_screening,
    outcome_description,
    numerator_screening,
    denominator_screening
  )

ovarian_surr7 <- rbind(UK_Pilot_Ovarian, UKCTOCS, PLCO)

ovarian_surr7 <- calculate_prop(ovarian_surr7)

finaldata_primary_7 <- final_table(mortality_data, ovarian_surr7)

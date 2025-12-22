### PROSTATE 

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
  filter_mortality_data(mortality_data, cancer = "Prostate")

mortality_data_1 <- mortality_data %>%
  #rate ratio > relative risk/HR, for the studies that report PY
  filter(type_indicator == "Rate ratio") %>%
  calculate_logRR() %>%
  replace_with_reported()

mortality_data_2 <- mortality_data %>%
  #relative risk > adjusted HR (Norrkoping); only relative risk available (ERSPC Pilot)
  filter(type_indicator == "Relative risk") %>%
  # already have rate ratio for UK CAP
  filter(trial_acronym != "UK CAP") %>%
  calculate_logrel() %>%
  replace_with_reported()

mortality_data <- rbind(mortality_data_1, mortality_data_2)


# Surrogate 4. -------------------------------------------------------------
surrogate_data_4 <-
  filter_surrogate_data(surrogate_data,
                        "Prostate",
                        "1./4. Absolute incidence of late / early stage cancer")

UKCAP <- surrogate_data_4 %>%
  filter(trial_acronym == "UK CAP") %>%
  #  6y FU - halfway point between end of screening and mortality
  filter(period_after_rnd == "<6 years") %>%
  # early stage definition: Gleason <=6 + 7 (reported: Localised = T1 or T2)
  filter(stage_category == "<=6" | stage_category == "7") %>%
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
  mutate(trial_acronym = "UK CAP", stage_category = "<=7")

#PLCO
PLCO <- surrogate_data_4 %>%
  filter(trial_acronym == "PLCO (Prostate)") %>%
  #  10y FU
  filter(period_after_rnd == "10") %>%
  # early stage definition: Gleason 2 to 6 + 7
  filter(stage_category == "2-4" |
           stage_category == "5-6" | stage_category == "7") %>%
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
  mutate(trial_acronym = "PLCO (Prostate)", stage_category = "<=7")


#Norrkoping
Norrkoping <- surrogate_data_4 %>%
  filter(trial_acronym == "Norrkoping") %>%
  #late stage definition: T1-2, N0/NX, and M0 (reported as loalised tumour definition)
  filter(stage_category == "T1-2, N0/NX, and M0") %>%
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


#ERSPC
ERSPC <- surrogate_data_4 %>%
  filter(trial_acronym == "ERSPC") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "ERSPC")]) %>%
  # same subgroup as mortality
  filter(subgroup_details == mortality_data$subgroup_details[(mortality_data$trial_acronym ==
                                                                "ERSPC")]) %>%
  # early stage definition: Gleason 2 to 6 + 7
  filter(stage_category == "2-6" | stage_category == "7") %>%
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
  mutate(trial_acronym = "ERSPC", stage_category = "<=7")

#ERSPC pilot 
ERSPC_pilot <- surrogate_data_4 %>%
  filter(trial_acronym == "ERSPC Pilot 1") %>%
  #Gleason score <=6
  filter(stage_category == "6") %>%
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


prostate_surr4 <-
  rbind(UKCAP, PLCO, Norrkoping, ERSPC, ERSPC_pilot) %>%
  calculate_logrel()

finaldata_primary_4 <- final_table(mortality_data, prostate_surr4)


# Surrogate 5. ------------------------------------------------------------
surrogate_data_5 <-
  filter_surrogate_data(surrogate_data,
                        "Prostate",
                        "5. % of target cancers that are screen-detected")

#PLCO
PLCO <- surrogate_data_5 %>%
    filter(trial_acronym == "PLCO (Prostate)") %>%
    filter(`Use for primary analysis`=="yes")
##  filter(period_after_rnd == "Screening period (T0-T5)")

#Norrkoping
Norrkoping <- surrogate_data_5 %>%
  filter(trial_acronym == "Norrkoping")

#ERSPC
ERSPC <- surrogate_data_5 %>%
  filter(trial_acronym == "ERSPC")  %>%
  # same subgroup as mortality
  filter(subgroup_details == mortality_data$subgroup_details[(mortality_data$trial_acronym ==
                                                                "ERSPC")]) %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "ERSPC")]) 

prostate_surr5 <- rbind(PLCO, Norrkoping, ERSPC)

prostate_surr5 <- prostate_surr5 %>%
  select(
    trial_acronym,
    screening_test,
    comparator_other_screening,
    outcome_description,
    numerator_screening,
    denominator_screening,
    reported_rate_screening
  ) %>%
  calculate_prop() %>%
  #reported rate screening only for PLCO - not in update
##  mutate(yi=ifelse(trial_acronym=="PLCO (Prostate)", as.numeric(reported_rate_screening)/100, yi)) %>%
  select(-reported_rate_screening)


finaldata_primary_5 <- final_table(mortality_data, prostate_surr5)


# Surrogate 6. ------------------------------------------------------------
surrogate_data_6 <- filter_surrogate_data(surrogate_data, "Prostate","6. % of high-grade target cancers that are screen-detected")
surrogate_data_6 <- surrogate_data_6 %>%
  #ERSPC only has subgroup estimates
  filter(trial_acronym!="ERSPC")
N <- min_N_trials(surrogate_data_6)
N ### 2 trials

prostate_surr6 <- surrogate_data_6 %>%
  select(trial_acronym, screening_test, numerator_screening, denominator_screening) %>%
  calculate_prop() %>%
  mutate(yi=ifelse(numerator_screening==0, 0, yi))

finaldata_primary_6 <- final_table(mortality_data, prostate_surr6)


# Surrogate 7. ------------------------------------------------------------
surrogate_data_7 <-
  filter_surrogate_data(surrogate_data, "Prostate", "7. Diagnostic yield at screening")

#Norrkoping
Norrkoping <- surrogate_data_7 %>%
  filter(trial_acronym == "Norrkoping") %>%
  # all rounds of screening
  filter(`Effective number of screening rounds` == "4") %>%
  select(
    trial_acronym,
    screening_test,
    comparator_other_screening,
    outcome_description,
    numerator_screening,
    denominator_screening
  )

#PLCO
PLCO <- surrogate_data_7 %>%
  filter(trial_acronym == "PLCO (Prostate)") %>%
  select(
    trial_acronym,
    screening_test,
    comparator_other_screening,
    outcome_description,
    numerator_screening,
    denominator_screening
  )

#ERSPC
ERSPC <- surrogate_data_7 %>%
  filter(trial_acronym == "ERSPC")  %>%
  # same subgroup as mortality
  filter(subgroup_details == mortality_data$subgroup_details[(mortality_data$trial_acronym ==
                                                                "ERSPC")]) %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "ERSPC")]) %>%
  select(
    trial_acronym,
    screening_test,
    comparator_other_screening,
    outcome_description,
    numerator_screening,
    denominator_screening
  )

prostate_surr7 <- rbind(Norrkoping, PLCO, ERSPC)

prostate_surr7 <- calculate_prop(prostate_surr7)

finaldata_primary_7 <- final_table(mortality_data, prostate_surr7)

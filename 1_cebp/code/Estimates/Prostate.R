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


# Surrogate 1. ------------------------------------------------------
surrogate_data_1 <-
  filter_surrogate_data(surrogate_data,
                        "Prostate",
                        "1./4. Absolute incidence of late / early stage cancer")

# UK CAP
UKCAP <- surrogate_data_1 %>%
  filter(trial_acronym == "UK CAP") %>%
  #  6y FU - halfway point between end of screening and mortality
  filter(period_after_rnd == "<6 years") %>%
  # late stage definition: Gleason >=8 (reported: Advanced stage = T4, N1, or M1)
  filter(stage_category == ">=8") %>%
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

#PLCO
PLCO <- surrogate_data_1 %>%
  filter(trial_acronym == "PLCO (Prostate)") %>%
  #  10y FU
  filter(period_after_rnd == "10") %>%
  # late stage definition: Gleason 8-10
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


#Norrkoping
Norrkoping <- surrogate_data_1 %>%
  filter(trial_acronym == "Norrkoping") %>%
  #late stage definition: T3-4, N1, or MX/M1 (reported as advanced tumour definition)
  filter(stage_category == "T3–4, N1 or MX/M1") %>%
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

#ERSPC
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

#ERSPC pilot
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


prostate_surr1 <-
  rbind(UKCAP, PLCO, Norrkoping, ERSPC, ERSPC_pilot) %>%
  calculate_logrel()

finaldata_primary_1 <- final_table(mortality_data, prostate_surr1)


# Surrogate 3. -------------------------------------------------------------
surrogate_data_3 <-
  filter_surrogate_data(surrogate_data,
                        "Prostate",
                        "3. % of target cancer diagnosed at late stage")

#UKCAP
# calculate surrogate 3 using numbers from surrogate 1 file
UKCAP <- surrogate_data_1 %>%
  filter(trial_acronym == "UK CAP") %>%
  #  6y FU - halfway point between end of screening and mortality
  filter(period_after_rnd == "<6 years") %>%
  filter(staging_system == "Gleason score") %>%
  #calculate all cancers diagnosed by arm and use as denominator (all stages + missing)
  mutate(denominator_screening = sum(as.numeric(numerator_screening))) %>%
  mutate(denominator_comparator = sum(as.numeric(numerator_comparator))) %>%
  # late stage definition: Gleason >=8 (reported: Advanced stage = T4, N1, or M1)
  filter(stage_category == ">=8") %>%
  mutate(outcome_description = "%  prostate cancer with Gleason score >=8 at diagnosis (per 100 men);  <6 years of follow-up; Intervention vs control") %>%
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

#ERSPC
# calculate surrogate 3 using numbers from surrogate 1 file
ERSPC <- surrogate_data_1 %>%
  filter(trial_acronym == "ERSPC") %>%
  # same subgroup as mortality
  filter(subgroup_details == mortality_data$subgroup_details[(mortality_data$trial_acronym ==
                                                                "ERSPC")]) %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "ERSPC")]) %>%
  filter(staging_system == "Gleason score") %>%
  #calculate all cancers diagnosed by arm and use as denominator (all stages + missing)
  mutate(denominator_screening = sum(as.numeric(numerator_screening)) +
           as.numeric(missing_screening)) %>%
  mutate(denominator_comparator = sum(as.numeric(numerator_comparator)) +
           as.numeric(missing_comparator)) %>%
  #late stage definition: gleason 8-10 (reported)
  filter(stage_category == "8-10") %>%
  mutate(outcome_description = "% prostate cancers with Gleason score 8-10 at diagnosis in men aged 55-69 years; whole trial (excluding France); up to 31/12/2008") %>%
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


#ERSPC pilot
# calculate surrogate 3 using numbers from surrogate 1 file
ERSPC_pilot <- surrogate_data_1 %>%
  filter(trial_acronym == "ERSPC Pilot 1") %>%
  filter(staging_system == "T stage") %>%
  #calculate all cancers diagnosed by arm and use as denominator (all stages + missing)
  mutate(denominator_screening = sum(as.numeric(numerator_screening)) +
           as.numeric(missing_screening)) %>%
  mutate(denominator_comparator = sum(as.numeric(numerator_comparator)) +
           as.numeric(missing_comparator))

ERSPC_pilot_denominator_s <- ERSPC_pilot$denominator_screening[1]
ERSPC_pilot_denominator_c <- ERSPC_pilot$denominator_comparator[1]

ERSPC_pilot <- surrogate_data_1 %>%
  filter(trial_acronym == "ERSPC Pilot 1") %>%
  # late stage definition: gleason >=7 (3+4)
  filter(stage_category == ">=7 (3+4)") %>%
  mutate(denominator_screening = ERSPC_pilot_denominator_s) %>%
  mutate(denominator_comparator = ERSPC_pilot_denominator_c) %>%
  mutate(outcome_description = "% prostate cancers with Gleason score >=3+4 at diagnosis ; follow-up for median 19 years (IQR 12-24)") %>%
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

#PLCO
PLCO <- surrogate_data_3 %>%
  filter(trial_acronym == "PLCO (Prostate)") %>%
  #  10y FU
  filter(period_after_rnd == "10") %>%
  # late stage definition: Gleason 8-10
  filter(stage_category == "8-10") %>%
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

#Norrkoping
Norrkoping <- surrogate_data_3 %>%
  filter(trial_acronym == "Norrkoping") %>%
  #late stage definition: T3-4, N1, or MX/M1 (reported as advanced tumour definition)
  filter(stage_category == "T3–4, N1 or MX/M1") %>%
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


prostate_surr3 <- rbind(UKCAP, PLCO, Norrkoping, ERSPC, ERSPC_pilot)

prostate_surr3 <- calculate_logrel(prostate_surr3)

finaldata_primary_3 <- final_table(mortality_data, prostate_surr3)

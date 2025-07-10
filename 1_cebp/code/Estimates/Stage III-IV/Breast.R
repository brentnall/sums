### BREAST - Sensitivity (Stage III-IV) Estimates
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

mortality_data_1 <-
  filter_mortality_data(mortality_data, cancer = "Breast") %>%
  filter(trial_acronym != "CNBSS-1" & trial_acronym != "CNBSS-2" & trial_acronym != "UK Age") %>%
  #rate ratio > relative risk/HR, for the studies that report PY
  filter(type_indicator == "Rate ratio") %>%
  calculate_logRR() %>%
  replace_with_reported()

mortality_data_2 <-
  filter_mortality_data(mortality_data, cancer = "Breast") %>%
  # LUSI and the Czech study only report relative risks (so were not included above)
  filter(trial_acronym == "Russia WHO Breast") %>%
  # correct typo in extracted 95%CI Russia WHO breast
  mutate(
    reported_U95CI_unadjusted = ifelse(
      trial_acronym == "Russia WHO Breast",
      "NR",
      reported_U95CI_unadjusted
    )
  ) %>%
  calculate_logrel() %>%
  replace_with_reported()

#Shanghai
mortality_Shanghai <- mortality_data %>%
  filter(trial_acronym == "Shanghai Breast") %>%
  filter(mortality_endpoint == "Cancer-specific") %>%
  #FU up to 2000
  filter(end_fu_date == "36891") %>%
  #rate ratio
  filter(type_indicator == "Rate ratio") %>%
  calculate_logRR() %>%
  replace_with_reported()

#Gothernburg
mortality_Gothenburg <-
  filter_mortality_data(mortality_data, cancer = "Breast") %>%
  filter(trial_acronym == "Gothenburg Breast") %>%
  #cases diagnosed in screening period (blinded review of cause of death); FU until 31/12/1996
  filter(str_detect(outcome_description, "blinded")) %>%
  calculate_logrel() %>%
  replace_with_reported()

#CNBSS-1/2
mortality_CNBSS <- mortality_data %>%
  filter(str_detect(trial_acronym, "CNBSS")) %>%
  filter(mortality_endpoint == "Cancer-specific") %>%
  #7 years of fu
  filter(period_after_rnd == "7 years") %>%
  filter(type_indicator == "Relative risk") %>%
  calculate_logrel() %>%
  replace_with_reported()

#UK Age
mortality_UKAge <- mortality_data %>%
  filter(str_detect(trial_acronym, "UK Age")) %>%
  filter(mortality_endpoint == "Cancer-specific") %>%
  filter(period_after_rnd == "Up to 18.5 years") %>%
  calculate_logRR() %>%
  replace_with_reported()


mortality_data <-
  rbind(
    mortality_data_1,
    mortality_data_2,
    mortality_Shanghai,
    mortality_Gothenburg,
    mortality_CNBSS,
    mortality_UKAge
  )


# Surrogate 1. ------------------------------------------------------
surrogate_data_1 <-
  filter_surrogate_data(surrogate_data,
                        "Breast",
                        "1./4. Absolute incidence of late / early stage cancer")

surrogate_data_1 <- surrogate_data_1 %>%
  # no main mortality endpoint available for CNBSS 1 and 2 combined
  filter(trial_acronym != "CNBSS (1 and 2 combined)")


# Mumbai Breast Cervix
Mumbai <- surrogate_data_1 %>%
  filter(trial_acronym == "Mumbai Breast Cervix") %>%
  #whole trial estimates
  filter(subgroup == "Whole trial") %>%
  # rate ratio > relative risk
  filter(type_indicator == "Rate ratio") %>%
  #late stage definition: Stage III-IV (reported)
  filter(stage_category == "Stage III-IV") %>%
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

#Edinburgh
Edinburgh <- surrogate_data_1 %>%
  filter(trial_acronym == "Edinburgh") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "Edinburgh")]) %>%
  filter(period_after_rnd == mortality_data$period_after_rnd[mortality_data$trial_acronym ==
                                                               "Edinburgh"]) %>%
  # rate ratio > relative risk
  filter(type_indicator == "Rate ratio") %>%
  #late stage definition: stage III-IV (reported)
  filter(stage_category == "III-IV") %>%
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

#Stockholm Breast
Stockholm <- surrogate_data_1 %>%
  filter(trial_acronym == "Stockholm Breast") %>%
  # follow-up time = mortality outcome follow-up (31/12/1986)
  filter(str_detect(outcome_description, "1986")) %>%
  filter(author_year == "Frisell 1997") %>%
  #late stage definition: Stage II-IV (reported)
  filter(stage_category == "III-IV") %>%
  filter(!str_detect(outcome_description, "adjusted")) %>%
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

#Trivandrum Breast
Trivandrum <- surrogate_data_1 %>%
  filter(trial_acronym == "Trivandrum Breast") %>%
  #  follow-up time = mortality outcome follow-up
  filter(str_detect(outcome_description, "follow-up until 31/12/2019")) %>%
  # late stage definition: stage III-IV (reported)
  filter(stage_category == "Stage III-IV") %>%
  # surrogate 1
  filter(str_detect(outcome_description, "Absolute incidence")) %>%
  #rate ratio>relative risk
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

# MMST I
MMSTI <-  surrogate_data_1 %>%
  filter(trial_acronym == "MMST I") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "MMST I")]) %>%
  # late stage definition: stage II-IV (rerorted)
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
  mutate(trial_acronym = "MMST I", stage_category = "Stage III-IV")

breast_surr1 <-
  rbind(Mumbai, Edinburgh, Stockholm, Trivandrum, MMSTI)

breast_surr1_rate <- breast_surr1 %>%
  filter(type_indicator == "Rate ratio") %>%
  calculate_logRR()

breast_surr1_risk <- breast_surr1 %>%
  filter(type_indicator == "Relative risk") %>%
  calculate_logrel()

breast_surr1 <- rbind(breast_surr1_rate, breast_surr1_risk)

finaldata_primary_1 <- final_table(mortality_data, breast_surr1)


# Surrogate 3. -------------------------------------------------------------
surrogate_data_3 <-
  filter_surrogate_data(surrogate_data,
                        "Breast",
                        "3. % of target cancer diagnosed at late stage")

# Mumbai Breast Cervix
Mumbai <- surrogate_data_3 %>%
  filter(trial_acronym == "Mumbai Breast Cervix") %>%
  #whole trial estimates
  filter(subgroup == "Whole trial") %>%
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

#Edinburgh
Edinburgh <- surrogate_data_3 %>%
  filter(trial_acronym == "Edinburgh") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "Edinburgh")]) %>%
  filter(period_after_rnd == mortality_data$period_after_rnd[mortality_data$trial_acronym ==
                                                               "Edinburgh"]) %>%
  #late stage definition: stage III-IV (reported)
  filter(stage_category == "III-IV") %>%
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

#Stockholm Breast
# calculate surrogate 3 using numbers from surrogate 1 file
Stockholm <- surrogate_data_1 %>%
  filter(trial_acronym == "Stockholm Breast") %>%
  # follow-up time = mortality outcome follow-up (31/12/1986)
  filter(str_detect(outcome_description, "1986")) %>%
  filter(author_year == "Frisell 1997") %>%
  filter(stage_category != "II-IV") %>%
  #not include CIS for consistency with other trials
  filter(stage_category != "CIS") %>%
  #calculate all cancers diagnosed by arm and use as denominator (all stages + missing)
  mutate(denominator_screening = sum(as.numeric(numerator_screening))) %>%
  mutate(denominator_comparator = sum(as.numeric(numerator_comparator))) %>%
  # late stage defiiniion: node positive (to be consistent with other trials included)
  filter(stage_category == "III-IV") %>%
  mutate(outcome_description = "% stage III-IV, invasive breast cancers (1981-1986); per 10,000 women; screening vs control group") %>%
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


#Trivandrum Breast
Trivandrum <- surrogate_data_3 %>%
  filter(trial_acronym == "Trivandrum Breast") %>%
  #FU up to 2019
  filter(end_fu_date == "43830") %>%
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

# MMST I
# calculate surrogate 3 using numbers from surrogate 1 file
MMSTI <-  surrogate_data_1 %>%
  filter(trial_acronym == "MMST I") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "MMST I")]) %>%
  filter(stage_category != "II-IV") %>%
  #not include stage 0 for consistency with other trials
  filter(stage_category != "0") %>%
  #calculate all cancers diagnosed by arm and use as denominator (all stages + missing)
  mutate(denominator_screening = sum(as.numeric(numerator_screening)) +
           as.numeric(missing_screening)) %>%
  mutate(denominator_comparator = sum(as.numeric(numerator_comparator)) +
           as.numeric(missing_comparator))

MMSTI_denominator_s <- MMSTI$denominator_screening[1]
MMSTI_denominator_c <- MMSTI$denominator_comparator[1]

MMSTI <-  surrogate_data_1 %>%
  filter(trial_acronym == "MMST I") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "MMST I")]) %>%
  # late stage definition: stage II-IV (rerorted)
  filter(stage_category == "III" | stage_category == "IV") %>%
  summarise(
    screening_test = first(screening_test),
    numerator_screening = sum(as.numeric(numerator_screening)),
    numerator_comparator = sum(as.numeric(numerator_comparator)),
    denominator_screening = MMSTI_denominator_s,
    denominator_comparator = MMSTI_denominator_c,
    comparator_other_screening = first(comparator_other_screening),
    comparator_noscreening = first(comparator_noscreening),
    reported_definition_late = first(reported_definition_late),
    reported_definition_early = first(reported_definition_early)
  ) %>%
  mutate(
    trial_acronym = "MMST I",
    stage_category = "Stage III-IV",
    outcome_description = "% stage III (UICC) breast cancer; follow-up to 31/12/1986; Screening vs control"
  )


breast_surr3 <-
  rbind(Mumbai, Edinburgh, Stockholm, Trivandrum, MMSTI)

breast_surr3 <- calculate_logrel(breast_surr3)

finaldata_primary_3 <- final_table(mortality_data, breast_surr3)

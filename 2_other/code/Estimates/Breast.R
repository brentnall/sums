### BREAST 
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


# Surrogate 4. -------------------------------------------------------------
surrogate_data_4 <-
  filter_surrogate_data(surrogate_data,
                        "Breast",
                        "1./4. Absolute incidence of late / early stage cancer")
surrogate_data_4 <- surrogate_data_4 %>%
  # no main mortality endpoint available for CNBSS (1 and 2 combined)
  filter(trial_acronym != "CNBSS (1 and 2 combined)")


# Mumbai Breast Cervix
Mumbai <- surrogate_data_4 %>%
  filter(trial_acronym == "Mumbai Breast Cervix") %>%
  #whole trial estimates
  filter(subgroup == "Whole trial") %>%
  # rate ratio > relative risk
  filter(type_indicator == "Rate ratio") %>%
  #early stage definition: Stage I-II (reported)
  filter(stage_category == "Stage I-II") %>%
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

# UK Age
UKAge <- surrogate_data_4 %>%
  filter(trial_acronym == "UK Age") %>%
  # follow-up time = mortality outcome follow-up not available, closest time point used
  filter(period_after_rnd == "Intervention period") %>%
  # early stage definition: lymph-node negative
  filter(stage_category == "Node negative") %>%
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
Edinburgh <- surrogate_data_4 %>%
  filter(trial_acronym == "Edinburgh") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "Edinburgh")]) %>%
  filter(period_after_rnd == mortality_data$period_after_rnd[mortality_data$trial_acronym ==
                                                               "Edinburgh"]) %>%
  # rate ratio > relative risk
  filter(type_indicator == "Rate ratio") %>%
  #early stage definition: Stage 0 - II (stage III-IV reported as late)
  filter(
    stage_category == "0 (non-invasive disease)" |
      stage_category == "I" | stage_category == "II"
  ) %>%
  summarise(
    screening_test = first(screening_test),
    outcome_description = last(outcome_description),
    numerator_screening = sum(as.numeric(numerator_screening)),
    numerator_comparator = sum(as.numeric(numerator_comparator)),
    denominator_screening = mean(as.numeric(denominator_screening)),
    denominator_comparator = mean(as.numeric(denominator_comparator)),
    type_indicator = first(type_indicator),
    reported_definition_late = first(reported_definition_late),
    reported_definition_early = first(reported_definition_early)
  ) %>%
  mutate(trial_acronym = "Edinburgh", stage_category = "Stage 0-II")

#CNBSS-1
CNBSS1 <- surrogate_data_4 %>%
  filter(trial_acronym == "CNBSS-1") %>%
  #  7 years FU
  filter(period_after_rnd == "7 years") %>%
  # late stage defiiniion: node negative (to be consistent with other trials included)
  filter(stage_category == "Node-negative") %>%
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

#CNBSS-2
CNBSS2 <- surrogate_data_4 %>%
  filter(trial_acronym == "CNBSS-2") %>%
  #  7 years FU
  filter(period_after_rnd == "7 years") %>%
  # late stage defiiniion: node negative (to be consistent with other trials included)
  filter(stage_category == "Node-negative") %>%
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
Stockholm <- surrogate_data_4 %>%
  filter(trial_acronym == "Stockholm Breast") %>%
  # follow-up time = mortality outcome follow-up (31/12/1986)
  filter(str_detect(outcome_description, "1986")) %>%
  #first publication
  filter(author_year == "Frisell 1997") %>%
  #early stage definition: Stage I (reported) + CIS
  filter(stage_category == "CIS" | stage_category == "I") %>%
  summarise(
    screening_test = first(screening_test),
    outcome_description = last(outcome_description),
    numerator_screening = sum(as.numeric(numerator_screening)),
    numerator_comparator = sum(as.numeric(numerator_comparator)),
    denominator_screening = mean(as.numeric(denominator_screening)),
    denominator_comparator = mean(as.numeric(denominator_comparator)),
    type_indicator = first(type_indicator),
    reported_definition_late = first(reported_definition_late),
    reported_definition_early = first(reported_definition_early)
  ) %>%
  mutate(trial_acronym = "Stockholm Breast", stage_category = "CIS + Stage I")

#Trivandrum Breast
Trivandrum <- surrogate_data_4 %>%
  filter(trial_acronym == "Trivandrum Breast") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "Trivandrum Breast")]) %>%
  # early stage definition: stage I-II (reported)
  filter(stage_category == "Stage I-II") %>%
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

#Two-County
Twocounty <- surrogate_data_4 %>%
  filter(trial_acronym == "Two-County") %>%
  # same subgroup as mortality
  filter(subgroup_details == mortality_data$subgroup_details[(mortality_data$trial_acronym ==
                                                                "Two-County")]) %>%
  #  FU until 31/12/1984 (screening period)
  filter(end_fu_date == "31047") %>%
  filter(str_detect(outcome_description, "1984")) %>%
  #early stage definition: Stage I (reported)
  filter(stage_category == "I" |
           stage_category == "Ductal in situ") %>%
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
  mutate(trial_acronym = "Two-County", stage_category = "DCIS/Stage I")

#Shanghai Breast
Shanghai <-  surrogate_data_4 %>%
  filter(trial_acronym == "Shanghai Breast") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == "36891") %>%
  # early stage definition: Tis (in situ, stage 0) or T1 (<=2 cm) (reported)
  filter(stage_category == "Tis (in situ, stage 0) or T1 (<=2 cm)") %>%
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

#Gothenburg Breast
Gothenburg <-  surrogate_data_4 %>%
  filter(trial_acronym == "Gothenburg Breast") %>%
  # whole trial estimate
  filter(subgroup == "Whole trial") %>%
  filter(author_year == "Bjurstam 2016") %>%
  # early stage definition: node negative
  filter(stage_category == "Node-negative") %>%
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
MMSTI <-  surrogate_data_4 %>%
  filter(trial_acronym == "MMST I") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym == "MMST I")]) %>%
  # early stage definition: stage 0-I (stage II-IV rerorted as late)
  filter(stage_category == "0" | stage_category == "I") %>%
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
  mutate(trial_acronym = "MMST I", stage_category = "Stage 0-I")

# Russia WHO Breast
RussiaWHO <-  surrogate_data_4 %>%
  filter(trial_acronym == "Russia WHO Breast") %>%
  # same f-u as mortality not available, closest is 31/12/1994
  filter(end_fu_date == "34699") %>%
  # early stage definition: Node negative (N0)
  filter(stage_category == "Node negative (N0)") %>%
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

# HIP
HIP <-  surrogate_data_4 %>%
  filter(trial_acronym == "HIP") %>%
  # follow-up: 5 years since rnd
  filter(period_after_rnd == "5 years") %>%
  # early stage definition: node negative
  filter(stage_category == "Node-negative") %>%
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


breast_surr4 <- rbind(Mumbai, UKAge, Edinburgh, CNBSS1, CNBSS2, Stockholm, Trivandrum, Twocounty, Shanghai, Gothenburg, MMSTI, RussiaWHO, HIP)

breast_surr4_rate <- breast_surr4 %>%
  filter(type_indicator=="Rate ratio") %>%
  calculate_logRR()

breast_surr4_risk <- breast_surr4 %>%
  filter(type_indicator=="Relative risk") %>%
  calculate_logrel()

breast_surr4 <- rbind(breast_surr4_rate, breast_surr4_risk)

finaldata_primary_4 <- final_table(mortality_data, breast_surr4)


# Surrogate 5. ------------------------------------------------------------
surrogate_data_5 <-
  filter_surrogate_data(surrogate_data,
                        "Breast",
                        "5. % of target cancers that are screen-detected")

# Mumbai Breast Cervix
Mumbai <- surrogate_data_5 %>%
  filter(trial_acronym == "Mumbai Breast Cervix") 

# UK Age
UKAge <- surrogate_data_5 %>%
  filter(trial_acronym == "UK Age") %>%
  # follow-up time = mortality outcome follow-up not available, closest time point used
  filter(period_after_rnd == "Intervention period") %>%
  # same stimates (diff numbers) use same publication as surrogates 1-4
  filter(author_year == "Duffy 2020") 

#Edinburgh
Edinburgh <- surrogate_data_5 %>%
  filter(trial_acronym == "Edinburgh") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "Edinburgh")]) %>%
  filter(period_after_rnd == mortality_data$period_after_rnd[mortality_data$trial_acronym ==
                                                               "Edinburgh"]) 


#CNBSS-1
CNBSS1 <- surrogate_data_5 %>%
  filter(trial_acronym == "CNBSS-1") %>%
  filter(period_after_rnd == "7 years") 


#CNBSS-2
CNBSS2 <- surrogate_data_5 %>%
  filter(trial_acronym == "CNBSS-2") 

#Stockholm Breast
Stockholm <- surrogate_data_5 %>%
  filter(trial_acronym == "Stockholm Breast") %>%
  # follow-up time = mortality outcome follow-up (31/12/1986)
  filter(str_detect(outcome_description, "1986")) %>%
  filter(author_year == "Frisell 1997") %>%
  #include invasive + CIS
  filter(str_detect(numerator_definition, "CIS")) 


#Trivandrum Breast
Trivandrum <- surrogate_data_5 %>%
  filter(trial_acronym == "Trivandrum Breast") %>%
  #  follow-up time = mortality outcome follow-up // all 3 rounds
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "Trivandrum Breast")]) 


#Two-County
Twocounty <- surrogate_data_5 %>%
  filter(trial_acronym == "Two-County") %>%
  # same subgroup as mortality
  filter(subgroup_details == "40-74 years") %>%
  # first publication
  filter(author_year == "Tabar 1992")


#Shanghai Breast
Shanghai <-  surrogate_data_5 %>%
  filter(trial_acronym == "Shanghai Breast") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == "36891") 


#Gothenburg Breast
Gothenburg <-  surrogate_data_5 %>%
  filter(trial_acronym == "Gothenburg Breast") %>%
  # whole trial estimate
  filter(subgroup == "Whole trial") %>%
  # invasive + DCIS
  filter(str_detect(outcome_description, "DCIS")) %>%
  #first publication
  filter(author_year == "Bjurstam 2003") 


# MMST I
MMSTI <-  surrogate_data_5 %>%
  filter(trial_acronym == "MMST I") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "MMST I")]) 


# HIP
HIP <-  surrogate_data_5 %>%
  filter(trial_acronym == "HIP") %>%
  # whole trial estimates
  filter(subgroup == "Whole trial") %>%
  # follow-up: 5 years since rnd
  filter(period_after_rnd == "5 years") %>%
  ## first publication
  filter(author_year == "Shapiro 1977") 

breast_surr5 <-
  rbind(
    Mumbai,
    UKAge,
    Edinburgh,
    CNBSS1,
    CNBSS2,
    Stockholm,
    Trivandrum,
    Twocounty,
    Shanghai,
    Gothenburg,
    MMSTI,
    HIP
  )

breast_surr5 <- breast_surr5 %>%
  select(
    trial_acronym,
    screening_test,
    comparator_other_screening,
    outcome_description,
    numerator_screening,
    denominator_screening
  ) %>%
  calculate_prop()

finaldata_primary_5 <- final_table(mortality_data, breast_surr5)

# Surrogate 6. ------------------------------------------------------------
surrogate_data_6 <- filter_surrogate_data(surrogate_data, "Breast","6. % of high-grade target cancers that are screen-detected")
N <- min_N_trials(surrogate_data_6)
N ### 4 trials (3, two-county estimate is for different subgroup than mortality)

#UK Age
UKAge <- surrogate_data_6 %>%
  filter(trial_acronym=="UK Age") %>%
  # intervention period has ended, first publication from those
  filter(author_year=="Moss 2015") %>%
  select(trial_acronym, screening_test, numerator_screening, denominator_screening)
  

#Kopparberg County - no main mortality outcome?

#Ostergotland County - no main mortality outcome?

#breast_surr6 <- rbind(UKAge, Kopparberg, Ostergotland)

breast_surr6 <- calculate_prop(UKAge)

finaldata_primary_6 <- final_table(mortality_data, breast_surr6)

# Surrogate 7. ------------------------------------------------------------
surrogate_data_7 <-
  filter_surrogate_data(surrogate_data, "Breast", "7. Diagnostic yield at screening")

# Mumbai Breast Cervix
Mumbai <- surrogate_data_7 %>%
  filter(trial_acronym == "Mumbai Breast Cervix") %>%
  select(
    trial_acronym,
    screening_test,
    comparator_other_screening,
    outcome_description,
    numerator_screening,
    denominator_screening
  )

# UK Age
UKAge <- surrogate_data_7 %>%
  filter(trial_acronym == "UK Age") %>%
  # after completion of screening
  filter(period_after_rnd == "Intervention period") %>%
  select(
    trial_acronym,
    screening_test,
    comparator_other_screening,
    outcome_description,
    numerator_screening,
    denominator_screening
  )


#Edinburgh
Edinburgh <- surrogate_data_7 %>%
  filter(trial_acronym == "Edinburgh") %>%
  #subgroup same as mortality (cohort 1)
  filter(subgroup_details == mortality_data$subgroup_details[(mortality_data$trial_acronym ==
                                                                "Edinburgh")]) %>%
  select(
    trial_acronym,
    screening_test,
    comparator_other_screening,
    outcome_description,
    numerator_screening,
    denominator_screening
  )


#CNBSS-1
CNBSS1 <- surrogate_data_7 %>%
  filter(trial_acronym == "CNBSS-1") %>%
  filter(end_fu_date == "5 screening rounds") %>%
  select(
    trial_acronym,
    screening_test,
    outcome_description,
    comparator_other_screening,
    numerator_screening,
    denominator_screening
  )

#CNBSS-2
CNBSS2 <- surrogate_data_7 %>%
  filter(trial_acronym == "CNBSS-2") %>%
  select(
    trial_acronym,
    screening_test,
    comparator_other_screening,
    outcome_description,
    numerator_screening,
    denominator_screening
  )


#Stockholm Breast
Stockholm <- surrogate_data_7 %>%
  filter(trial_acronym == "Stockholm Breast") %>%
  # whole trial estimate
  filter(subgroup == "Whole trial") %>%
  # all screening rounds
  filter(`Effective number of screening rounds` == "2") %>%
  #CIS + invasive breast
  filter(str_detect(outcome_description, "CIS")) %>%
  #first publication
  filter(author_year == "Frisell 1989") %>%
  select(
    trial_acronym,
    screening_test,
    comparator_other_screening,
    outcome_description,
    numerator_screening,
    denominator_screening
  )


#Trivandrum Breast
Trivandrum <- surrogate_data_7 %>%
  filter(trial_acronym == "Trivandrum Breast") %>%
  # all 3 screening rounds
  filter(`Effective number of screening rounds` == "3") %>%
  select(
    trial_acronym,
    screening_test,
    comparator_other_screening,
    outcome_description,
    numerator_screening,
    denominator_screening
  )


#Two-County
Twocounty <- surrogate_data_7 %>%
  filter(trial_acronym == "Two-County") %>%
  # same subgroup as mortality
  filter(subgroup_details == mortality_data$subgroup_details[(mortality_data$trial_acronym ==
                                                                "Two-County")]) %>%
  #first publication
  filter(author_year == "Tabar 1992") %>%
  select(
    trial_acronym,
    screening_test,
    comparator_other_screening,
    outcome_description,
    numerator_screening,
    denominator_screening
  )


#Shanghai Breast
Shanghai <-  surrogate_data_7 %>%
  filter(trial_acronym == "Shanghai Breast") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == "34699") %>%
  select(
    trial_acronym,
    screening_test,
    comparator_other_screening,
    outcome_description,
    numerator_screening,
    denominator_screening
  )


#Gothenburg Breast
Gothenburg <-  surrogate_data_7 %>%
  filter(trial_acronym == "Gothenburg Breast") %>%
  # whole trial estimate
  filter(subgroup == "Whole trial") %>%
  # invasive + DCIS
  filter(str_detect(outcome_description, "DCIS")) %>%
  #first publication
  filter(author_year == "Bjurstam 2003") %>%
  select(
    trial_acronym,
    screening_test,
    outcome_description,
    comparator_other_screening,
    numerator_screening,
    denominator_screening
  )

# MMST I
MMSTI <-  surrogate_data_7 %>%
  filter(trial_acronym == "MMST I") %>%
  #  follow-up time = mortality outcome follow-up; same number of screening rounds
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "MMST I")]) %>%
  select(
    trial_acronym,
    screening_test,
    comparator_other_screening,
    outcome_description,
    numerator_screening,
    denominator_screening
  )

# HIP
HIP <-  surrogate_data_7 %>%
  filter(trial_acronym == "HIP") %>%
  # whole trial estimates
  filter(subgroup == "Whole trial") %>%
  # follow-up: 5 years since rnd (4 rounds)
  filter(period_after_rnd == "5 years") %>%
  #first publication
  filter(author_year == "Shapiro 1985") %>%
  select(
    trial_acronym,
    screening_test,
    comparator_other_screening,
    outcome_description,
    numerator_screening,
    denominator_screening
  )

#Russia WHO
RussiaWHO <-  surrogate_data_7 %>%
  filter(trial_acronym == "Russia WHO Breast") %>%
  #whole intervention phase
  filter(end_fu_date == "33238") %>%
  select(
    trial_acronym,
    screening_test,
    comparator_other_screening,
    outcome_description,
    numerator_screening,
    denominator_screening
  )


breast_surr7 <-
  rbind(
    Mumbai,
    UKAge,
    Edinburgh,
    CNBSS1,
    CNBSS2,
    Stockholm,
    Trivandrum,
    Twocounty,
    Shanghai,
    Gothenburg,
    MMSTI,
    RussiaWHO,
    HIP
  )

breast_surr7 <- calculate_prop(breast_surr7)

finaldata_primary_7 <- final_table(mortality_data, breast_surr7)


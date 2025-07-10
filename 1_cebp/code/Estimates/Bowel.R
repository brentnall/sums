### BOWEL
setwd("../../code")
source("Allcancer_functions.R")

# Load data ---------------------------------------------------------------
setwd("../data")
mortality_data <- read_excel("SUMS_Results_Master_final.xlsx", sheet="Master_Mortality_Results")
surrogate_data <- read_excel("SUMS_Results_Master_final.xlsx", sheet="Master_Surrogate_Results") %>%
  #rename for consistency
  mutate(trial_acronym=ifelse(trial_acronym=="PLCO (Bowel)", "PLCO (CRC)", trial_acronym))

mortality_data <- tibble::as_tibble(mortality_data)
surrogate_data <- tibble::as_tibble(surrogate_data)


# PRIMARY ANALYSIS - Estimate Selection -----------------------------------

# Mortality  ---------------------------------------------------------------

mortality_data1 <-
  filter_mortality_data(mortality_data, cancer = "Bowel") %>%
  filter(type_indicator == "Rate ratio") %>%
  filter(trial_acronym != "NORCCAP")

NordICC <-  filter_mortality_data(mortality_data, cancer = "Bowel") %>%
  filter(trial_acronym == "NordICC") %>%
  filter(str_detect(outcome_description, "Kaplan-Meier estimate"))

UKFSST <- filter_mortality_data(mortality_data, cancer = "Bowel") %>%
  filter(trial_acronym == "UKFSST") %>%
  filter(type_indicator == "Hazard ratio") %>%
  # change to correspond to estimate used
  mutate(type_indicator = "Rate ratio")

NORCCAP <- mortality_data %>%
  filter(trial_acronym == "NORCCAP") %>%
  filter(mortality_endpoint == "Cancer-specific") %>%
  #FU 15 years
  filter(end_fu_date == "42369") %>%
  filter(subgroup == "Age group") %>%
  filter(type_indicator == "Hazard ratio") %>%
  mutate(type_indicator = "Rate ratio")

Telemark <- filter_mortality_data(mortality_data, cancer = "Bowel") %>%
  filter(trial_acronym == "Telemark Polyp I")

mortality_data <-
  rbind(mortality_data1, NordICC, UKFSST, NORCCAP, Telemark)

# calculate effect estimates or replace with reported if available
mortality_data_rate <- mortality_data %>%
  filter(type_indicator == "Rate ratio") %>%
  calculate_logRR() %>%
  replace_with_reported()

mortality_data_risk <- mortality_data %>%
  filter(type_indicator == "Relative risk") %>%
  calculate_logrel() %>%
  replace_with_reported()


mortality_data <- rbind(mortality_data_rate, mortality_data_risk)




# Surrogate 1. ------------------------------------------------------
surrogate_data_1 <-
  filter_surrogate_data(surrogate_data,
                        "Bowel",
                        "1./4. Absolute incidence of late / early stage cancer")

#NordICC
NordICC <- surrogate_data_1 %>%
  filter(trial_acronym == "NordICC") %>%
  #late stage definition: Dukes’ stage C or D (reported)
  filter(stage_category == "Dukes' C or D") %>%
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

#SCORE
SCORE <- surrogate_data_1 %>%
  filter(trial_acronym == "SCORE") %>%
  # same definition of colorectal cancer as mortality (distal + proximal)
  filter(str_detect(outcome_description, "distal")) %>%
  filter(str_detect(outcome_description, "proximal")) %>%
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
    type_indicator,
    reported_definition_late,
    reported_definition_early,
    stage_category
  )

#UKFSST
UKFSST <- surrogate_data_1 %>%
  filter(trial_acronym == "UKFSST") %>%
  #whole trial estimate
  filter(subgroup == "Whole trial") %>%
  # same definition of colorectal cancer as mortality (distal + proximal)
  filter(str_detect(outcome_description, "all sites")) %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "UKFSST")]) %>%
  # rate ratio > relative risk
  filter(type_indicator == "Hazard ratio") %>%
  mutate(type_indicator = "Rate ratio") %>%
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

#Finnish RHS CRC
FinnishRHS <- surrogate_data_1 %>%
  filter(trial_acronym == "Finnish RHS CRC") %>%
  #whole trial estimate
  filter(subgroup == "Whole trial") %>%
  #  follow-up time = mortality outcome follow-up not available, closest is 31/12/2011
  filter(end_fu_date == "40908") %>%
  # late stage definition: node positive (N>=1)
  filter(stage_category == "N>=1") %>%
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

#Funen
Funen <- surrogate_data_1 %>%
  filter(trial_acronym == "Funen") %>%
  #  follow-up time = mortality outcome follow-up (August 1995)
  filter(end_fu_date == "34912") %>%
  #late stage definition: stage C, distant spread, and no classification
  filter(stage_category == "C, distant sprad or no classification") %>%
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

#Gothenburg CRC
Gothenburg <- surrogate_data_1 %>%
  filter(trial_acronym == "Gothenburg CRC") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "Gothenburg CRC")]) %>%
  filter(period_after_rnd == mortality_data$period_after_rnd[mortality_data$trial_acronym ==
                                                               "Gothenburg CRC"]) %>%
  # late stage definition: Duke's D (reported in previous publications)
  filter(stage_category == "D") %>%
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
    type_indicator,
    reported_definition_late,
    reported_definition_early,
    stage_category
  )

#Minnesota
Minnesota <- surrogate_data_1 %>%
  filter(trial_acronym == "Minnesota") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "Minnesota")]) %>%
  filter(period_after_rnd == mortality_data$period_after_rnd[mortality_data$trial_acronym ==
                                                               "Minnesota"]) %>%
  # late stage definition: Duke's D
  filter(stage_category == "D") %>%
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
    type_indicator,
    reported_definition_late,
    reported_definition_early,
    stage_category
  )

#NORCCAP
NORCCAP <- surrogate_data_1 %>%
  filter(trial_acronym == "NORCCAP") %>%
  # late stage definition: Dukes' C or D (reported)
  filter(str_detect(stage_category, "Dukes' C")) %>%
  #FU 31/12/2006 (median 7 years)
  filter(end_fu_date == "39082") %>%
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

#Nottingham
Nottingham <- surrogate_data_1 %>%
  filter(trial_acronym == "Nottingham") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "Nottingham")]) %>%
  filter(period_after_rnd == mortality_data$period_after_rnd[mortality_data$trial_acronym ==
                                                               "Nottingham"]) %>%
  # late stage definition: Dukes' C or D (reported)
  filter(stage_category == "C or D") %>%
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
    type_indicator,
    reported_definition_late,
    reported_definition_early,
    stage_category
  )

#	Telemark Polyp I
Telemark <- surrogate_data_1 %>%
  filter(trial_acronym == "Telemark Polyp I") %>%
  #  follow-up time until 1993
  filter(end_fu_date == "1993") %>%
  # late stage definition: Dukes' C or D
  filter(stage_category == "C") %>%
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

#PLCO (CRC)
PLCO <- surrogate_data_1 %>%
  filter(trial_acronym == "PLCO (CRC)") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "PLCO (CRC)")]) %>%
  # same definition of colorectal cancer as mortality (distal + proximal)
  filter(!str_detect(outcome_description, "Distal")) %>%
  filter(!str_detect(outcome_description, "Proximal")) %>%
  # late stage definition: stage III/IV
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
  mutate(trial_acronym = "PLCO (CRC)", stage_category = "Stage III-IV")


crc_surr1 <-
  rbind(
    NordICC,
    SCORE,
    FinnishRHS,
    Funen,
    Gothenburg,
    Minnesota,
    NORCCAP,
    Nottingham,
    Telemark,
    PLCO,
    UKFSST
  )

crc_surr1_rate <- crc_surr1 %>%
  filter(type_indicator == "Rate ratio") %>%
  calculate_logRR()

crc_surr1_risk <- crc_surr1 %>%
  filter(type_indicator == "Relative risk") %>%
  calculate_logrel()

crc_surr1 <- rbind(crc_surr1_rate, crc_surr1_risk)

finaldata_primary_1 <- final_table(mortality_data, crc_surr1)


# Surrogate 3. -------------------------------------------------------------
surrogate_data_3 <-
  filter_surrogate_data(surrogate_data,
                        "Bowel",
                        "3. % of target cancer diagnosed at late stage")

#NordICC
NordICC <- surrogate_data_3 %>%
  filter(trial_acronym == "NordICC") %>%
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

#SCORE
SCORE <- surrogate_data_3 %>%
  filter(trial_acronym == "SCORE") %>%
  # same definition of colorectal cancer as mortality (distal + proximal)
  filter(str_detect(outcome_description, "distal")) %>%
  filter(str_detect(outcome_description, "proximal")) %>%
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

#Finnish RHS
# calculate surrogate 3 using numbers from surrogate 1 file
FinnishRHS <- surrogate_data_1 %>%
  filter(trial_acronym == "Finnish RHS CRC") %>%
  #whole trial estimate
  filter(subgroup == "Whole trial") %>%
  #  follow-up time = mortality outcome follow-up not available, closest is 31/12/2011
  filter(end_fu_date == "40908") %>%
  #  node status
  filter(staging_system == "Node status") %>%
  #calculate all cancers diagnosed by arm and use as denominator (all stages + missing)
  mutate(denominator_screening = sum(as.numeric(numerator_screening)) +
           as.numeric(missing_screening)) %>%
  mutate(denominator_comparator = sum(as.numeric(numerator_comparator)) +
           as.numeric(missing_comparator)) %>%
  #late stage definition: N>=1
  filter(stage_category == "N>=1") %>%
  mutate(outcome_description = "% CRC with N>=1 (men + women; left + right colon); follow-up until 31/12/2011") %>%
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

#Funen
Funen <- surrogate_data_3 %>%
  filter(trial_acronym == "Funen") %>%
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

#Gotthenburg
Gothenburg <- surrogate_data_3 %>%
  filter(trial_acronym == "Gothenburg CRC") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "Gothenburg CRC")]) %>%
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

#NORCCAP
NORCCAP <- surrogate_data_3 %>%
  filter(trial_acronym == "NORCCAP") %>%
  #FU 31/12/2006 (median 7 years)
  filter(end_fu_date == "39082") %>%
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

#Nottingham
Nottingham <- surrogate_data_3 %>%
  filter(trial_acronym == "Nottingham") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "Nottingham")]) %>%
  filter(period_after_rnd == mortality_data$period_after_rnd[mortality_data$trial_acronym ==
                                                               "Nottingham"]) %>%
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
  filter(trial_acronym == "PLCO (CRC)") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "PLCO (CRC)")]) %>%
  # late stage definition: stage III/IV
  filter(stage_category == "III" | stage_category == "IV") %>%
  summarise(
    screening_test = first(screening_test),
    outcome_description = first(outcome_description),
    numerator_screening = sum(as.numeric(numerator_screening)),
    numerator_comparator = sum(as.numeric(numerator_comparator)),
    denominator_screening = mean(as.numeric(denominator_screening)),
    denominator_comparator = mean(as.numeric(denominator_comparator)),
    reported_definition_late = first(reported_definition_late),
    reported_definition_early = first(reported_definition_early),
    comparator_other_screening = first(comparator_other_screening),
    comparator_noscreening = first(comparator_noscreening)
  ) %>%
  mutate(trial_acronym = "PLCO (CRC)", stage_category = "Stage III-IV")

#Minnesota
# calculate surrogate 3 using numbers from surrogate 1 file
Minnesota_annual <- surrogate_data_1 %>%
  filter(trial_acronym == "Minnesota") %>%
  filter(str_detect(outcome_description, "Annual")) %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "Minnesota")]) %>%
  filter(period_after_rnd == mortality_data$period_after_rnd[mortality_data$trial_acronym ==
                                                               "Minnesota"]) %>%
  #rate ratio > relative risk
  filter(type_indicator == "Rate ratio") %>%
  #calculate all cancers diagnosed by arm and use as denominator (all stages + missing)
  mutate(denominator_screening = sum(as.numeric(numerator_screening)) +
           as.numeric(missing_screening)) %>%
  mutate(denominator_comparator = sum(as.numeric(numerator_comparator)) +
           as.numeric(missing_comparator)) %>%
  # late stage defiiniion: node positive (to be consistent with other trials included)
  filter(stage_category == "D") %>%
  mutate(outcome_description = "% Dukes' D CRC (per 1,000 PY); 13 years of follow-up; Annual FOBT vs contro") %>%
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

Minnesota_biennial <- surrogate_data_1 %>%
  filter(trial_acronym == "Minnesota") %>%
  filter(str_detect(outcome_description, "Biennial")) %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "Minnesota")]) %>%
  filter(period_after_rnd == mortality_data$period_after_rnd[mortality_data$trial_acronym ==
                                                               "Minnesota"]) %>%
  #rate ratio > relative risk
  filter(type_indicator == "Rate ratio") %>%
  #calculate all cancers diagnosed by arm and use as denominator (all stages + missing)
  mutate(denominator_screening = sum(as.numeric(numerator_screening)) +
           as.numeric(missing_screening)) %>%
  mutate(denominator_comparator = sum(as.numeric(numerator_comparator)) +
           as.numeric(missing_comparator)) %>%
  # late stage defiiniion: node positive (to be consistent with other trials included)
  filter(stage_category == "D") %>%
  mutate(outcome_description = "% Dukes' D CRC (per 1,000 PY); 13 years of follow-up; Biennial FOBT vs contro") %>%
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

#Telemark
# calculate surrogate 3 using numbers from surrogate 1 file
Telemark <- surrogate_data_1 %>%
  filter(trial_acronym == "Telemark Polyp I") %>%
  #  follow-up time until 1993
  filter(end_fu_date == "1993") %>%
  filter(str_detect(outcome_description, "Duke")) %>%
  #calculate all cancers diagnosed by arm and use as denominator (all stages + missing)
  mutate(denominator_screening = sum(as.numeric(numerator_screening)) +
           as.numeric(missing_screening)) %>%
  mutate(denominator_comparator = sum(as.numeric(numerator_comparator)) +
           as.numeric(missing_comparator)) %>%
  # late stage definition: Dukes' C or D
  filter(stage_category == "C") %>%
  mutate(outcome_description = "% Dukes' C CRC; follow-up for 10 years (1983-1993); Screening vs control") %>%
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


crc_surr3 <-
  rbind(
    NordICC,
    SCORE,
    FinnishRHS,
    Funen,
    Gothenburg,
    NORCCAP,
    Nottingham,
    PLCO,
    Minnesota_annual,
    Minnesota_biennial,
    Telemark
  )

crc_surr3 <- calculate_logrel(crc_surr3)

finaldata_primary_3 <- final_table(mortality_data, crc_surr3)

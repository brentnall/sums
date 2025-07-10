### BOWEL - Sensitivity (Stage III-IV) Estimates
setwd("../../../code")
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

crc_surr1 <- rbind(SCORE, PLCO) %>%
  calculate_logRR()

finaldata_primary_1 <- final_table(mortality_data, crc_surr1)


# Surrogate 3. -------------------------------------------------------------
surrogate_data_3 <-
  filter_surrogate_data(surrogate_data,
                        "Bowel",
                        "3. % of target cancer diagnosed at late stage")

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

crc_surr3 <- rbind(SCORE, PLCO)

crc_surr3 <- calculate_logrel(crc_surr3)

finaldata_primary_3 <- final_table(mortality_data, crc_surr3)

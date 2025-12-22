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


# Surrogate 4. -------------------------------------------------------------
surrogate_data_4 <- filter_surrogate_data(surrogate_data, "Bowel","1./4. Absolute incidence of late / early stage cancer")

surrogate_data_4 <- surrogate_data_4 %>%
  #SCORE - only late stage data available
  filter(trial_acronym!="SCORE") %>%
  # UKFSST - only incidence of CRC available
  filter(trial_acronym!="UKFSST")

#NordICC
NordICC <- surrogate_data_4 %>%
  filter(trial_acronym == "NordICC") %>%
  #early stage definition: Dukes’ stage A or B (reported)
  filter(stage_category == "Dukes' A or B") %>%
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
FinnishRHS <- surrogate_data_4 %>%
  filter(trial_acronym == "Finnish RHS CRC") %>%
  #whole trial estimate
  filter(subgroup == "Whole trial") %>%
  #  follow-up 31/12/2011
  filter(end_fu_date == "40908") %>%
  # early stage definition: N=0
  filter(stage_category == "N0") %>%
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
Funen <- surrogate_data_4 %>%
  filter(trial_acronym == "Funen") %>%
  #  follow-up August 1995
  filter(end_fu_date == "34912") %>%
  # rate ratio > relative risk
  filter(type_indicator == "Rate ratio") %>%
  #early stage definition: A+B (reported advanced=stage C, distant spread, and no classification )
  filter(stage_category == "A" | stage_category == "B") %>%
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
  mutate(trial_acronym = "Funen", stage_category = "Dukes' A or B")

#Gothenburg CRC
Gothenburg <- surrogate_data_4 %>%
  filter(trial_acronym == "Gothenburg CRC") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "Gothenburg CRC")]) %>%
  filter(period_after_rnd == mortality_data$period_after_rnd[mortality_data$trial_acronym ==
                                                               "Gothenburg CRC"]) %>%
  #rate ratio > relative risk
  filter(type_indicator == "Rate ratio") %>%
  # early stage definition: Duke's A-B
  filter(stage_category == "A" | stage_category == "B") %>%
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
  mutate(trial_acronym = "Gothenburg CRC", stage_category = "Dukes' A or B")

#Minnesota
Minnesota <- surrogate_data_4 %>%
  filter(trial_acronym == "Minnesota") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "Minnesota")]) %>%
  filter(period_after_rnd == mortality_data$period_after_rnd[mortality_data$trial_acronym ==
                                                               "Minnesota"]) %>%
  # early stage definition: Duke's A (reported?)
  filter(stage_category == "A") %>%
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
NORCCAP <- surrogate_data_4 %>%
  filter(trial_acronym == "NORCCAP") %>%
  #follow-up 31/12/2006
  filter(end_fu_date == "39082") %>%
  # early stage definition: Dukes' A or B (reported)
  filter(stage_category == "Dukes' A or B") %>%
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
Nottingham <- surrogate_data_4 %>%
  filter(trial_acronym == "Nottingham") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "Nottingham")]) %>%
  filter(period_after_rnd == mortality_data$period_after_rnd[mortality_data$trial_acronym ==
                                                               "Nottingham"]) %>%
  #rate ratio > relative risk
  filter(type_indicator == "Rate ratio") %>%
  # early stage definition: Dukes' A or B (C or D reported as late stage)
  filter(stage_category == "A" | stage_category == "B") %>%
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
  mutate(trial_acronym = "Nottingham", stage_category = "Dukes' A or B")

#	Telemark Polyp I
Telemark <- surrogate_data_4 %>%
  filter(trial_acronym == "Telemark Polyp I") %>%
  #  follow-up time until 1993
  filter(end_fu_date == "1993") %>%
  # early stage definition: Dukes' A or B
  filter(stage_category == "A" | stage_category == "B") %>%
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
  mutate(trial_acronym = "Telemark Polyp I", stage_category = "Dukes' A or B")

#PLCO (CRC)
PLCO <- surrogate_data_4 %>%
  filter(trial_acronym == "PLCO (CRC)") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "PLCO (CRC)")]) %>%
  # same definition of colorectal cancer as mortality (distal + proximal)
  filter(!str_detect(outcome_description, "Distal")) %>%
  filter(!str_detect(outcome_description, "Proximal")) %>%
  # early stage definition: stage I/II
  filter(stage_category == "I" | stage_category == "II") %>%
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
  mutate(trial_acronym = "PLCO (CRC)", stage_category = "Stage I-II")

crc_surr4 <-
  rbind(NordICC,
        FinnishRHS,
        Funen,
        Gothenburg,
        Minnesota,
        NORCCAP,
        Nottingham,
        Telemark,
        PLCO)

crc_surr4_rate <- crc_surr4 %>%
  filter(type_indicator == "Rate ratio") %>%
  calculate_logRR()

crc_surr4_risk <- crc_surr4 %>%
  filter(type_indicator == "Relative risk") %>%
  calculate_logrel()

crc_surr4 <- rbind(crc_surr4_rate, crc_surr4_risk)

finaldata_primary_4 <- final_table(mortality_data, crc_surr4)


# # Surrogate 5. ------------------------------------------------------------
surrogate_data_5 <-
  filter_surrogate_data(surrogate_data,
                        "Bowel",
                        "5. % of target cancers that are screen-detected")

#NordICC
NordICC <- surrogate_data_5 %>%
  filter(trial_acronym == "NordICC") 

#SCORE
SCORE <- surrogate_data_5 %>%
  filter(trial_acronym == "SCORE") %>%
  # same definition of colorectal cancer as mortality (distal + proximal)
  filter(str_detect(outcome_description, "distal")) %>%
  filter(str_detect(outcome_description, "proximal")) 

#UKFSST
UKFSST <- surrogate_data_5 %>%
  filter(trial_acronym == "UKFSST") %>%
  #whole trial estimate
  filter(subgroup == "Whole trial") %>%
  # same definition of colorectal cancer as mortality (distal + proximal)
  filter(!str_detect(outcome_description, "distal")) %>%
  filter(!str_detect(outcome_description, "proximal")) 

#Finnish RHS CRC
FinnishRHS <- surrogate_data_5 %>%
  filter(trial_acronym == "Finnish RHS CRC") %>%
  #whole trial estimate
  filter(subgroup == "Whole trial") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "Finnish RHS CRC")]) 


#Funen
Funen <- surrogate_data_5 %>%
  filter(trial_acronym == "Funen") %>%
  #  follow-up time = mortality outcome follow-up (August 1995)
  filter(end_fu_date == "34912") 

#Gothenburg CRC
Gothenburg <- surrogate_data_5 %>%
  filter(trial_acronym == "Gothenburg CRC") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "Gothenburg CRC")]) %>%
  filter(period_after_rnd == mortality_data$period_after_rnd[mortality_data$trial_acronym ==
                                                               "Gothenburg CRC"])

#Minnesota
Minnesota <- surrogate_data_5 %>%
  filter(trial_acronym == "Minnesota")

# separate annual/biennial in two rows
screening1 <- Minnesota %>%
  mutate(comparator_other_screening="Not applicable")


screening2 <- Minnesota %>%
  mutate(numerator_screening = numerator_comparator,
         denominator_screening = denominator_comparator) %>%
  mutate(screening_test = "FOBT (biennial)") %>%
  mutate(comparator_other_screening="Not applicable")

Minnesota <- rbind(screening1, screening2)

#NORCCAP
NORCCAP <- surrogate_data_5 %>%
  filter(trial_acronym == "NORCCAP") %>%
  #  follow-up time - 2006
  filter(end_fu_date == "39082") %>%
  filter(screening_test == "Flexible sigmoidoscopy (+/- FOBT)")

#Nottingham
Nottingham <- surrogate_data_5 %>%
  filter(trial_acronym == "Nottingham") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "Nottingham")]) %>%
  filter(period_after_rnd == mortality_data$period_after_rnd[mortality_data$trial_acronym ==
                                                               "Nottingham"]) 

#	Telemark Polyp I
Telemark <- surrogate_data_5 %>%
  filter(trial_acronym == "Telemark Polyp I") 

#PLCO (CRC)
PLCO <- surrogate_data_5 %>%
  filter(trial_acronym == "PLCO (CRC)") %>%
  # total (baseline + screening rounds at 3y/5y)
  filter(str_detect(outcome_description, "total")) 

crc_surr5 <-
  rbind(
    NordICC,
    SCORE,
    UKFSST,
    FinnishRHS,
    Funen,
    Gothenburg,
    Minnesota,
    NORCCAP,
    Nottingham,
    Telemark,
    PLCO
  )

crc_surr5 <- crc_surr5 %>%
  select(
    trial_acronym,
    screening_test,
    comparator_other_screening,
    outcome_description,
    numerator_screening,
    denominator_screening
  ) %>%
  calculate_prop()

finaldata_primary_5 <- final_table(mortality_data, crc_surr5)


# Surrogate 6. ------------------------------------------------------------
surrogate_data_6 <- filter_surrogate_data(surrogate_data, "Bowel","6. % of high-grade target cancers that are screen-detected")
N <- min_N_trials(surrogate_data_6)
N ### 1 trial

crc_surr6 <- surrogate_data_6 %>%
  #  all screening rounds
  filter(`Effective number of screening rounds`==mortality_data$Effective.number.of.screening.rounds[(mortality_data$trial_acronym=="Nottingham")]) %>%
  select(trial_acronym, screening_test, outcome_description, numerator_screening, denominator_screening) %>%
  calculate_prop()

finaldata_primary_6 <- final_table(mortality_data, crc_surr6)


# Surrogate 7. ------------------------------------------------------------
surrogate_data_7 <-
  filter_surrogate_data(surrogate_data, "Bowel", "7. Diagnostic yield at screening")

#NordICC
NordICC <- surrogate_data_7 %>%
  filter(trial_acronym == "NordICC") %>%
  #adenoma + CRC
  summarise(
    screening_test = first(screening_test),
    comparator_other_screening = first(comparator_other_screening),
    numerator_screening = sum(as.numeric(numerator_screening)),
    denominator_screening = mean(as.numeric(denominator_screening))
  ) %>%
  mutate(trial_acronym = "NordICC", outcome_description = "Diagnostic yield at screening - CRC + Adenomas (per 100 analysed); once-only colonoscopy; Screening arm only")

#SCORE
SCORE <- surrogate_data_7 %>%
  filter(trial_acronym == "SCORE") %>%
  # same definition of colorectal cancer as mortality (distal + proximal)
  filter(str_detect(outcome_description, "distal")) %>%
  filter(str_detect(outcome_description, "proximal")) %>%
  select(
    trial_acronym,
    screening_test,
    comparator_other_screening,
    outcome_description,
    numerator_screening,
    denominator_screening
  )

#UKFSST
UKFSST <- surrogate_data_7 %>%
  filter(trial_acronym == "UKFSST") %>%
  #whole trial estimate
  filter(subgroup == "Whole trial") %>%
  # same definition of colorectal cancer as mortality (distal + proximal)
  filter(!str_detect(outcome_description, "distal")) %>%
  filter(!str_detect(outcome_description, "proximal")) %>%
  # high risk pops + CRC
  summarise(
    screening_test = first(screening_test),
    comparator_other_screening = first(comparator_other_screening),
    numerator_screening = sum(as.numeric(numerator_screening)),
    denominator_screening = mean(as.numeric(denominator_screening))
  ) %>%
  mutate(trial_acronym = "UKFSST", outcome_description = "Diagnostic yield at screening (high-risk polyps + CRC); per 100 people analysed; one round of flexible sigmoidoscopy")

#Finnish RHS CRC
FinnishRHS <- surrogate_data_7 %>%
  filter(trial_acronym == "Finnish RHS CRC") %>%
  #whole trial estimate
  filter(subgroup == "Whole trial") %>%
  # include all screening rounds (31/12/2012)
  filter(str_detect(numerator_definition, "2012")) %>%
  select(
    trial_acronym,
    screening_test,
    comparator_other_screening,
    outcome_description,
    numerator_screening,
    denominator_screening
  )
  

#Funen
Funen <- surrogate_data_7 %>%
  filter(trial_acronym == "Funen") %>%
  #  follow-up time = mortality outcome follow-up (August 1995)
  filter(end_fu_date == "34912") %>%
  #adenoma + CRC
  summarise(
    screening_test = first(screening_test),
    comparator_other_screening = first(comparator_other_screening),
    numerator_screening = sum(as.numeric(numerator_screening)),
    denominator_screening = mean(as.numeric(denominator_screening))
  ) %>%
  mutate(trial_acronym = "Funen", outcome_description = "Diagnostic yield at screening - CRC + Adenoma >=10 mm (per 1,000 people randomised);  5 screening rounds (August 1985 - August 1995)")

#Gothenburg CRC
Gothenburg <- surrogate_data_7 %>%
  filter(trial_acronym == "Gothenburg CRC") %>%
  #  number of screening rounds = mortality outcome
  filter(
    `Effective number of screening rounds` == mortality_data$Effective.number.of.screening.rounds[(mortality_data$trial_acronym ==
                                                                                                     "Gothenburg CRC")]
  ) %>%
  # all adenomas
  filter(!str_detect(numerator_definition, ">=1")) %>%
  summarise(
    screening_test = first(screening_test),
    comparator_other_screening = first(comparator_other_screening),
    numerator_screening = sum(as.numeric(numerator_screening)),
    denominator_screening = mean(as.numeric(denominator_screening))
  ) %>%
  mutate(trial_acronym = "Gothenburg CRC", outcome_description = "Diagnostic yield at screening, CRC + adenomas (2-3 rounds); per 1,000 people randomised; all cohorts (DOB 1918-1931); up to the end of 2nd rescreening")


#Minnesota
Minnesota <- surrogate_data_7 %>%
  filter(trial_acronym == "Minnesota")

# separate annual/biennial in two rows
screening1 <- Minnesota %>%
  select(
    trial_acronym,
    screening_test,
    comparator_other_screening,
    outcome_description,
    numerator_screening,
    denominator_screening
  ) %>%
  mutate(comparator_other_screening="Not applicable")

screening2 <- Minnesota %>%
  mutate(numerator_screening = numerator_comparator,
         denominator_screening = denominator_comparator) %>%
  mutate(screening_test = "FOBT (biennial)") %>%
  select(
    trial_acronym,
    screening_test,
    comparator_other_screening,
    outcome_description,
    numerator_screening,
    denominator_screening
  )%>%
  mutate(comparator_other_screening="Not applicable")


Minnesota <- rbind(screening1, screening2)

#NORCCAP
NORCCAP <- surrogate_data_7 %>%
  filter(trial_acronym == "NORCCAP") %>%
  #whole trial
  filter(subgroup == "Whole trial") %>%
  # comaprator group = control (as mortality)
  filter(comparator_other_screening == "Not applicable") %>%
  #all adenoma + CRC
  filter(!str_detect(outcome_description, "advanced adenoma")) %>%
  summarise(
    screening_test = first(screening_test),
    comparator_other_screening = first(comparator_other_screening),
    numerator_screening = sum(as.numeric(numerator_screening)),
    denominator_screening = mean(as.numeric(denominator_screening))
  ) %>%
  mutate(trial_acronym = "NORCCAP", outcome_description = "Diagnostic yield at screening (CRC + any adenoma); per 1,000 people analysed; 1 screening round; FS +/- FOBT")


#Nottingham
Nottingham <- surrogate_data_7 %>%
  filter(trial_acronym == "Nottingham") %>%
  #3-6 screening rounds
  filter(`Effective number of screening rounds` == "3-6") %>%
  # all adenomas + CRC
  filter(
    !str_detect(outcome_description, "cm") 
  ) %>%
  summarise(
    screening_test = first(screening_test),
    comparator_other_screening = first(comparator_other_screening),
    numerator_screening = sum(as.numeric(numerator_screening)),
    denominator_screening = mean(as.numeric(denominator_screening))
  ) %>%
  mutate(trial_acronym = "Nottingham", outcome_description = "Diagnostic yield at screening, CRC + Adenomas (per 1,000 people analysed; 3-6 screening rounds offered; up to 30/06/1995 FOBT only")

#	Telemark Polyp I
Telemark <- surrogate_data_7 %>%
  filter(trial_acronym == "Telemark Polyp I") %>%
  #  number of screening rounds = mortality outcome (1-2 rounds)
  filter(
    `Effective number of screening rounds` == mortality_data$Effective.number.of.screening.rounds[(mortality_data$trial_acronym ==
                                                                                                     "Telemark Polyp I")]
  ) %>%
  filter(str_detect(outcome_description, "1983 and 1985")) %>%
  filter(author_year == "Hoff 1996") %>%
  filter(
      str_detect(outcome_description, "total neoplasias")
  ) %>%
  select(
    trial_acronym,
    screening_test,
    comparator_other_screening,
    outcome_description,
    numerator_screening,
    denominator_screening
  )
 
#PLCO (CRC)
PLCO <- surrogate_data_7 %>%
  filter(trial_acronym == "PLCO (CRC)") %>%
  # total (baseline + screening rounds at 3y/5y)
  filter(str_detect(outcome_description, "total")) %>%
  select(
    trial_acronym,
    screening_test,
    comparator_other_screening,
    outcome_description,
    numerator_screening,
    denominator_screening
  )


crc_surr7 <-
  rbind(
    NordICC,
    SCORE,
    UKFSST,
    FinnishRHS,
    Funen,
    Gothenburg,
    Minnesota,
    NORCCAP,
    Nottingham,
    Telemark,
    PLCO
  )

crc_surr7 <- calculate_prop(crc_surr7)

finaldata_primary_7 <- final_table(mortality_data, crc_surr7)

finaldata_primary_7 <- finaldata_primary_7

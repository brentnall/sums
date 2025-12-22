### CERVICAL 
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
  filter_mortality_data(mortality_data, cancer = "Cervical")

mortality_data_1 <- mortality_data %>%
  #rate ratio > relative risk
  filter(trial_acronym != "Tamil Nadu/Dindigul") %>%
  filter(type_indicator == "Rate ratio")

#Tamil Nadu/Dindigul only reports HR or relative risk - calculated rate ratio from crude numbers
mortality_data_2 <- mortality_data %>%
  #rate ratio > relative risk
  filter(trial_acronym == "Tamil Nadu/Dindigul") %>%
  filter(type_indicator == "Hazard ratio")

mortality_data <- rbind(mortality_data_1, mortality_data_2)
mortality_data <- mortality_data %>%
  calculate_logRR() %>%
  replace_with_reported()

# Surrogate 4. -------------------------------------------------------------
surrogate_data_4 <-
  filter_surrogate_data(surrogate_data,
                        "Cervical",
                        "1./4. Absolute incidence of late / early stage cancer")

# Finnish RHS Cervix
Finnish_RHS_Cervix <- surrogate_data_4 %>%
  filter(trial_acronym == "Finnish RHS Cervix") %>%
  #rate ratio > relative risk
  filter(type_indicator == "Rate ratio") %>%
  #early stage definition: CIN3/AIS (reported)
  filter(stage_category == "Non-invasive (CIN3/AIS)") %>%
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
  )

# Mumbai Breast Cervix
Mumbai_Breast_Cervix <- surrogate_data_4 %>%
  filter(trial_acronym == "Mumbai Breast Cervix") %>%
  #rate ratio > relative risk
  filter(type_indicator == "Rate ratio") %>%
  #late stage definition: Stage <IIB (reported)
  filter(stage_category == "Stage <IIB") %>% 
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
  )

# Osmanabad Cervix
Osmanabad_Cervix <- surrogate_data_4 %>%
  filter(trial_acronym == "Osmanabad Cervix") %>%
  # late stage definition: Stage I
  filter(stage_category == "Stage I") %>%
  # rate ratio > relative risk (and adjusted hazard ratio)
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
  )


# Tamil Nadu/Dindigul
Tamil_Nadu_Dindigul <- surrogate_data_4 %>%
  filter(trial_acronym == "Tamil Nadu/Dindigul") %>%
  # whole trial estimates
  filter(subgroup == "Whole trial") %>%
  # follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "Tamil Nadu/Dindigul")]) %>%
  # late stage definition: 	Preclinical = Stage I (reported)
  filter(stage_category == "Stage I") %>%
  # rate ratio > relative risk (and adjusted hazard ratio)
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
  )

cervical_surr4 <-
  rbind(Finnish_RHS_Cervix,
        Mumbai_Breast_Cervix,
        Osmanabad_Cervix,
        Tamil_Nadu_Dindigul)

cervical_surr4 <- calculate_logRR(cervical_surr4)

finaldata_primary_4 <- final_table(mortality_data, cervical_surr4)


# Surrogate 5. ------------------------------------------------------------
surrogate_data_5 <-
  filter_surrogate_data(surrogate_data,
                        "Cervical",
                        "5. % of target cancers that are screen-detected")

# Osmanabad Cervix
Osmanabad_Cervix <- surrogate_data_5 %>%
  filter(trial_acronym == "Osmanabad Cervix") %>%
  #FU 31/12/2007
  filter(end_fu_date == "39447") 

# Tamil Nadu/Dindigul
Tamil_Nadu_Dindigul <- surrogate_data_5 %>%
  filter(trial_acronym == "Tamil Nadu/Dindigul") %>%
  # follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "Tamil Nadu/Dindigul")]) 


cervical_surr5 <- rbind(Osmanabad_Cervix, Tamil_Nadu_Dindigul)

cervical_surr5 <- cervical_surr5 %>%
  select(
    trial_acronym,
    screening_test,
    comparator_other_screening,
    outcome_description,
    numerator_screening,
    denominator_screening
  ) %>%
  calculate_prop()

finaldata_primary_5 <- final_table(mortality_data, cervical_surr5)


# Surrogate 6. ------------------------------------------------------------
surrogate_data_6 <- filter_surrogate_data(surrogate_data, "Cervical","6. % of high-grade target cancers that are screen-detected")
N <- min_N_trials(surrogate_data_6)
N ### 0 trials


# Surrogate 7. ------------------------------------------------------------
surrogate_data_7 <-
  filter_surrogate_data(surrogate_data, "Cervical", "7. Diagnostic yield at screening")

# Osmanabad Cervix
Osmanabad_Cervix <- surrogate_data_7 %>%
  filter(trial_acronym == "Osmanabad Cervix") %>%
  # first publication
  filter(end_fu_date == "November 2003") %>%
  filter(
    str_detect(outcome_description, "invasive cancer") &
      str_detect(outcome_description, "CIN2/3")
  ) %>%
  select(
    trial_acronym,
    screening_test,
    comparator_other_screening,
    outcome_description,
    numerator_screening,
    denominator_screening
  )


# Tamil Nadu/Dindigul
Tamil_Nadu_Dindigul <- surrogate_data_7 %>%
  filter(trial_acronym == "Tamil Nadu/Dindigul") %>%
  # whole trial estimates
  filter(subgroup == "Whole trial") %>%
  # first publication
  filter(author_year == "Sankar 2004") %>%
  # CIN2/3 + invasive cancers
  filter(
    str_detect(outcome_description, "invasive") &
      str_detect(outcome_description, "CIN2/3")
  ) %>%
  select(
    trial_acronym,
    screening_test,
    comparator_other_screening,
    outcome_description,
    numerator_screening,
    denominator_screening
  )

cervical_surr7 <- rbind(Osmanabad_Cervix, Tamil_Nadu_Dindigul)

cervical_surr7 <- calculate_prop(cervical_surr7)

finaldata_primary_7 <- final_table(mortality_data, cervical_surr7)

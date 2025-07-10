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



# Surrogate 1. ------------------------------------------------------
surrogate_data_1 <-
  filter_surrogate_data(surrogate_data,
                        "Cervical",
                        "1./4. Absolute incidence of late / early stage cancer")


# Finnish RHS Cervix
Finnish_RHS_Cervix <- surrogate_data_1 %>%
  filter(trial_acronym == "Finnish RHS Cervix") %>%
  #rate ratio > relative risk
  filter(type_indicator == "Rate ratio") %>%
  #late stage definition: invasive
  filter(stage_category == "Invasive") %>%
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
Mumbai_Breast_Cervix <- surrogate_data_1 %>%
  filter(trial_acronym == "Mumbai Breast Cervix") %>%
  #rate ratio > relative risk
  filter(type_indicator == "Rate ratio") %>%
  #late stage definition: Stage b	%IIB (reported)
  filter(stage_category == "Stage >=IIB") %>%
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
Osmanabad_Cervix <- surrogate_data_1 %>%
  filter(trial_acronym == "Osmanabad Cervix") %>%
  # late stage definition: Stage II+ (reported)
  filter(stage_category == "Stage II+") %>%
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
Tamil_Nadu_Dindigul <- surrogate_data_1 %>%
  filter(trial_acronym == "Tamil Nadu/Dindigul") %>%
  # whole trial estimates
  filter(subgroup == "Whole trial") %>%
  # follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "Tamil Nadu/Dindigul")]) %>%
  # late stage definition: Stage II+ (reported)
  filter(stage_category == "Stage II+") %>%
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

cervical_surr1 <-
  rbind(Finnish_RHS_Cervix,
        Mumbai_Breast_Cervix,
        Osmanabad_Cervix,
        Tamil_Nadu_Dindigul)

cervical_surr1 <- calculate_logRR(cervical_surr1)

finaldata_primary_1 <- final_table(mortality_data, cervical_surr1)


# Surrogate 3. -------------------------------------------------------------
surrogate_data_3 <-
  filter_surrogate_data(surrogate_data,
                        "Cervical",
                        "3. % of target cancer diagnosed at late stage")

# Mumbai Breast Cervix
Mumbai_Breast_Cervix <- surrogate_data_3 %>%
  filter(trial_acronym == "Mumbai Breast Cervix") %>%
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
    comparator_noscreening,
  )

# Osmanabad Cervix
Osmanabad_Cervix <- surrogate_data_3 %>%
  filter(trial_acronym == "Osmanabad Cervix") %>%
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
    comparator_noscreening,
  )

# Tamil Nadu/Dindigul
Tamil_Nadu_Dindigul <- surrogate_data_3 %>%
  filter(trial_acronym == "Tamil Nadu/Dindigul") %>%
  # follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "Tamil Nadu/Dindigul")]) %>%
  filter(period_after_rnd == mortality_data$period_after_rnd[mortality_data$trial_acronym ==
                                                               "Tamil Nadu/Dindigul"]) %>%
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
    comparator_noscreening,
  )

cervical_surr3 <-
  rbind(Mumbai_Breast_Cervix,
        Osmanabad_Cervix,
        Tamil_Nadu_Dindigul)

cervical_surr3 <- calculate_logrel(cervical_surr3)

finaldata_primary_3 <- final_table(mortality_data, cervical_surr3)

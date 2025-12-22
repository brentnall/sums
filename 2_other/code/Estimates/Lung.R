### LUNG 
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
  filter_mortality_data(mortality_data, cancer = "Lung") %>%
  filter(trial_acronym != "Mayo Lung Project") %>%
  #rate ratio > relative risk/HR, for the studies that report it
  filter(type_indicator == "Rate ratio") %>%
  calculate_logRR() %>%
  replace_with_reported()

#Czech study
mortality_Czech <- mortality_data %>%
  filter(str_detect(trial_acronym, "Czech study")) %>%
  filter(mortality_endpoint == "Cancer-specific") %>%
  filter(period_after_rnd == "15 years") %>%
  filter(!str_detect(outcome_description, "diagnosed within")) %>%
  calculate_logrel() %>%
  replace_with_reported()

#LUSI
mortality_LUSI <- mortality_data %>%
  filter(str_detect(trial_acronym, "LUSI")) %>%
  filter(mortality_endpoint == "Cancer-specific") %>%
  filter(subgroup == "Whole trial") %>%
  filter(str_detect(outcome_description, "until 30 April 2018")) %>%
  filter(type_indicator == "Relative risk") %>%
  calculate_logrel() %>%
  replace_with_reported()

#Mayo Lung Project
mortality_Mayo <- mortality_data %>%
  filter(str_detect(trial_acronym, "Mayo Lung Project")) %>%
  filter(mortality_endpoint == "Cancer-specific") %>%
  #FU until 1983
  filter(end_fu_date == "30498") %>%
  filter(type_indicator == "Relative risk") %>%
  filter(!str_detect(outcome_description, "death certificate")) %>%
  calculate_logrel() %>%
  replace_with_reported()


mortality_data <-
  rbind(mortality_data_1,
        mortality_Czech,
        mortality_LUSI,
        mortality_Mayo)
  
# Surrogate 4. -------------------------------------------------------------
surrogate_data_4 <-
  filter_surrogate_data(surrogate_data,
                        "Lung",
                        "1./4. Absolute incidence of late / early stage cancer")

#Czech study - only estimate available screening ongoing
czech_study <- surrogate_data_4 %>%
  filter(trial_acronym == "Czech study") %>%
  filter(stage_category=="Stage I") %>%
  select(
    trial_acronym,
    screening_test,
    outcome_description,
    numerator_screening,
    numerator_comparator,
    denominator_screening,
    type_indicator,
    denominator_comparator,
    reported_definition_late,
    reported_definition_early,
    stage_category
  )
  


#DLCST
DLCST <- surrogate_data_4 %>%
  filter(trial_acronym == "DLCST") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "DLCST")]) %>%
  # early stage definition: Stage I + II (reported)
  filter(stage_category == "Stage I + II") %>%
  select(
    trial_acronym,
    screening_test,
    outcome_description,
    numerator_screening,
    numerator_comparator,
    denominator_screening,
    type_indicator,
    denominator_comparator,
    reported_definition_late,
    reported_definition_early,
    stage_category
  )


#	Johns Hopkins  
Johns_Hopkins <- surrogate_data_4 %>%
  filter(trial_acronym == "Johns Hopkins") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "Johns Hopkins")]) %>%
  filter(period_after_rnd == mortality_data$period_after_rnd[mortality_data$trial_acronym ==
                                                               "Johns Hopkins"]) %>%
  # early stage definition: stage 0-1 (not reported as 'late stage', but estimate available)
  filter(stage_category == "Stage 0" |
           stage_category == "Stage 1") %>%
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
  mutate(trial_acronym = "Johns Hopkins", stage_category = "Stage 0-1")

#PLCO
PLCO <- surrogate_data_4 %>%
  filter(trial_acronym == "PLCO (Lung)") %>%
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
  mutate(outcome_description = "% non-small cell lung cancers diagnosed at stage I-II; 13 years of follow-up or December 31st, 2009") %>%
  mutate(trial_acronym = "PLCO (Lung)", stage_category = "Stage I-II (non-small)")

#UKLS
UKLS <- surrogate_data_4 %>%
  filter(trial_acronym == "UKLS") %>%
  # rate ratio > relative risk
  filter(rate_unit == "Per 1,000 PY") %>%
  #correct typo
  mutate(type_indicator = "Rate ratio") %>%
  filter(stage_category == "Stage I" |
           stage_category == "Stage II") %>%
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
  mutate(trial_acronym = "UKLS", stage_category = "Stage I-II")


#NELSON 
NELSON <- surrogate_data_4 %>%
  filter(trial_acronym == "NELSON") %>%
  # early stage definition: Stage I-II (reported late stage = stage III+)
  filter(stage_category == "Stage I" |
           stage_category == "Stage II") %>%
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
  mutate(trial_acronym = "NELSON", stage_category = "Stage I-II")

#Lung Screening Study (LSS)
LSS <- surrogate_data_4 %>%
  filter(trial_acronym == "Lung Screening Study (LSS)") %>%
  # definition: baseline, 1 eound of screening and interval cancers
  filter(`Effective number of screening rounds` == "2") %>%
  #early stage definition: stage I-II (stage III-IV reported as late)
  filter(stage_category == "Stage I" |
           stage_category == "Stage II") %>%
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
  mutate(trial_acronym = "Lung Screening Study (LSS)", stage_category =
           "Stage I-II")

#Memorial Sloan-Kettering
MSK <- surrogate_data_4 %>%
  filter(trial_acronym == "Memorial Sloan-Kettering") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "Memorial Sloan-Kettering")]) %>%
  filter(period_after_rnd == mortality_data$period_after_rnd[mortality_data$trial_acronym ==
                                                               "Memorial Sloan-Kettering"]) %>%
  #early stage definition: stage 0-1 (stage 2-4 only late stage available for surrogate 3)
  filter(stage_category == "Stage 0" |
           stage_category == "Stage 1") %>%
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
  mutate(trial_acronym = "Memorial Sloan-Kettering", stage_category = "Stage 0-1")

#DANTE
DANTE <- surrogate_data_4 %>%
  filter(trial_acronym == "DANTE") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "DANTE")]) %>%
  filter(period_after_rnd == mortality_data$period_after_rnd[mortality_data$trial_acronym ==
                                                               "DANTE"]) %>%
  # early stage definition: stage I (reported)
  filter(stage_category == "Stage I") %>%
  select(
    trial_acronym,
    screening_test,
    outcome_description,
    numerator_screening,
    numerator_comparator,
    denominator_screening,
    type_indicator,
    denominator_comparator,
    reported_definition_late,
    reported_definition_early,
    stage_category
  )

#LUSI
LUSI <- surrogate_data_4 %>%
  filter(trial_acronym == "LUSI") %>%
  filter(subgroup == "Whole trial") %>%
  # 7 years FU - halfway point between end of screening and mortality
  filter(period_after_rnd == "7 years") %>%
  # late stage definition: stage I (reported)
  filter(stage_category == "Stage I") %>%
  select(
    trial_acronym,
    screening_test,
    outcome_description,
    numerator_screening,
    numerator_comparator,
    denominator_screening,
    type_indicator,
    denominator_comparator,
    reported_definition_late,
    reported_definition_early,
    stage_category
  )


#ITALUNG
ITALUNG <- surrogate_data_4 %>%
  filter(trial_acronym == "ITALUNG") %>%
  # follow-up time = mortality outcome follow-up not available, closest time point (prior to mortality) is 31/12/2013
  filter(end_fu_date == "41639") %>%
  #early stage defintion: Stage I-II (late stage III_IV reported)
  filter(stage_category == "Stage I" |
           stage_category == "Stage II") %>%
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
  mutate(trial_acronym = "ITALUNG", stage_category = "Stage I-II")

#	Mayo Lung Project
Mayo <- surrogate_data_4 %>%
  filter(trial_acronym == "Mayo Lung Project") %>%
  # follow-up time = mortality outcome follow-up not available, closest time point is 01/07/1983
  filter(end_fu_date == "30498") %>%
  #early stage definition: Stage 0-II(reported)
  filter(stage_category == "Stage 0-II") %>%
  select(
    trial_acronym,
    screening_test,
    outcome_description,
    numerator_screening,
    numerator_comparator,
    denominator_screening,
    type_indicator,
    denominator_comparator,
    reported_definition_late,
    reported_definition_early,
    stage_category
  )

# MILD - only subgroup
MILD <- surrogate_data_4 %>%
  filter(trial_acronym == "MILD") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "MILD")]) %>%
  filter(period_after_rnd == mortality_data$period_after_rnd[mortality_data$trial_acronym ==
                                                               "MILD"]) %>%
  # same comparator as mortality
  filter(comparator_noscreening == "Usual care") %>%
  # rate ratio > relative risk
  filter(type_indicator == "Rate ratio") %>%
  #early stage definition: stage I
  filter(stage_category == "Stage I") %>%
  select(
    trial_acronym,
    screening_test,
    outcome_description,
    numerator_screening,
    numerator_comparator,
    denominator_screening,
    type_indicator,
    denominator_comparator,
    reported_definition_late,
    reported_definition_early,
    stage_category
  )


#NLST
NLST <- surrogate_data_4 %>%
  filter(trial_acronym == "NLST") %>%
  filter(subgroup == "Whole trial") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "NLST")]) %>%
  filter(period_after_rnd == mortality_data$period_after_rnd[mortality_data$trial_acronym ==
                                                               "NLST"]) %>%
  #stage definition: stage III-IV (TBD)
  filter(stage_category == "Stage I" |
           stage_category == "Stage II") %>%
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
  mutate(trial_acronym = "NLST", stage_category = "Stage I-II")


##!add czech study when estimate is decided
lung_surr4 <-
  rbind(czech_study,
    DANTE,
        DLCST,
        ITALUNG,
        Johns_Hopkins,
        LSS,
        LUSI,
        Mayo,
        MILD,
        MSK,
        NELSON,
        NLST,
        PLCO,
        UKLS)

lung_surr4_rate <- lung_surr4 %>%
  filter(type_indicator == "Rate ratio") %>%
  calculate_logRR()

lung_surr4_risk <- lung_surr4 %>%
  filter(type_indicator == "Relative risk") %>%
  calculate_logrel()

lung_surr4 <- rbind(lung_surr4_rate, lung_surr4_risk)

finaldata_primary_4 <- final_table(mortality_data, lung_surr4)

# Surrogate 5. ------------------------------------------------------------
surrogate_data_5 <- filter_surrogate_data(surrogate_data, "Lung","5. % of target cancers that are screen-detected")

#Czech study
czech_study <- surrogate_data_5 %>%
  filter(trial_acronym=="Czech study") %>%
  #screening period
  filter(str_detect(outcome_description, "6 year")) 

#DLCST 
# DLCST_denominator <- surrogate_data_4 %>%
#   filter(trial_acronym=="DLCST") %>%
#   filter(`Achieved follow-up time at analysis, Comparator`=="Median 9.80 years") %>%
#   filter(stage_category=="Stage I + II" | stage_category=="Stage III-IV" ) %>%
#   summarise(denominator_screening=sum(as.numeric(numerator_screening)))

# DLCST_denominator <- surrogate_data_4 %>%
#      filter(trial_acronym=="DLCST") %>%
#   filter(end_fu_date=="40268") %>%
#   filter(stage_category=="Stage I" |
#            stage_category=="Stage II" |
#            stage_category=="Stage III" |
#            stage_category=="Stage IV") %>%
#   summarise(denominator_screening=sum(as.numeric(numerator_screening)))
  

DLCST <- surrogate_data_5 %>%
    filter(trial_acronym=="DLCST") %>%
    filter(`Use for primary analysis`=="yes")

#Johns Hopkins
Johns_Hopkins <- surrogate_data_5 %>%
  filter(trial_acronym=="Johns Hopkins")

#PLCO
PLCO <- surrogate_data_5 %>%
  filter(trial_acronym=="PLCO (Lung)") %>%
    filter(`Use for primary analysis`=="yes")


#UKLS
UKLS <- surrogate_data_5 %>%
  filter(trial_acronym=="UKLS") 

#NELSON
NELSON <- surrogate_data_5 %>%
  filter(trial_acronym=="NELSON") 


#Lung Screening Study (LSS)
LSS <- surrogate_data_5 %>%
  filter(trial_acronym=="Lung Screening Study (LSS)") %>%
  filter(`Effective number of screening rounds`=="2") 

#Memorial Sloan-Kettering
MSK <- surrogate_data_5 %>%
  filter(trial_acronym=="Memorial Sloan-Kettering") %>%
  #  follow-up time = mortality outcome follow-up 
  filter(end_fu_date=="1982")

#DANTE  
DANTE <- surrogate_data_5 %>%
  filter(trial_acronym=="DANTE") %>%
  #  follow-up time = mortality outcome follow-up 
  filter(end_fu_date==mortality_data$end_fu_date[(mortality_data$trial_acronym=="DANTE")]) %>%
  filter(period_after_rnd==mortality_data$period_after_rnd[mortality_data$trial_acronym=="DANTE"]) 

#LUSI
LUSI <- surrogate_data_5 %>%
  filter(trial_acronym=="LUSI") %>%
#  filter(subgroup=="Whole trial") %>%
#  filter(str_detect(outcome_description, "screening period")) %>%
    filter(`Use for primary analysis`=="yes")


#ITALUNG
ITALUNG <-surrogate_data_5 %>%
  filter(trial_acronym=="ITALUNG") %>%
  # follow-up time = mortality outcome follow-up not available, closest time point (prior to mortality) is 31/12/2013
  filter(end_fu_date=="41639")

#	Mayo Lung Project
Mayo <- surrogate_data_5 %>%
  filter(trial_acronym=="Mayo Lung Project") %>%
  # follow-up time = mortality outcome follow-up not available, closest time point is 01/07/1983
  filter(end_fu_date=="30498") 

# MILD - only subgroup
MILD <- surrogate_data_5 %>%
    filter(trial_acronym=="MILD") %>%
    filter(`Use for primary analysis`=="yes")


#NLST
# NLST <- surrogate_data_5  %>%
#   filter(trial_acronym=="NLST")
# #calculate from surrogate 1/4
# NLST <- surrogate_data_4 %>%
#   filter(trial_acronym=="NLST") %>%
#   filter(end_fu_date=="1-year since T2") %>%
#   filter(stage_category=="Stage I" | stage_category=="Stage II" | stage_category=="Stage III" | stage_category=="Stage IV" ) %>%
#   summarise(trial_acronym=first(trial_acronym),
#             screening_test=first(screening_test),
#             outcome_description="am",
#             #7 missing in screening arm
#             denominator_screening=sum(as.numeric(numerator_screening))+7,
#             numerator_screening=123
#            )

lung_surr5 <- rbind(czech_study, DANTE, DLCST, ITALUNG, Johns_Hopkins, LSS, LUSI, Mayo, MILD, MSK, NELSON, PLCO, UKLS)

lung_surr5 <- lung_surr5 %>%
  select(
    trial_acronym,
    screening_test,
    comparator_other_screening,
    outcome_description,
    numerator_screening,
    denominator_screening
  ) %>%
  calculate_prop()

finaldata_primary_5 <- final_table(mortality_data, lung_surr5)
  
# Surrogate 6. ------------------------------------------------------------
surrogate_data_6 <- filter_surrogate_data(surrogate_data, "Lung","6. % of high-grade target cancers that are screen-detected")
N <- min_N_trials(surrogate_data_6)
N ### 0 trials


# Surrogate 7. ------------------------------------------------------------
surrogate_data_7 <-
  filter_surrogate_data(surrogate_data, "Lung", "7. Diagnostic yield at screening")

#Czech study
czech_study <- surrogate_data_7 %>%
  filter(trial_acronym == "Czech study") 

#DLCST
DLCST <- surrogate_data_7 %>%
  filter(trial_acronym == "DLCST") %>%
  #all screening rounds
  filter(str_detect(numerator_definition, "after 5 screening rounds")) 

#	Johns Hopkins
Johns_Hopkins <- surrogate_data_7 %>%
  filter(trial_acronym == "Johns Hopkins") %>%
  #follow-up time = mortality outcome not available; closest is SD cancers 5-7y since rnd
  filter(period_after_rnd == "5-7 years") 

#PLCO
PLCO <- surrogate_data_7 %>%
  filter(trial_acronym == "PLCO (Lung)") 

#UKLS
UKLS <- surrogate_data_7 %>%
  filter(trial_acronym == "UKLS") 

#NELSON
NELSON <- surrogate_data_7 %>%
  filter(trial_acronym == "NELSON") 

#Lung Screening Study (LSS)
LSS <- surrogate_data_7 %>%
  filter(trial_acronym == "Lung Screening Study (LSS)") %>%
  filter(str_detect(outcome_description, "Baseline plus Year 1")) 

#Memorial Sloan-Kettering
MSK <- surrogate_data_7 %>%
  filter(trial_acronym == "Memorial Sloan-Kettering") %>%
  #baseline + incident
  filter(str_detect(outcome_description, "baseline")) %>%
  filter(str_detect(outcome_description, "incidence")) 

#DANTE
DANTE <- surrogate_data_7 %>%
  filter(trial_acronym == "DANTE") %>% ###only baseline
  filter(author_year == "Infante 2008") 

#LUSI
LUSI <- surrogate_data_7 %>%
  filter(trial_acronym == "LUSI")

#ITALUNG
ITALUNG <- surrogate_data_7 %>%
  filter(trial_acronym == "ITALUNG") 

#	Mayo Lung Project
Mayo <- surrogate_data_7 %>%
  filter(trial_acronym == "Mayo Lung Project") %>%
  # total screening period
  filter(`Effective number of screening rounds` == "18") 

# MILD - only subgroup
MILD <- surrogate_data_7 %>%
  filter(trial_acronym == "MILD") %>%
  filter(author_year == "Infante 2017") %>%
  # same comparator as mortality (combined)
  filter(comparator_other_screening == "Not applicable") 

#NLST
NLST <- surrogate_data_7 %>%
  filter(trial_acronym == "NLST") 

lung_surr7 <-
  rbind(
    czech_study,
    DANTE,
    DLCST,
    ITALUNG,
    Johns_Hopkins,
    LSS,
    LUSI,
    Mayo,
    MILD,
    MSK,
    NELSON,
    NLST,
    PLCO,
    UKLS
  )

lung_surr7 <- lung_surr7 %>%
  select(
    trial_acronym,
    screening_test,
    comparator_other_screening,
    outcome_description,
    numerator_screening,
    denominator_screening
  ) %>%
  calculate_prop()

finaldata_primary_7 <- final_table(mortality_data, lung_surr7)

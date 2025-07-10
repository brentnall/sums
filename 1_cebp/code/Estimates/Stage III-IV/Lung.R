### LUNG- Sensitivity (Stage III-IV) Estimates
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



# Surrogate 1. ------------------------------------------------------
surrogate_data_1 <-
  filter_surrogate_data(surrogate_data,
                        "Lung",
                        "1./4. Absolute incidence of late / early stage cancer")

#Czech study
czech_study <- surrogate_data_1 %>%
  filter(trial_acronym == "Czech study") %>%
  # late stage definition: stage III (reported)
  filter(stage_category == "Stage III") %>%
  # follow-up time = mortality outcome follow-up not available, closest time point is 6y from rnd
  filter(period_after_rnd == "6 years") %>%
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

#DLCST
DLCST <- surrogate_data_1 %>%
  filter(trial_acronym == "DLCST") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "DLCST")]) %>%
  # late stage definition: Stage III-IV (reported)
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


#	Johns Hopkins
Johns_Hopkins <- surrogate_data_1 %>%
  filter(trial_acronym == "Johns Hopkins") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "Johns Hopkins")]) %>%
  filter(period_after_rnd == mortality_data$period_after_rnd[mortality_data$trial_acronym ==
                                                               "Johns Hopkins"]) %>%
  filter(stage_category == "Stage 3" |
           stage_category == "Stage 4") %>%
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
  mutate(trial_acronym = "Johns Hopkins", stage_category = "Stage III-IV")

#UKLS
UKLS <- surrogate_data_1 %>%
  filter(trial_acronym == "UKLS") %>%
  # 4y FU halfway timepoint from randomisation to mortality timepoint
  filter(period_after_rnd == "4 years") %>%
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

#NELSON
NELSON <- surrogate_data_1 %>%
  filter(trial_acronym == "NELSON") %>%
  # late stage definition: Stage III-IV (reported)
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

#Lung Screening Study (LSS)
LSS <- surrogate_data_1 %>%
  filter(trial_acronym == "Lung Screening Study (LSS)") %>%
  # definition: baseline, 1 round of screening and interval cancers
  filter(`Effective number of screening rounds` == "2") %>%
  #late stage definition: stage III-IV (reported)
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

#Memorial Sloan-Kettering
MSK <- surrogate_data_1 %>%
  filter(trial_acronym == "Memorial Sloan-Kettering") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "Memorial Sloan-Kettering")]) %>%
  filter(period_after_rnd == mortality_data$period_after_rnd[mortality_data$trial_acronym ==
                                                               "Memorial Sloan-Kettering"]) %>%
  filter(stage_category == "Stage 3" |
           stage_category == "Stage 4") %>%
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
  mutate(trial_acronym = "Memorial Sloan-Kettering", stage_category = "Stage III-IV")

#DANTE
DANTE <- surrogate_data_1 %>%
  filter(trial_acronym == "DANTE") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "DANTE")]) %>%
  filter(
    stage_category == "Stage IIIA" |
      stage_category == "Stage IIIB" | stage_category == "Stage IV"
  ) %>%
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
  mutate(trial_acronym = "DANTE", stage_category = "Stage III-IV")

#LUSI
LUSI <- surrogate_data_1 %>%
  filter(trial_acronym == "LUSI") %>%
  filter(subgroup == "Whole trial") %>%
  # 7 years (halfway point between end of screening and mortality)
  filter(period_after_rnd == "7 years") %>%
  filter(stage_category == "Stage III-IV") %>%
  #unadjusted relative risk > adjusted HR
  filter(type_indicator == "Relative risk") %>%
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



#ITALUNG
ITALUNG <- surrogate_data_1 %>%
  filter(trial_acronym == "ITALUNG") %>%
  # follow-up time = mortality outcome follow-up not available, closest time point (prior to mortality) is 31/12/2013
  filter(end_fu_date == "41639") %>%
  #late stage defintion: Stage III-IV (reported)
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


#	Mayo Lung Project
Mayo <- surrogate_data_1 %>%
  filter(trial_acronym == "Mayo Lung Project") %>%
  # follow-up time = mortality outcome follow-up not available, closest time point is 01/07/1983
  filter(end_fu_date == "30498") %>%
  #late stage definition: stage III-IV (reported)
  filter(stage_category == "Stage III-IV") %>%
  filter(!str_detect(outcome_description, "stage III-IV at diagnosis or")) %>%
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

# MILD
MILD <- surrogate_data_1 %>%
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
  filter(stage_category == "Stage III" |
           stage_category == "Stage IV") %>%
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
  mutate(trial_acronym = "MILD", stage_category = "Stage III-IV")


#NLST
NLST <- surrogate_data_1 %>%
  filter(trial_acronym == "NLST") %>%
  filter(subgroup == "Whole trial") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "NLST")]) %>%
  filter(period_after_rnd == mortality_data$period_after_rnd[mortality_data$trial_acronym ==
                                                               "NLST"]) %>%
  #stage definition: stage III-IV
  filter(stage_category == "Stage III" |
           stage_category == "Stage IV") %>%
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
  mutate(trial_acronym = "NLST", stage_category = "Stage III-IV")

#PLCO
PLCO <- surrogate_data_1 %>%
  filter(trial_acronym == "PLCO (Lung)") %>%
  #late stage definition: stage III+IV non-small lung cancer; 7years FU
  filter(period_after_rnd == "7")

lung_surr1 <-
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
    UKLS
  )

lung_surr1_rate <- lung_surr1 %>%
  filter(type_indicator == "Rate ratio") %>%
  calculate_logRR()

lung_surr1_risk <- lung_surr1 %>%
  filter(type_indicator == "Relative risk") %>%
  calculate_logrel()

PLCO <- PLCO %>%
  calculate_logrel() %>%
  replace_with_reported() %>%
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
    yi,
    vi,
    sei,
    zi,
    pval,
    ci.lb,
    ci.ub
  )



lung_surr1 <- rbind(lung_surr1_rate, lung_surr1_risk, PLCO)

finaldata_primary_1 <- final_table(mortality_data, lung_surr1)


# Surrogate 3. -------------------------------------------------------------
surrogate_data_3 <-
  filter_surrogate_data(surrogate_data,
                        "Lung",
                        "3. % of target cancer diagnosed at late stage")

#Czech study
czech_study <- surrogate_data_3 %>%
  filter(trial_acronym == "Czech study") %>%
  # follow-up time = mortality outcome follow-up not available, closest time point is 6y from rnd
  filter(period_after_rnd == "6 years") %>%
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
    type_indicator
  )

#DLCST
DLCST <- surrogate_data_3 %>%
  filter(trial_acronym == "DLCST") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "DLCST")]) %>%
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
    type_indicator
  )

#	Johns Hopkins
# calculate surrogate 3 using numbers from surrogate 1 file
Johns_Hopkins <- surrogate_data_1 %>%
  filter(trial_acronym == "Johns Hopkins") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "Johns Hopkins")]) %>%
  filter(period_after_rnd == mortality_data$period_after_rnd[mortality_data$trial_acronym ==
                                                               "Johns Hopkins"]) %>%
  filter(stage_category != "Stage 2-4") %>%
  #calculate all cancers diagnosed by arm and use as denominator (all stages + missing)
  mutate(denominator_screening = sum(as.numeric(numerator_screening)) +
           as.numeric(missing_screening)) %>%
  mutate(denominator_comparator = sum(as.numeric(numerator_comparator)) +
           as.numeric(missing_comparator))

Johns_Hopkins_denominator_s <-
  Johns_Hopkins$denominator_screening[1]
Johns_Hopkins_denominator_c <-
  Johns_Hopkins$denominator_comparator[1]

Johns_Hopkins <-  surrogate_data_1 %>%
  filter(trial_acronym == "Johns Hopkins") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "Johns Hopkins")]) %>%
  filter(period_after_rnd == mortality_data$period_after_rnd[mortality_data$trial_acronym ==
                                                               "Johns Hopkins"]) %>%
  # late stage definition: stage III-IV (not reported as 'late stage', but estimate available)
  filter(stage_category == "Stage 3" |
           stage_category == "Stage 4") %>%
  summarise(
    screening_test = first(screening_test),
    numerator_screening = sum(as.numeric(numerator_screening)),
    numerator_comparator = sum(as.numeric(numerator_comparator)),
    denominator_screening = Johns_Hopkins_denominator_s,
    denominator_comparator = Johns_Hopkins_denominator_c,
    type_indicator = first(type_indicator),
    reported_definition_late = first(reported_definition_late),
    reported_definition_early = first(reported_definition_early)
  ) %>%
  mutate(
    trial_acronym = "Johns Hopkins",
    stage_category = "Stage III-IV",
    outcome_description = "% stage 3-4 lung cancers;  9 years of follow-up after study entry;  Screening vs another type of screening"
  )


#PLCO
PLCO <- surrogate_data_3 %>%
  filter(trial_acronym == "PLCO (Lung)") %>%
  #late stage definition: stage III+IV non-small lung cancer
  filter(stage_category == "III" | stage_category == "IV") %>%
  summarise(
    screening_test = first(screening_test),
    numerator_screening = sum(as.numeric(numerator_screening)),
    numerator_comparator = sum(as.numeric(numerator_comparator)),
    denominator_screening = mean(as.numeric(denominator_screening)),
    denominator_comparator = mean(as.numeric(denominator_comparator)),
    type_indicator = first(type_indicator),
    reported_definition_late = first(reported_definition_late),
    reported_definition_early = first(reported_definition_early)
  ) %>%
  mutate(outcome_description = "% non-small cell lung cancers diagnosed at stage III-IV; 13 years of follow-up or December 31st, 2009") %>%
  mutate(trial_acronym = "PLCO (Lung)", stage_category = "Stage III-IV (non-small)")


#UKLS
UKLS <- surrogate_data_3 %>%
  filter(trial_acronym == "UKLS") %>%
  # 4y FU halfway timepoint from randomisation to mortality timepoint
  filter(period_after_rnd == "4 years") %>%
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
    type_indicator
  )


#NELSON
NELSON <- surrogate_data_3 %>%
  filter(trial_acronym == "NELSON") %>%
  # late stage definition: Stage III-IV (reported)
  filter(stage_category == "Stage III-IV") %>%
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
    type_indicator
  )


#Lung Screening Study (LSS)
LSS <- surrogate_data_3 %>%
  filter(trial_acronym == "Lung Screening Study (LSS)") %>%
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
    type_indicator
  )


#Memorial Sloan-Kettering
# calculate surrogate 3 using numbers from surrogate 1 file
MSK <- surrogate_data_1 %>%
  filter(trial_acronym == "Memorial Sloan-Kettering") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "Memorial Sloan-Kettering")]) %>%
  filter(period_after_rnd == mortality_data$period_after_rnd[mortality_data$trial_acronym ==
                                                               "Memorial Sloan-Kettering"]) %>%
  filter(stage_category != "Stage 2-4") %>%
  #calculate all cancers diagnosed by arm and use as denominator (all stages + missing)
  mutate(denominator_screening = sum(as.numeric(numerator_screening)) +
           as.numeric(missing_screening)) %>%
  mutate(denominator_comparator = sum(as.numeric(numerator_comparator)) +
           as.numeric(missing_comparator))

MSK_denominator_s <- MSK$denominator_screening[1]
MSK_denominator_c <- MSK$denominator_comparator[1]

MSK <-  surrogate_data_1 %>%
  filter(trial_acronym == "Memorial Sloan-Kettering") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "Memorial Sloan-Kettering")]) %>%
  filter(period_after_rnd == mortality_data$period_after_rnd[mortality_data$trial_acronym ==
                                                               "Memorial Sloan-Kettering"]) %>%
  filter(stage_category == "Stage 3" |
           stage_category == "Stage 4") %>%
  summarise(
    screening_test = first(screening_test),
    numerator_screening = sum(as.numeric(numerator_screening)),
    numerator_comparator = sum(as.numeric(numerator_comparator)),
    denominator_screening = MSK_denominator_s,
    denominator_comparator = MSK_denominator_c,
    type_indicator = first(type_indicator),
    reported_definition_late = first(reported_definition_late),
    reported_definition_early = first(reported_definition_early)
  ) %>%
  mutate(
    trial_acronym = "Memorial Sloan-Kettering",
    stage_category = "Stage III-IV",
    outcome_description = "% Absolute incidence of stage 4 lung cancers; follow-up until 9 years since study entry; Dual screen vs X-ray only"
  )


#DANTE
# calculate surrogate 3 using numbers from surrogate 1 file
DANTE <- surrogate_data_1 %>%
  filter(trial_acronym == "DANTE") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "DANTE")]) %>%
  filter(stage_category != "Stage IA") %>%
  filter(stage_category != "Stage IB") %>%
  filter(stage_category != "Stage II-IV") %>%
  #calculate all cancers diagnosed by arm and use as denominator (all stages + missing)
  mutate(denominator_screening = sum(as.numeric(numerator_screening)) +
           as.numeric(missing_screening)) %>%
  mutate(denominator_comparator = sum(as.numeric(numerator_comparator)) +
           as.numeric(missing_comparator))


DANTE_denominator_s <- DANTE$denominator_screening[1]
DANTE_denominator_c <- DANTE$denominator_comparator[1]

DANTE <-  surrogate_data_1 %>%
  filter(trial_acronym == "DANTE") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "DANTE")]) %>%
  # late stage definition: stage III-IV (not reported as 'late stage', but estimate available)
  filter(
    stage_category == "Stage IIIA" |
      stage_category == "Stage IIIB" | stage_category == "Stage IV"
  ) %>%
  summarise(
    screening_test = first(screening_test),
    numerator_screening = sum(as.numeric(numerator_screening)),
    numerator_comparator = sum(as.numeric(numerator_comparator)),
    denominator_screening = DANTE_denominator_s,
    denominator_comparator = DANTE_denominator_c,
    type_indicator = first(type_indicator),
    reported_definition_late = first(reported_definition_late),
    reported_definition_early = first(reported_definition_early)
  ) %>%
  mutate(
    trial_acronym = "DANTE",
    stage_category = "Stage III-IV",
    outcome_description = "% stage III-IV lung cancers diagnosed up to 15 May 2013;  screening vs another type of screening"
  )



#LUSI
# calculate surrogate 3 using numbers from surrogate 1 file
LUSI <- surrogate_data_1 %>%
  filter(trial_acronym == "LUSI") %>%
  filter(subgroup == "Whole trial") %>%
  # 7 years (halfway point between end of screening and mortality)
  filter(period_after_rnd == "7 years") %>%
  filter(stage_category != "Stage II-IV") %>%
  filter(stage_category != "Stage I") %>%
  #calculate all cancers diagnosed by arm and use as denominator (all stages + missing)
  mutate(denominator_screening = sum(as.numeric(numerator_screening)) +
           as.numeric(missing_screening)) %>%
  mutate(denominator_comparator = sum(as.numeric(numerator_comparator))) %>%
  # late stage definition: stage II+ (reported)
  filter(stage_category == "Stage III-IV") %>%
  mutate(outcome_description = "% advanced stage lung cancer (stage III-IV) (per 100 people);  7 years after randomisation; Screening vs controls") %>%
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
    type_indicator
  )


#ITALUNG
ITALUNG <- surrogate_data_3 %>%
  filter(trial_acronym == "ITALUNG") %>%
  # follow-up time = mortality outcome follow-up not available, closest time point (prior to mortality) is 31/12/2013
  filter(end_fu_date == "41639") %>%
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
    type_indicator
  )

#	Mayo Lung Project
Mayo <- surrogate_data_3 %>%
  filter(trial_acronym == "Mayo Lung Project") %>%
  # follow-up time = mortality outcome follow-up not available, closest time point is 01/07/1983
  filter(end_fu_date == "30498") %>%
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
    type_indicator
  )


# MILD - only subgroup
# calculate surrogate 3 using numbers from surrogate 1 file
MILD <- surrogate_data_1 %>%
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
  filter(stage_category != "Stage II-IV") %>%
  #calculate all cancers diagnosed by arm and use as denominator (all stages + missing)
  mutate(denominator_screening = sum(as.numeric(numerator_screening)) +
           as.numeric(missing_screening)) %>%
  mutate(denominator_comparator = sum(as.numeric(numerator_comparator)) +
           as.numeric(missing_comparator))


MILD_denominator_s <- MILD$denominator_screening[1]
MILD_denominator_c <- MILD$denominator_comparator[1]

MILD <-  surrogate_data_1 %>%
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
  filter(stage_category == "Stage III" |
           stage_category == "Stage IV") %>%
  summarise(
    screening_test = first(screening_test),
    numerator_screening = sum(as.numeric(numerator_screening)),
    numerator_comparator = sum(as.numeric(numerator_comparator)),
    denominator_screening = MILD_denominator_s,
    denominator_comparator = MILD_denominator_c,
    type_indicator = first(type_indicator),
    reported_definition_late = first(reported_definition_late),
    reported_definition_early = first(reported_definition_early)
  ) %>%
  mutate(
    trial_acronym = "MILD",
    stage_category = "Stage III-IV",
    outcome_description = "% stage III-IV lung cancers; follow-up until June 2018; LDCT (annual and biennial combined) vs control"
  )



#NLST
# calculate surrogate 3 using numbers from surrogate 1 file
NLST <- surrogate_data_1 %>%
  filter(trial_acronym == "NLST") %>%
  filter(subgroup == "Whole trial") %>%
  #  follow-up time = mortality outcome follow-up
  filter(end_fu_date == mortality_data$end_fu_date[(mortality_data$trial_acronym ==
                                                      "NLST")]) %>%
  filter(period_after_rnd == mortality_data$period_after_rnd[mortality_data$trial_acronym ==
                                                               "NLST"]) %>%
  filter(!str_detect(stage_category, "A")) %>%
  filter(!str_detect(stage_category, "B")) %>%
  #calculate all cancers diagnosed by arm and use as denominator (all stages + missing)
  mutate(denominator_screening = sum(as.numeric(numerator_screening)) +
           as.numeric(missing_screening)) %>%
  mutate(denominator_comparator = sum(as.numeric(numerator_comparator)) +
           as.numeric(missing_comparator)) %>%
  filter(stage_category == "Stage III" |
           stage_category == "Stage IV") %>%
  summarise(
    screening_test = first(screening_test),
    numerator_screening = sum(as.numeric(numerator_screening)),
    numerator_comparator = sum(as.numeric(numerator_comparator)),
    denominator_screening = mean(as.numeric(denominator_screening)),
    denominator_comparator = mean(as.numeric(denominator_comparator)),
    type_indicator = first(type_indicator),
    reported_definition_late = first(reported_definition_late),
    reported_definition_early = first(reported_definition_early)
  ) %>%
  mutate(
    trial_acronym = "NLST",
    stage_category = "Stage III-IV",
    outcome_description = "%stage III-IV lung cancers (per 100 people randomised); follow-up until 31/12/2009; LDCT vs CXR"
  )




lung_surr3 <-
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
    PLCO,
    UKLS,
    NLST
  )

lung_surr3 <- calculate_logrel(lung_surr3)

finaldata_primary_3 <- final_table(mortality_data, lung_surr3)

  
setwd("../../../code")
source("Allcancer_functions.R")

# Load data ---------------------------------------------------------------
setwd("../data")
mortality_data <- read_excel("SUMS_Results_Master_final.xlsx", sheet="Master_Mortality_Results")
surrogate_data <- read_excel("SUMS_Results_Master_final.xlsx", sheet="Master_Surrogate_Results")

mortality_data <- tibble::as_tibble(mortality_data)
surrogate_data <- tibble::as_tibble(surrogate_data)

setwd("../code/Estimates")
source("Breast.R")

# PRIMARY ANALYSIS ESTIMATES ----------------------------------------------
#select trials that have a midpoint estimate (earlier than mortality FU)
mid_timepoint <- finaldata_primary_1 %>%
  filter(
    trial_acro_sm == "HIP" |
      trial_acro_sm == "Stockholm Breast" |
      trial_acro_sm == "Russia WHO Breast" |
      trial_acro_sm == "Gothenburg Breast" |
      trial_acro_sm == "UK Age"
  ) %>%
  mutate(timepoint = "Earlier than mortality")


# SAME FU AS MORTALITY ----------------------------------------------------

# HIP # no estimate available

# Stockholm # no estimate available

# RussiaWHO # no estimate available

# Gothenburg # no estimate available

# UKAge #no estimate available

#for all other trials we used surrogate measured at time of mortality
main_timepoint <- finaldata_primary_1 %>%
  filter(
    trial_acro_sm != "HIP" &
      trial_acro_sm != "Stockholm Breast" &
      trial_acro_sm != "Russia WHO Breast" &
      trial_acro_sm != "Gothenburg Breast" &
      trial_acro_sm != "UK Age"
  ) %>%
  mutate(timepoint = "Mortality")


# ONGOING SCREENING --------------------------------------------------------
HIP <-
  filter_surrogate_data(surrogate_data,
                        "Breast",
                        "1./4. Absolute incidence of late / early stage cancer") %>%
  filter(trial_acronym == "HIP") %>%
  filter(period_after_rnd == "1 year") %>%
  filter(stage_category == "Node-positive") %>%
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


Stockholm <-
  filter_surrogate_data(surrogate_data,
                        "Breast",
                        "1./4. Absolute incidence of late / early stage cancer") %>%
  filter(trial_acronym == "Stockholm Breast") %>%
  filter(stage_category == "II-IV") %>%
  filter(str_detect(outcome_description, "1981-1984")) %>%
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


RussiaWHO <-
  filter_surrogate_data(surrogate_data,
                        "Breast",
                        "1./4. Absolute incidence of late / early stage cancer") %>%
  filter(trial_acronym == "Russia WHO Breast") %>%
  # late stage definition: 	Node positive (N1-N2)
  filter(stage_category == "Node positive (N1-N2)") %>%
  # follow-up up to 1990
  filter(end_fu_date == "33238") %>%
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

# Gothenburg # no estimate available

UKAge <-
  filter_surrogate_data(surrogate_data,
                        "Breast",
                        "1./4. Absolute incidence of late / early stage cancer") %>%
  filter(trial_acronym == "UK Age") %>%
  # follow-up time up to 1999
  filter(end_fu_date == "36525") %>%
  filter(type_indicator == "Rate ratio") %>%
  # late stage definition: lymph-node positive
  filter(stage_category == "Node positive") %>%
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

# Mumbai #no estimate available

# Edinburgh #no earlier timepoint than chosen (mortality measured when screening still ongoing)
#
# Trivandrum #no estimate available

Twocounty <-
  filter_surrogate_data(surrogate_data,
                        "Breast",
                        "1./4. Absolute incidence of late / early stage cancer") %>%
  filter(trial_acronym == "Two-County") %>%
  filter(staging_system == "Numerical stage") %>%
  filter(str_detect(outcome_description, "3 year")) %>%
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

# Shanghai #no estimate available

CNBSS1 <-
  filter_surrogate_data(surrogate_data,
                        "Breast",
                        "1./4. Absolute incidence of late / early stage cancer") %>%
  filter(trial_acronym == "CNBSS-1") %>%
  #  follow-up time = mortality outcome follow-up (5 yrs from rnd)
  filter(period_after_rnd == "5 years") %>%
  # late stage defiiniion: Stage III-IV (reported)
  filter(stage_category == "Node-positive") %>%
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

CNBSS2 <-
  filter_surrogate_data(surrogate_data,
                        "Breast",
                        "1./4. Absolute incidence of late / early stage cancer") %>%
  filter(trial_acronym == "CNBSS-2") %>%
  #  follow-up time = mortality outcome follow-up (5 yrs from rnd)
  filter(period_after_rnd == "5 years") %>%
  # late stage defiiniion: Stage III-IV (reported)
  filter(stage_category == "Node-positive") %>%
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

MMSTI <-
  filter_surrogate_data(surrogate_data,
                        "Breast",
                        "1./4. Absolute incidence of late / early stage cancer") %>%
  filter(trial_acronym == "MMST I") %>%
  filter(str_detect(outcome_description, "5 years")) %>%
  calculate_logrel() %>%
  mutate(yi = log(as.numeric(calculated_estimate_unadjusted))) %>%
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


breast_surr1 <-
  rbind(HIP, Stockholm, RussiaWHO, UKAge, Twocounty, CNBSS1, CNBSS2)

breast_surr1_rate <- breast_surr1 %>%
  filter(type_indicator == "Rate ratio") %>%
  calculate_logRR()

breast_surr1_risk <- breast_surr1 %>%
  filter(type_indicator == "Relative risk") %>%
  calculate_logrel()

breast_surr1 <- rbind(breast_surr1_rate, breast_surr1_risk, MMSTI)

screening_timepoint <- final_table(mortality_data, breast_surr1) %>%
  mutate(timepoint = "Intervention phase")

# COMBINE ALL ESTIMATES ---------------------------------------------------
finaldata_subgroup_1 <-
  rbind(mid_timepoint, main_timepoint, screening_timepoint)

### CERVICAL - Sensitivity (Stage III-IV) Estimates
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
surrogate_data_1 <- filter_surrogate_data(surrogate_data, "Cervical","1./4. Absolute incidence of late / early stage cancer")

## No trials with stage III-IV as late stage definition.


# Surrogate 3. -------------------------------------------------------------
surrogate_data_3 <- filter_surrogate_data(surrogate_data, "Cervical", "3. % of target cancer diagnosed at late stage")

## No trials with stage III-IV as late stage definition.
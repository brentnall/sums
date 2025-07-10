### ALL CANCER FUNCTIONS

# Load packages -----------------------------------------------------------
library(readxl)
library(tidyverse)
library(metafor)


# Filter mortality data -------------------------------------------------------------
filter_mortality_data <- function(mortality_data, cancer) {
  #cancer type
  mortality_data <- mortality_data %>%
    filter(cancer_type == cancer) %>%
  #main timepoint for mortality
    filter(main_timepoint_RoB == "yes")
  
  return(mortality_data)
}



# Filter surrogate data ------------------------------------------------------------
filter_surrogate_data <-
  function(surrogate_data, cancer, surrogate_x) {
    surrogate_data <- surrogate_data %>%
      #cancer type
      filter(cancer_type == cancer) %>%
      #surrogate
      filter(surrogate == surrogate_x)
    
    return(surrogate_data)
    
  }

#Count N trials available (need minimum 3)
min_N_trials <- function(surrogate_data) {
  N <- n_distinct(surrogate_data$trial_acronym)
  return(N)
}


# Effect Estimates --------------------------------------------------------
# Calculate log rate ratio and 95% confidence intervals
calculate_logRR <- function(data){
  data <- data %>%
    mutate(numerator_screening=as.numeric(numerator_screening)) %>%
    mutate(numerator_comparator=as.numeric(numerator_comparator)) %>%
    mutate(denominator_screening=as.numeric(denominator_screening)) %>%
    mutate(denominator_comparator=as.numeric(denominator_comparator))
  
  data <- escalc(data=data, 
              x1i=numerator_screening, x2i=numerator_comparator, 
              t1i=denominator_screening, t2i=denominator_comparator, 
              measure="IRR")
  
  data <- summary.escalc(data)
  
  return(data)
}

# Calculate log relative risk and 95% confidence intervals
calculate_logrel <- function(data){
  data <- data %>%
    mutate(numerator_screening=as.numeric(numerator_screening)) %>%
    mutate(numerator_comparator=as.numeric(numerator_comparator)) %>%
    mutate(denominator_screening=as.numeric(denominator_screening)) %>%
    mutate(denominator_comparator=as.numeric(denominator_comparator))
  
  data <- escalc(data=data, 
                 ai=numerator_screening, ci=numerator_comparator, 
                 n1i=denominator_screening, n2i=denominator_comparator, 
                 measure="RR")
  
  data <- summary.escalc(data)
  
  return(data)
}


# Calculate proportion and 95% confidence intervals
calculate_prop <- function(data){
  data <- data %>%
    mutate(numerator_screening=as.numeric(numerator_screening)) %>%
    mutate(denominator_screening=as.numeric(denominator_screening)) 
    
  data <- escalc(data=data, 
                 xi=numerator_screening, 
                 ni=denominator_screening, 
                 measure="PR")
  
  data <- summary.escalc(data)
  
  return(data)
}

#Replace estimate with reported estimate if available
replace_with_reported <- function(data){
  data <- data%>%
    mutate(yi = ifelse(reported_estimate_unadjusted!="NR", log(as.numeric(reported_estimate_unadjusted)), yi)) %>%
    mutate(ci.lb = ifelse(reported_L95CI_unadjusted!="NR", log(as.numeric(reported_L95CI_unadjusted)), ci.lb)) %>%
    mutate(ci.ub = ifelse(reported_U95CI_unadjusted!="NR", log(as.numeric(reported_U95CI_unadjusted)), ci.ub)) %>%
    #obtain SE from 95% CIs
    mutate(sei = (ci.ub-yi)/1.96) %>%
    # calculate variance
    mutate(vi = sei^2)
  
  return(data)
}


# Create final table for analysis with mortality and surrogate estimates
final_table <- function(mortality_data, surrogate_data){
  mortality_data <- mortality_data %>%
    select(trial_acronym, screening_test, outcome_description, numerator_screening, numerator_comparator, denominator_screening, denominator_comparator, yi, vi, sei, ci.lb, ci.ub)
  
  data <- full_join(mortality_data, surrogate_data, by=c("trial_acronym", "screening_test"))
  
  colnames(data) <- sub(".x", "_m", colnames(data))
  colnames(data) <- sub(".y", "_s", colnames(data))
  
  return(data)
}


# Analysis ----------------------------------------------------------------
# - regression (through origin), weights: inverse variance
# - weighted correlation
# - returns line for results table

analysis <- function(finaldata_primary, cancer, surrogate, analysis, bootfunction){

  finaldata_primary <- finaldata_primary %>%
    filter(!is.na(yi_s))
    
  # Regression (inverse variance weights)
  m_reg <- lm(formula = yi_m ~ yi_s, data = finaldata_primary, weights = weight_var)
  m_reg_sum <- summary(m_reg)
  
  cis <- bootfunction(finaldata_primary)
  
  # correlation
  correlation <- weights::wtd.cor(finaldata_primary$yi_m, finaldata_primary$yi_s, weight=finaldata_primary$weight_var)
  
  results_surrogate <- data.frame(Cancer=cancer,
                                  Surrogate=surrogate,
                                  Analysis=analysis,
                                  Ntrials=finaldata_primary$N_trials[1],
                                  Intercept = as.numeric(m_reg_sum$coefficients[1]),
                                  Intercept_lb_ci = cis[1,1], 
                                  Intercept_ub_ci = cis[2,1],
                                  Slope_lb_ci = cis[1,2],
                                  Slope_ub_ci = cis[2,2],
                                  Slope = as.numeric(m_reg_sum$coefficients[2]),
                                  Correlation=correlation[1],
                                  Correlation_lb_ci = cis[1,3],
                                  Correlation_ub_ci = cis[2,3],
                                  Rsquared=m_reg_sum$r.squared,
                                  Rsquared_lb_ci = cis[1,4],
                                  Rsquared_ub_ci = cis[2,4]
    )
  
  #round results to 2 decimal places
  results_surrogate <- results_surrogate %>%
    mutate(across(where(is.numeric), ~ round(., 2)))
}
  
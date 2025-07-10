library(boot)


# Bootstrap where trials are unit of resample - Main Analysis -------------

bootstrap_trial <- function(finaldata_primary_1, myR = 50000) {
  
  ## stats function
  mysim_cluster <- function(finaldata_primary_1, indices) {
    mysimdta <- finaldata_primary_1[indices, ]
    
    m_reg <-
      summary(lm(
        formula = yi_m ~ yi_s,
        data = mysimdta,
        weights = weight_var
      ))
    
    # weighted correlation(weights=inverse variance)
    correlation <-
      tryCatch(
        weights::wtd.cor(mysimdta$yi_m, mysimdta$yi_s, weight = mysimdta$weight_var),
        error = function(e) {
          return(NA)
        }
      )
    
    # return line for results table
    myres <- c(
      Intercept = as.numeric(m_reg$coefficients[1]),
      Slope = as.numeric(m_reg$coefficients[2]),
      Correlation = correlation[1],
      Rsquared = m_reg$r.squared
    )
    
    myres
    
  }
  
  
  mysim_cluster2 <-
    function(intrial,
             indices,
             finaldata_primary_1,
             inlist) {
      indices2 <- unlist(inlist[indices])
      
      mysimdta <- finaldata_primary_1[indices2,]
      
      m_reg <-
        summary(lm(
          formula = yi_m ~ yi_s,
          data = mysimdta,
          weights = weight_var
        ))
      
      # weighted correlation(weights=inverse variance)
      correlation <-
        tryCatch(
          weights::wtd.cor(mysimdta$yi_m, mysimdta$yi_s, weight = mysimdta$weight_var),
          error = function(e) {
            return(NA)
          }
        )
      
      # return line for results table
      myres <- c(
        Intercept = as.numeric(m_reg$coefficients[1]),
        Slope = as.numeric(m_reg$coefficients[2]),
        Correlation = correlation[1],
        Rsquared = m_reg$r.squared
      )
      
      myres
      
    }
  
    
  set.seed(20231124)
  

  ## row numbers by trial in a list length n trials
  ## triallist<-list(c(1),c(2,3),c(4))
  ## trial numbers
  trialnum <- 1:(length(unique(finaldata_primary_1$trial_acro_sm)))
  
  if (length(trialnum) < nrow(finaldata_primary_1)) {
    triallist <-
      tapply(1:nrow(finaldata_primary_1), finaldata_primary_1$trial_acro_sm, function(ind)
        list(ind))
    
    mybootstrap <-
      boot(
        trialnum,
        mysim_cluster2,
        R = myR,
        finaldata_primary_1 = finaldata_primary_1,
        inlist = triallist
      )
    
  }    else {
    mybootstrap <- boot(finaldata_primary_1, mysim_cluster, R = myR)
  }
  
  
  # for bootstraps where correlation=NA (because too few trials) give slope a value of either -inf or inf, and correlation a value of -1 or 1
  slope_NA <- c(-99999, 99999)
  corr_NA <- c(-1, 1)
  
  mybootstrap$t[, 2] <-
    ifelse(is.na(mybootstrap$t[, 3]),
           sample(slope_NA, size = sum(is.na(mybootstrap$t[, 3])), replace = T),
           mybootstrap$t[, 2])
  
  mybootstrap$t[, 3] <-
    ifelse(is.na(mybootstrap$t[, 3]),
           sample(corr_NA, size = sum(is.na(mybootstrap$t[, 3])), replace = T),
           mybootstrap$t[, 3])
  
  
  #intercept
  ci_intercept <- boot.ci(mybootstrap, index = 1)
  #slope
  ci_slope <- boot.ci(mybootstrap, index = 2)
  #correlation
  ci_correlation <-
    tryCatch(
      boot.ci(mybootstrap, index = 3),
      error = function(e) {
        return(c(-1, 1))
      }
    )
  #R^2
  ci_rsquared <-
    tryCatch(
      boot.ci(mybootstrap, index = 4),
      error = function(e) {
        return(c(0, 1))
      }
    )
  
  ##Percantile CIs
  myci <- matrix(nrow = 2, ncol = 4)
  
  #intercept
  myci[1, 1] <- ci_intercept$percent[4]
  myci[2, 1] <- ci_intercept$percent[5]
  
  #slope
  myci[1, 2] <- ci_slope$percent[4]
  myci[2, 2] <- ci_slope$percent[5]
  
  #correlation
  myci[1, 3] <-
    ifelse(is(ci_correlation, "bootci"),
           ci_correlation$percent[4],
           ci_correlation[1])
  myci[2, 3] <-
    ifelse(is(ci_correlation, "bootci"),
           ci_correlation$percent[5],
           ci_correlation[2])
  
  #R^2
  myci[1, 4] <-
    ifelse(is(ci_rsquared, "bootci"), ci_rsquared$percent[4], ci_rsquared[1])
  myci[2, 4] <-
    ifelse(is(ci_rsquared, "bootci"), ci_rsquared$percent[5], ci_rsquared[2])
  
  ##print results
  return(myci)
}





# Bootstrap where subgroup (cancer/test type) is unit of resample - Exploratory FE Meta-Analysis -------------


#Surrogate 1
bootstrap_1 <- function(finaldata_primary_1, myR = 50000) {
  ## stats function
  mysim_cluster <- function(finaldata_primary_1, indices) {
    mysimdta <- finaldata_primary_1[indices, ]
    
    m_reg <-
      summary(
        lm(
          formula = Effect_Mortality ~ `Effect_Surrogate 1 (vi_m)`,
          data = mysimdta,
          weights = weight_var
        )
      )
    
    # weighted correlation(weights=inverse variance)
    correlation <-
      tryCatch(
        weights::wtd.cor(
          mysimdta$Effect_Mortality,
          mysimdta$`Effect_Surrogate 1 (vi_m)`,
          weight = mysimdta$weight_var
        ),
        error = function(e) {
          return(NA)
        }
      )
    
    # return line for results table
    myres <- c(
      Intercept = as.numeric(m_reg$coefficients[1]),
      Slope = as.numeric(m_reg$coefficients[2]),
      Correlation = correlation[1],
      Rsquared = m_reg$r.squared
    )
    
    myres
    
  }
  
  mysim_cluster2 <-
    function(intrial,
             indices,
             finaldata_primary_1,
             inlist) {
      indices2 <- unlist(inlist[indices])
      
      mysimdta <- finaldata_primary_1[indices2, ]
      
      m_reg <-
        summary(
          lm(
            formula = Effect_Mortality ~ `Effect_Surrogate 1 (vi_m)`,
            data = mysimdta,
            weights = weight_var
          )
        )
      
      # weighted correlation(weights=inverse variance)
      correlation <-
        tryCatch(
          weights::wtd.cor(
            mysimdta$Effect_Mortality,
            mysimdta$`Effect_Surrogate 1 (vi_m)`,
            weight = mysimdta$weight_var
          ),
          error = function(e) {
            return(NA)
          }
        )
      
      # return line for results table
      myres <- c(
        Intercept = as.numeric(m_reg$coefficients[1]),
        Slope = as.numeric(m_reg$coefficients[2]),
        Correlation = correlation[1],
        Rsquared = m_reg$r.squared
      )
      
      myres
      
    }
  
  
  set.seed(20231124)
  
  trialnum <- 1:(length(unique(finaldata_primary_1$Subgroup)))
  
  if (length(trialnum) < nrow(finaldata_primary_1)) {
    triallist <-
      tapply(1:nrow(finaldata_primary_1), finaldata_primary_1$Subgroup, function(ind)
        list(ind))
    
    mybootstrap <-
      boot(
        trialnum,
        mysim_cluster2,
        R = myR,
        finaldata_primary_1 = finaldata_primary_1,
        inlist = triallist
      )
    
  }    else {
    mybootstrap <- boot(finaldata_primary_1, mysim_cluster, R = myR)
  }
  
  # for bootstraps where correlation=NA (because too few estimates) give slope a value of either -inf or inf, and correlation a value of -1 or 1
  slope_NA <- c(-99999, 99999)
  corr_NA <- c(-1, 1)
  
  mybootstrap$t[, 2] <-
    ifelse(is.na(mybootstrap$t[, 3]),
           sample(slope_NA, size = sum(is.na(mybootstrap$t[, 3])), replace = T),
           mybootstrap$t[, 2])
  
  mybootstrap$t[, 3] <-
    ifelse(is.na(mybootstrap$t[, 3]),
           sample(corr_NA, size = sum(is.na(mybootstrap$t[, 3])), replace = T),
           mybootstrap$t[, 3])
  
  
  #intercept
  ci_intercept <- boot.ci(mybootstrap, index = 1)
  #slope
  ci_slope <- boot.ci(mybootstrap, index = 2)
  #correlation
  ci_correlation <-
    tryCatch(
      boot.ci(mybootstrap, index = 3),
      error = function(e) {
        return(c(-1, 1))
      }
    )
  #R^2
  ci_rsquared <-
    tryCatch(
      boot.ci(mybootstrap, index = 4),
      error = function(e) {
        return(c(0, 1))
      }
    )
  
  ## Percantile CIs
  myci <- matrix(nrow = 2, ncol = 4)
  
  #intercept
  myci[1, 1] <- ci_intercept$percent[4]
  myci[2, 1] <- ci_intercept$percent[5]
  
  #slope
  myci[1, 2] <- ci_slope$percent[4]
  myci[2, 2] <- ci_slope$percent[5]
  
  #correlation
  myci[1, 3] <-
    ifelse(is(ci_correlation, "bootci"),
           ci_correlation$percent[4],
           ci_correlation[1])
  myci[2, 3] <-
    ifelse(is(ci_correlation, "bootci"),
           ci_correlation$percent[5],
           ci_correlation[2])
  
  #R^2
  myci[1, 4] <-
    ifelse(is(ci_rsquared, "bootci"),
           ci_rsquared$percent[4],
           ci_rsquared[1])
  myci[2, 4] <-
    ifelse(is(ci_rsquared, "bootci"),
           ci_rsquared$percent[5],
           ci_rsquared[2])
  
  ##print results
  return(myci)
}


#Surrogate 3
bootstrap_3 <- function(finaldata_primary_1, myR = 50000) {
  ## stats function
  mysim_cluster <- function(finaldata_primary_1, indices) {
    mysimdta <- finaldata_primary_1[indices, ]
    
    m_reg <-
      summary(
        lm(
          formula = Effect_Mortality ~ `Effect_Surrogate 3 (vi_m)`,
          data = mysimdta,
          weights = weight_var
        )
      )
    
    # weighted correlation(weights=inverse variance)
    correlation <-
      tryCatch(
        weights::wtd.cor(
          mysimdta$Effect_Mortality,
          mysimdta$`Effect_Surrogate 3 (vi_m)`,
          weight = mysimdta$weight_var
        ),
        error = function(e) {
          return(NA)
        }
      )
    
    # return line for results table
    myres <- c(
      Intercept = as.numeric(m_reg$coefficients[1]),
      Slope = as.numeric(m_reg$coefficients[2]),
      Correlation = correlation[1],
      Rsquared = m_reg$r.squared
    )
    
    myres
    
  }
  
  mysim_cluster2 <-
    function(intrial,
             indices,
             finaldata_primary_1,
             inlist) {
      indices2 <- unlist(inlist[indices])
      
      mysimdta <- finaldata_primary_1[indices2, ]
      
      m_reg <-
        summary(
          lm(
            formula = Effect_Mortality ~ `Effect_Surrogate 3 (vi_m)`,
            data = mysimdta,
            weights = weight_var
          )
        )
      
      # weighted correlation(weights=inverse variance)
      correlation <-
        tryCatch(
          weights::wtd.cor(
            mysimdta$Effect_Mortality,
            mysimdta$`Effect_Surrogate 3 (vi_m)`,
            weight = mysimdta$weight_var
          ),
          error = function(e) {
            return(NA)
          }
        )
      
      # return line for results table
      myres <- c(
        Intercept = as.numeric(m_reg$coefficients[1]),
        Slope = as.numeric(m_reg$coefficients[2]),
        Correlation = correlation[1],
        Rsquared = m_reg$r.squared
      )
      
      myres
      
    }
  
  
  set.seed(20231124)
  
  trialnum <- 1:(length(unique(finaldata_primary_1$Subgroup)))
  
  if (length(trialnum) < nrow(finaldata_primary_1)) {
    triallist <-
      tapply(1:nrow(finaldata_primary_1), finaldata_primary_1$Subgroup, function(ind)
        list(ind))
    
    mybootstrap <-
      boot(
        trialnum,
        mysim_cluster2,
        R = myR,
        finaldata_primary_1 = finaldata_primary_1,
        inlist = triallist
      )
    
  }    else {
    mybootstrap <- boot(finaldata_primary_1, mysim_cluster, R = myR)
  }
  
  # for bootstraps where correlation=NA (because too few estimates) give slope a value of either -inf or inf, and correlation a value of -1 or 1
  slope_NA <- c(-99999, 99999)
  corr_NA <- c(-1, 1)
  
  mybootstrap$t[, 2] <-
    ifelse(is.na(mybootstrap$t[, 3]),
           sample(slope_NA, size = sum(is.na(mybootstrap$t[, 3])), replace = T),
           mybootstrap$t[, 2])
  
  mybootstrap$t[, 3] <-
    ifelse(is.na(mybootstrap$t[, 3]),
           sample(corr_NA, size = sum(is.na(mybootstrap$t[, 3])), replace = T),
           mybootstrap$t[, 3])
  
  
  #intercept
  ci_intercept <- boot.ci(mybootstrap, index = 1)
  #slope
  ci_slope <- boot.ci(mybootstrap, index = 2)
  #correlation
  ci_correlation <-
    tryCatch(
      boot.ci(mybootstrap, index = 3),
      error = function(e) {
        return(c(-1, 1))
      }
    )
  #R^2
  ci_rsquared <-
    tryCatch(
      boot.ci(mybootstrap, index = 4),
      error = function(e) {
        return(c(0, 1))
      }
    )
  
  ##Percantile CIs
  myci <- matrix(nrow = 2, ncol = 4)
  
  #intercept
  myci[1, 1] <- ci_intercept$percent[4]
  myci[2, 1] <- ci_intercept$percent[5]
  
  #slope
  myci[1, 2] <- ci_slope$percent[4]
  myci[2, 2] <- ci_slope$percent[5]
  
  #correlation
  myci[1, 3] <-
    ifelse(is(ci_correlation, "bootci"),
           ci_correlation$percent[4],
           ci_correlation[1])
  myci[2, 3] <-
    ifelse(is(ci_correlation, "bootci"),
           ci_correlation$percent[5],
           ci_correlation[2])
  
  #R^2
  myci[1, 4] <-
    ifelse(is(ci_rsquared, "bootci"),
           ci_rsquared$percent[4],
           ci_rsquared[1])
  myci[2, 4] <-
    ifelse(is(ci_rsquared, "bootci"),
           ci_rsquared$percent[5],
           ci_rsquared[2])
  
  ##print results
  return(myci)
}

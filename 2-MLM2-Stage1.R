#!/usr/bin/env Rscript
cat("[BP-DNAme] Running Rscript 2-MLM2-Stage1...\n")
################################################################################
# BP-DNAme -  Summarize detrended BP values and merge with clinical covariates
################################################################################
# Author: Lucile
# Date of creation: September 2021
# Note: 
# Run Rscript 0-* and 1-* before
################################################################################
# Description:
# 
#------------------------------------------------------------------------------#
# Method:
# 
#------------------------------------------------------------------------------#
# Note:
# 
################################################################################
# Paths and parameters
#------------------------------------------------------------------------------#

# Input file (centered individual BP values)#
# generated using Rscript 1-BP-Centering.R
BPfile <-  "Data/N2002_BP_centered.rds"

# Covariates, generated using Rscript 0-Covariates-Preprocessing.R
covFile <- "Data/N2002_covariates.rds"

# Cell-type estimates (ilr-transformed, ready for adjustment in LM)
CCFile <- "Data/EDEN_CC.planet_rpc_ilr.rds"

windows <- c("pregnancy", "earlyP", "lateP")
limit <- 27 # week of gestation that separates earlyP from lateP
time <- "GA" # "GA" (fit autocorrelation structure) or NULL

################################################################################
# R packages
#------------------------------------------------------------------------------#

library(magrittr)

################################################################################
# R function
#------------------------------------------------------------------------------#

MLM2_stage1 <- function(exposure, time = NULL, data, PLOT_AR1 = FALSE){
  
  #--------------------------------------------------------------------------#
  # Two-stage model: first stage
  #--------------------------------------------------------------------------#
  ## Stage one: mixed model 
  ## (assuming 'constant' exposure, repeatedly measured)
  
  if( is.null(time) ){
    
    ### ML regression formula: exposure ~ 1 + (1 | id)
    formula.repeat <- stats::as.formula(paste(exposure, " ~ (1 | id)"))
    
    ### Fit mixed linear model
    cat("Fit mixed linear model without correlation structure \n")
    fit <- lmerTest::lmer(formula = formula.repeat, data = data)
    ### Gather fixed and random effects
    fixed <- as.vector(summary(fit)$coefficients["(Intercept)","Estimate"])
    random <- lme4::ranef(fit) %>% as.data.frame()
    
    # conditional standard deviation:
    # condsd <- random$condsd
    ### Compute individual coefficient (pop fixed + indiv random)
    exposure.estim <- data.frame(id = random$grp,
                                 exposure = random$condval+fixed,
                                 condsd = random$condsd)
    colnames(exposure.estim)[3] <- paste0(exposure,".condsd")
    
  }else{
    ### Fit mixed linear model with time autocorrelation
    cat("Fit mixed linear model with time autocorrelation \n")
    fit <- nlme::lme(fixed = stats::as.formula(paste(exposure, " ~ 1")), 
                     random = ~1|id, 
                     correlation = nlme::corAR1(form = stats::as.formula(paste("~",time))),
                     data = data)
    
    if( PLOT_AR1 ) print( plot(nlme::ACF(fit, resType = "normalized")) )
    ### Gather fixed and random effects
    fixed <- as.numeric( summary(fit)$coefficients$fixed )
    random <- summary(fit)$coefficients$random$id 
    ### Compute individual coefficient (pop fixed + indiv random)
    exposure.estim <- data.frame(id = rownames(random),
                                 exposure = as.numeric(random)+fixed)
    
  }
  
  colnames(exposure.estim)[2] <- exposure

  return( exposure.estim )
  
}

################################################################################
#------------------------------------------------------------------------------#
################################################################################
# Run over windows
################################################################################

for(window in windows){
  
  cat("Time window:", window, "\n")
  # Output file
  outFile <- paste0("Data/N668_BP1+covariates_", window, ".rds")
  print(outFile)
  
  ##############################################################################
  # Compute individual average values for each time window, using a mixed LM
  #----------------------------------------------------------------------------#
  BP <- readRDS(file = BPfile)
  
  # Restrict to a specific window, or not.
  if(window == "earlyP") BP <- BP %>% dplyr::filter(GA<=limit)
  if(window == "lateP") BP <- BP %>% dplyr::filter(GA>limit)
  
  BP <- BP %>% dplyr::group_by(id, GA) %>% 
               dplyr::mutate(SBP = mean(SBP),
                             DBP = mean(DBP),
                             MAP = mean(MAP),
                             PP = mean(PP)) %>%
              dplyr::distinct(id, GA, .keep_all = T) %>%
              data.frame()
  dim(BP)
  
  ##############################################################################
  ## MLM - Stage one: estimating individual mean levels 
  ## (autocorrelation structure - ...why not?)
  #----------------------------------------------------------------------------#
  
  BPhat <- MLM2_stage1(exposure = "SBP", time = time, data = BP)
  
  BPhat <- merge(BPhat,
                 MLM2_stage1(exposure = "DBP", time = time, data = BP),
                 by = "id")
  
  BPhat <- merge(BPhat,
                 MLM2_stage1(exposure = "MAP", time = time, data = BP),
                 by = "id")
  
  BPhat <- merge(BPhat,
                 MLM2_stage1(exposure = "PP", time = time, data = BP),
                 by = "id")
  
  dim(BPhat) %>% print()
  
  ################################################################################
  # Merge with covariates
  #------------------------------------------------------------------------------#
  
  BPhat <- merge(BPhat, readRDS(file = covFile), by = "id", all = F)
  
  CC <- readRDS(file = CCFile)
  
  BPhat <- merge(BPhat, data.frame(id = row.names(CC), CC), all = F)
  
  dim(BPhat) %>% print()
  
  ################################################################################
  # Save
  #------------------------------------------------------------------------------#
  
  saveRDS(BPhat, file = outFile)
  
}

################################################################################
rm(list = ls())
################################################################################
#!/usr/bin/env Rscript
cat("[BP-DNAme] Running Rscript 1-BP-Centering...\n")
################################################################################
# BP-DNAme -  Detrending of Blood pressure values
################################################################################
# Author: Lucile
# Date: August 2021
# Update: 18th November 2021
# Removal of subjects with known hypertensive disorder prior to mean computations
#
# Note: 
# - 0-BP-Preprocessing before
# - 0-Covariates-Preprocessing.R
################################################################################
# Description:
# *BP* trajectories during pregnancy are represented as a non-stationary time serie:
# \[ BP_{t} = \mu_{t} + y_{t} \]
# where $\mu_{t}$ is the non-linear trend 
# and $y_{t} = BP_{t} - \mu_{t}$ is a stationary time series $BP$
#------------------------------------------------------------------------------#
# Method:
# Starting from raw data extracted using Rscript 0-PreprocessBP.R
# 1. Remove time-trend, by centering BP measurements according to the EDEN cohort 
# (robust) average *BP level, computed over a month-wide window
# 2. Estimate individual average *BP level (over a given time period) 
# using a mixed linear model assuming no time trend: BP ~ 1 + (1|id) 
# accounting/not accounting for time autocorrelation (AR1)
#------------------------------------------------------------------------------#
# Note:
# Estimations from the whole EDEN cohort (N=668+1334)
################################################################################
# Paths and parameters
#------------------------------------------------------------------------------#

# Input file, generated using Rscript 0-BP-Preprocessing 
BPfile <- "Data/N2002_BP.rds"

# Covariates, generated using Rscript 0-Covariates-Preprocessing.R
covFile <- "Data/N2002_covariates.rds"

# Ouput file (centered individual BP values)
outFile <-  "Data/N2002_BP_centered.rds"

alpha <- 0.05 # Trimmed mean with alpha% most extreme values excluded (both sides)
window_width <- 2.5

################################################################################
# R packages
#------------------------------------------------------------------------------#

library(magrittr)

################################################################################
# R functions
#------------------------------------------------------------------------------#

movingAverage <- function(timepoint, timepoints, values, alpha, window_width){
  
  indices <- which(timepoints>=timepoint-window_width &
                     timepoints<=timepoint+window_width)
  mean <- ifelse(length(indices>0), mean(values[indices], trim = alpha), NA)
  
  return(mean)
}

################################################################################
#------------------------------------------------------------------------------#
################################################################################
# Blood pressure
################################################################################

BP <- readRDS(file = BPfile)
BP <- BP[complete.cases(BP),] %>% dplyr::arrange(id, GA)
BP$month <- tapply(BP$GA, 1:nrow(BP), function(ga) ceiling(ga/5))

dim(BP); BP.filtered <- merge(BP, readRDS(covFile), by = "id"); dim(BP.filtered)
BP.filtered <- BP.filtered %>% dplyr::filter(h_HTA == 0 & HDP == 0)

################################################################################
# Recentrage par mois
# Moyennes calculées sur la pop. pour chaque mois de grossesse
################################################################################

#SBP
BP$SBP.meansw <- tapply(BP$GA, 1:nrow(BP),
                        function(x) movingAverage(timepoint = x, 
                                                  timepoints = BP.filtered$GA, 
                                                  values = BP.filtered$SBP,
                                                  alpha = alpha, 
                                                  window_width = 2.5))

#DBP
BP$DBP.meansw <- tapply(BP$GA, 1:nrow(BP),
                        function(x) movingAverage(timepoint = x, 
                                                  timepoints = BP.filtered$GA, 
                                                  values = BP.filtered$DBP,
                                                  alpha = alpha, 
                                                  window_width = 2.5))


#MAP
BP$MAP.meansw <- tapply(BP$GA, 1:nrow(BP),
                        function(x) movingAverage(timepoint = x, 
                                                  timepoints = BP.filtered$GA, 
                                                  values = BP.filtered$MAP,
                                                  alpha = alpha, 
                                                  window_width = 2.5))

#PP
BP$PP.meansw <- tapply(BP$GA, 1:nrow(BP),
                       function(x) movingAverage(timepoint = x, 
                                                 timepoints = BP.filtered$GA, 
                                                 values = BP.filtered$PP,
                                                 alpha = alpha, 
                                                 window_width = 2.5))

BP$SBP <- BP$SBP-BP$SBP.meansw
BP$DBP <- BP$DBP-BP$DBP.meansw
BP$MAP <- BP$MAP-BP$MAP.meansw
BP$PP <- BP$PP-BP$PP.meansw

################################################################################
cat("[BP-DNAme] Centered BP values saved in file:", outFile, "\n")
#------------------------------------------------------------------------------#

saveRDS(BP, file = outFile)

################################################################################
rm(list = ls())
################################################################################
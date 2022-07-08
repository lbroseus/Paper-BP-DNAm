#!/usr/bin/env Rscript
rm(list = ls())
################################################################################
# EWAS: BP ~ methylation 
################################################################################
# Author: Lucile
# Date of Creation: April 2021
# Update: June 2021 -> have it run on the server
# Update: Novembre 2021
# Update: January 2021 
################################################################################
# Notes:
# EDEN data (N = 668 - 2)
# Parameters to be adapted depending on the analysis
################################################################################

library(magrittr)

################################################################################
# Paths and parameters
#------------------------------------------------------------------------------#

exposure <- "MAP"
cat("Tested exposure:", exposure)

# for BP indicators (SBP, DBP, MAP, PP)
clinical_confounders <- c("center", "age", "parity", "sex", "BMI", "education", "smoke")

# Type of adjustment (whether to adjust on CC and LF)
adjustmentType <- "I" # must be in c("I", "II", "III")
cat("Adjustment type:", adjustmentType)

# Window of exposure
timePeriod <- "earlyP" # must be in c("pregnancy", "earlyP", "lateP")
cat("Time period:", timePeriod)

# Modelling scheme
method <- "mlm2"

################################################################################
# Sensitivity analyses (one at a time)
#------------------------------------------------------------------------------#

# Adjustment on gestational duration
adjOnGestDuration <- FALSE
# Pre-eclampsia (PE)
includePE <- TRUE
# Prematurity (preterm)
includePreTerm <- TRUE
# Gestational diabetes (GDM)
includeGDM <- TRUE
# Smoking status (smoke)
includeSmokers <- TRUE

#------------------------------------------------------------------------------#

dataPath <- "Data"

resultPath <- "Results"

#------------------------------------------------------------------------------#
# Shared parameters (assumed to be fixed for all analyses..?)

rm_samples <- c("7327")

methFile <- paste0(dataPath, "/bmiq_processed_noOutliers.rds")
covFile <- paste0(dataPath, "/metaData.rds")
BPFile <- paste0(dataPath, "/N668_BP1+covariates_", timePeriod, ".rds") #Depends on the window of exposure
stopifnot( file.exists(BPFile) )

cellComp <- paste("CC", 1:5, sep = "")
hidden_confounders <- paste("LF", 1:5, sep = "")
known_technical_confounders <- c("batch", "chip", "plate")

transformToMvalue <- TRUE
maxit <- 400  #For robust regressions
ncores <- min(2,bigstatsr::nb_cores())

################################################################################
# R packages

#devtools::load_all(paste0(customPackagesPath, "/LuBReg/"))

if( !require(LuBReg) ) devtools::install_github("lbroseus/LuBReg")
library(LuBReg)

################################################################################
# R functions

################################################################################
# Data management
#------------------------------------------------------------------------------#

#--------------------#
# Methylation data
#--------------------#

Bvalues <- readRDS(file = methFile)

#--------------------#
# BP and exposures
#--------------------#

BP <- readRDS(BPFile)

length(unique(BP$id))

################################################################################
# Sensitivity analysis (?)
#------------------------------------------------------------------------------#

crit <- sum(as.numeric(c(!includeGDM, !includePE, !includeSmokers, !includePreTerm, adjOnGestDuration)))

if( crit > 1 ){
  stop("Please provide only one filter for sensitivity analysis \n")
}else if(crit == 0 ){
  subDir <- "Main"
  #All individuals will be included
}else if( !includePE ){
  cat("Sensitivity analysis on PE cases (excluded). \n")
  subDir <- "PE"
  BP <- BP %>% dplyr::filter(HDP != 2)
}else if( !includeGDM ){
  cat("Sensitivity analysis on GDM cases (excluded). \n")
  subDir <- "GDM"
  BP <- BP %>% dplyr::filter(GDM != 1)
}else if( !includePreTerm ){
  cat("Sensitivity analysis on pre-term births \n")
  subDir <- "PreTerm"
  BP <- BP %>% dplyr::filter(gestDuration > 37)
}else if( !includeSmokers ){
  cat("Sensitivity analysis on smokers (excluded) \n")
  subDir <- "Smoke"
  BP <- BP %>% dplyr::filter(smoke != "continuedSmoker")
}else{
  subDir <- "gestDuration"
  clinical_confounders <- c(clinical_confounders, "gestDuration")
  cat("Sensitivity analysis: additional adjustment on gestational duration.\n")
}

################################################################################
# Categorical covariates as factor (here, only needed for technical confounders)
#------------------------------------------------------------------------------#

BP[,known_technical_confounders] <- apply(BP[,known_technical_confounders], 2, as.factor)

################################################################################
# Set technical confounders to adjust for
#------------------------------------------------------------------------------#

technical_confounders <- known_technical_confounders
if(adjustmentType %in% c("II", "III")) technical_confounders <- c(technical_confounders, cellComp)
if(adjustmentType %in% c("III")) technical_confounders <- c(technical_confounders, hidden_confounders)
cat("Technical confounders:", technical_confounders, "\n")

################################################################################
# Set output files and directories
#------------------------------------------------------------------------------#

file_name <- paste0(exposure, ".", method, ".", timePeriod)
print(file_name)

path <- paste0(resultPath, "/")
if( !dir.exists(path)) dir.create(path)  

path <- paste0(path, "/", paste("EWAS",adjustmentType, sep = ""), sep = "")
if( !dir.exists(path)) dir.create(path)

path <- paste0(path, "/", subDir)
if( !dir.exists(path)) dir.create(path)

print(path)

################################################################################
# Synchronize data sets
#------------------------------------------------------------------------------#

mmList <- list(Bvalues, BP)

#----------------------------------------#
# Processing of mmList[["methylData"]]
#----------------------------------------#

names(mmList) <- c("methylData", "metaData")

# rm specified samples
if(!is.null(rm_samples) ) mmList[["methylData"]] <- mmList[["methylData"]][, !(colnames(mmList[["methylData"]]) %in% rm_samples)]
#rm samples not in the meta data set
mmList[["methylData"]] <- mmList[["methylData"]][,colnames(mmList[["methylData"]]) %in% mmList[["metaData"]]$id]
#order samples by name
mmList[["methylData"]] <- mmList[["methylData"]][, order(colnames(mmList[["methylData"]]))]

cat("Methylation dataset has", 
    ncol(mmList[["methylData"]]), "individuals and spans",
    nrow(mmList[["methylData"]]), "sites.\n")

#----------------------------------------#
# Processing of mmList[["metaData"]]
#----------------------------------------#

mmList[["metaData"]] <- mmList[["metaData"]] %>% 
  dplyr::filter(id %in% colnames(mmList[["methylData"]])) %>%
  dplyr::arrange(id)

################################################################################
# Two-stage linear regressions with random effect
#------------------------------------------------------------------------------#
# methyl ~ DBP/SBP/MAP/PP
#------------------------------------------------------------------------------#

system.time(
  
  results <- LuBReg::LocusWiseRLM(meth_data = mmList[[ 1 ]],
                                  covariates = mmList[[ 2 ]],
                                  exposure,
                                  clinical_confounders,
                                  technical_confounders,
                                  maxit = maxit,
                                  transformToMvalue = transformToMvalue,
                                  ncores = ncores,
                                  path,
                                  file_name)
  
)
################################################################################

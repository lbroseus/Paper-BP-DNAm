#!/usr/bin/env Rscript
cat("[BP-DNAme] Running Rscript 1-BP-fPCA...\n")
################################################################################
# EDEN - Perform function PCA on BP data
################################################################################
# Author: Lucile
# Date of creation: April 2021
################################################################################
# Note:
# Estimations from the whole EDEN cohort (N=668+1334)
################################################################################
# Paths and parameters
#------------------------------------------------------------------------------#

# Input file, generated using Rscript 0-BP-Preprocessing 
BPfile <- "Data/N2002_BP.rds"

# Ouput file (centered individual BP values)
outFile <-  "Data/N2002_BP_fPCA.RData"

# Min number of measurement for including a subject
minNbObs <- 1

# Bounds for time (GA)
timeLow <- 5
timeUp <- 45

################################################################################
# R packages
#------------------------------------------------------------------------------#

library(magrittr)
library(fdapace)

################################################################################
# R functions
#------------------------------------------------------------------------------#

# Function performFPCA()
# Expected column id
# varName: points a numerical variable
# timeVarName: numeric column
# returns fpca object

performFPCA <- function(df, varName, timeVarName, minNbObs = 1){
  
  #--------------------------------------#
  # Check input parameters
  #--------------------------------------#
  
  if( ! "id" %in% colnames(df) ) stop("Missing column 'id'.\n")
  if( ! varName %in% colnames(df) ) stop(paste("Missing column",  varName, ".\n"))
  if( ! timeVarName %in% colnames(df) ) stop(paste("Missing column",  timeVarName, ".\n"))
  
  #--------------------------------------#
  # Prepare input data
  #--------------------------------------#
  
  # Standardize colnames
  colnames(df)[which(colnames(df) == varName)] <- "measurement"
  colnames(df)[which(colnames(df) == timeVarName)] <- "time"
  
  # Remove missing values
  df <- df %>% dplyr::filter( !(is.na( measurement ) ) ) 
  
  # Compute # obs per individual and filter
  df <- df %>% dplyr::group_by(id) %>% dplyr::mutate(nbObs = n()) %>% data.frame()
  df <- df[which(df$nbObs >= minNbObs),]
  
  cat("There are", length(unique(df$id)), "individuals with at least", minNbObs, "measurements \n")
  
  cat("Overall, there are", nrow(df), "data points to be used\n")
  
  # If needed, aggregate simultaneous obs. by their mean
  df <- df %>% dplyr::group_by(id, time) %>% 
    dplyr::summarize(measurement = mean(measurement)) %>%
    data.frame() %>% 
    dplyr::arrange(id, time)
  
  #--------------------------------------#
  # Shape for fdapace functions
  #--------------------------------------#
  
  ## Observations
  yList <- split(df[, c("measurement")],f = factor(df[, c("id")]))
  names(yList) <- levels( factor(df[, c("id")]) )
  length(yList)
  
  ## Time
  tList <- split(df[, c("time")],f = factor(df[, c("id")]))
  names(tList) <- levels( factor(df[, c("id")]) )
  
  #--------------------------------------#
  # Run FPCA
  #--------------------------------------#
  
  fpca <- fdapace::FPCA(Ly = yList, Lt = tList, list(plot = FALSE))
  rownames(fpca$xiEst) <- names(yList)
  
  return( fpca )
  
}

################################################################################
# Perform fPCA for each indicator
#------------------------------------------------------------------------------#

BP <- readRDS(BPfile)

BP <- BP %>% dplyr::filter(GA >= timeLow & GA <= timeUp)
  
fpca_SBP <- performFPCA(df = BP, 
                        varName = "SBP", 
                        timeVarName = "GA", 
                        minNbObs = minNbObs)
  
fpca_DBP <- performFPCA(df = filter(BP, !is.na(DBP)), 
                        varName = "DBP", 
                        timeVarName = "GA", 
                        minNbObs = minNbObs)
  
fpca_MAP <- performFPCA(df = filter(BP, !is.na(MAP)), 
                        varName = "MAP", 
                        timeVarName = "GA", 
                        minNbObs = minNbObs)
  
fpca_PP <- performFPCA(df = filter(BP, !is.na(PP)), 
                       varName = "PP", 
                       timeVarName = "GA", 
                       minNbObs = minNbObs)
  

################################################################################
cat("[BP-DNAme] fPCA results saved in file:", outFile, "\n")
#------------------------------------------------------------------------------#

save(fpca_DBP, fpca_SBP, fpca_MAP, fpca_PP, minNbObs, file = outFile)

################################################################################
rm(list = ls())
################################################################################
rm(list = ls())
################################################################################
# BP-DNAme: GAMP Test
# Association test between global methylation profiles and
# average maternal bloog pressure indicators using GAMP
################################################################################
# Author: Lucile
# Date of creation: November 2021 
# Update: June 2021, added GAMP analysis with adjutement for cell comp.
################################################################################
# Note:
################################################################################
# Paths and parameters
#------------------------------------------------------------------------------#

#Meth data generated using Rscript 0_DNAme-Preprocessing.R
methFile.outgone <- "Data/bmiq_processed_noOutliers.rds" 

windows <- c("pregnancy", "earlyP", "lateP")
exposures <- c("SBP", "DBP", "MAP", "PP")
# Clinical Covariates for adjustment (0-Covariates_Preprocessing.R)
phenoFilePrefix <- "~/Work/BP-DNAme/Data/N668_BP1+covariates_"

rmSamples <- c("7327")

################################################################################
# R packages
#------------------------------------------------------------------------------#

library(magrittr)
#library(data.table)
library(GAMP)

################################################################################
# Load methylation data
#------------------------------------------------------------------------------#

meth <- readRDS(methFile.outgone) 
meth <- meth[, order(colnames(meth))]

################################################################################
# GAMP analysis with main model
#------------------------------------------------------------------------------#

clinical_confounders <- c("center", "age", "sex", "BMI", "parity", "education", "smoke")
technical_confounders <- c("batch", "chip", "plate")
cellTypes <- 0  # set cellTypes <- 0 for not adjusting

outFile <- "~/Work/BP-DNAme/Results/Table/BP-GAMP.xlsx"

#------------------------------------------------------------------------------#
# Prepare result data.frame
#------------------------------------------------------------------------------#

analyses <- data.frame(window = rep(windows, length(exposures))) %>% dplyr::arrange(window)
analyses <- cbind.data.frame(analyses,
                             data.frame(exposure = rep(exposures, length(windows))))
analyses$main_exposition <- paste(analyses$exposure, analyses$window, sep = "_")

analyses$pval.densities <- NA
analyses$pval.cdf <- NA

analyses

#------------------------------------------------------------------------------#
# Apply GAMP
#------------------------------------------------------------------------------#

for( e in seq_along(analyses$main_exposition)){ 
  
  exposure <- analyses$exposure[e]
  phenoFile <- paste0(phenoFilePrefix, analyses$window[e], ".rds")
  
  PHENO <- readRDS( phenoFile )
  PHENO[, technical_confounders] <- apply(PHENO[, technical_confounders], 2, as.factor)
  
  PHENO <- PHENO %>% dplyr::mutate(batch = ifelse(batch>0, 1, 0),
                                   bmi = ifelse(BMI %in% c("underweight", "normal"), 0, 1),
                                   smoke = ifelse(smoke == "continuedSmoker", 1, 0))
  
  ################################################################################
  # rm mixed-up samples
  #------------------------------------------------------------------------------#
  
  PHENO <- PHENO[!(as.character(PHENO$id) %in% rmSamples),]
  PHENO <- PHENO[order(as.character(PHENO$id)),]
  dim(PHENO)
  
  Z <- as.matrix(meth[,colnames(meth) %in% PHENO$id])
  y <- as.vector(PHENO[,exposure])
  
  X <-  as.matrix(PHENO[,c(clinical_confounders, technical_confounders)])

  analyses$pval.densities[e] <- GAMP::TestDensities(Z = Z,
                                                    y = y,
                                                    X = X, hideProgress = T)
  
  analyses$pval.cdf[e] <- GAMP::TestCDFs(Z = Z,
                                         y = y,
                                         X = X, hideProgress = T)
  
}

analyses

################################################################################
# Save - main model
#------------------------------------------------------------------------------#

xlsx::write.xlsx(analyses, file = outFile)

################################################################################
# GAMP analysis adjusting for CC
#------------------------------------------------------------------------------#

clinical_confounders <- c("center", "age", "sex", "BMI", "parity", "education", "smoke")
technical_confounders <- c("batch", "chip", "plate")
cellTypes <- paste0("CC", 1:5)  # set cellTypes <- 0 for not adjusting

outFile <- "Results/Tables/BP-GAMP-CCadj.xlsx"

#------------------------------------------------------------------------------#
# Prepare result data.frame
#------------------------------------------------------------------------------#

analyses <- data.frame(window = rep(windows, length(exposures))) %>% dplyr::arrange(window)
analyses <- cbind.data.frame(analyses,
                             data.frame(exposure = rep(exposures, length(windows))))
analyses$main_exposition <- paste(analyses$exposure, analyses$window, sep = "_")

analyses$pval.densities <- NA
analyses$pval.cdf <- NA

analyses

#------------------------------------------------------------------------------#
# Apply GAMP
#------------------------------------------------------------------------------#

for( e in seq_along(analyses$main_exposition)){ 
  
  exposure <- analyses$exposure[e]
  phenoFile <- paste0(phenoFilePrefix, analyses$window[e], ".rds")
  
  PHENO <- readRDS( phenoFile )
  PHENO[, technical_confounders] <- apply(PHENO[, technical_confounders], 2, as.factor)
  
  PHENO <- PHENO %>% dplyr::mutate(batch = ifelse(batch>0, 1, 0),
                                   bmi = ifelse(BMI %in% c("underweight", "normal"), 0, 1),
                                   smoke = ifelse(smoke == "continuedSmoker", 1, 0))
  
  ################################################################################
  # rm mixed-up samples
  #------------------------------------------------------------------------------#
  
  PHENO <- PHENO[!(as.character(PHENO$id) %in% rmSamples),]
  PHENO <- PHENO[order(as.character(PHENO$id)),]
  dim(PHENO)
  
  Z <- as.matrix(meth[,colnames(meth) %in% PHENO$id])
  y <- as.vector(PHENO[,exposure])
  
  X <-  as.matrix(PHENO[,c(clinical_confounders, technical_confounders, cellTypes)])
  
  analyses$pval.densities[e] <- GAMP::TestDensities(Z = Z,
                                                    y = y,
                                                    X = X, hideProgress = T)
  
  analyses$pval.cdf[e] <- GAMP::TestCDFs(Z = Z,
                                         y = y,
                                         X = X, hideProgress = T)
  
}

analyses

################################################################################
# Save - Model CC adjusted
#------------------------------------------------------------------------------#

xlsx::write.xlsx(analyses, file = outFile)

################################################################################

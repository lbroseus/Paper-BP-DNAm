#!/usr/bin/env Rscript
cat("[BP-DNAme] Running Rscript 4-TestGlobalMeth.R...\n")
################################################################################
# BP-DNAme: Test association with several summaries of global methylation
#           (LINE, Alu, median meth)
# Association test between global methylation levels (LINE/Alu) and
# average maternal bloog pressure indicators
################################################################################
# Author: Lucile
# Date of creation: November 2021 
################################################################################
# Notes:
################################################################################
# Paths and parameters
#------------------------------------------------------------------------------#

#Design information:
methglobFile <- "~/Data/methglobsd2017_04_06.RData" # object methglobsd

#Meth data generated using Rscript 0-EDEN-DNAmePreprocessing.R
methFile.outgone <- "~/Work/bmiq_processed_noOutliers.rds" 

clinical_confounders <- c("center", "age", "sex", "BMI", "parity", "education", "smoke")
technical_confounders <- c("batch", "puce", "plaque")
cellTypes <- 0  # set cellTypes <- 0 for not adjusting

windows <- c("pregnancy", "earlyP", "lateP")
exposures <- c("SBP", "DBP", "MAP", "PP")
# Clinical Covariates for adjustment (0-Covariates_Preprocessing.R)
phenoFilePrefix <- "~/Work/BP-DNAme/Data/N668_BP1+covariates_"
# Output file
outFile <- "Results/Tables/BP-globmeth.xlsx"

rmSamples <- c("7327")

################################################################################
# R packages
#------------------------------------------------------------------------------#

library(magrittr)
#library(data.table)

################################################################################
# Fonctions source (copied from file)
#------------------------------------------------------------------------------#

source("Rscripts/Rfunctions.R")

################################################################################
# Prepare result data.frame
#------------------------------------------------------------------------------#

analyses <- data.frame(window = rep(windows, length(exposures))) %>% dplyr::arrange(window)
analyses <- cbind.data.frame(analyses,
                             data.frame(exposure = rep(exposures, length(windows))))
analyses$main_exposition <- paste(analyses$exposure, analyses$window, sep = "_")

analyses$Estimate.line <- analyses$pval.line <- NA
analyses$Estimate.alu <- analyses$pval.alu <- NA
analyses$Estimate.median <- analyses$pval.median <- NA

analyses

################################################################################
# Load needed columns from the PHENO dataframe 
# (renamed covariates + cell components)
#------------------------------------------------------------------------------#

setwd(workDir)

################################################################################
# Load global methylation data
#------------------------------------------------------------------------------#

load( methglobFile )

################################################################################
# LINE-1 / Alu / Global median
#------------------------------------------------------------------------------#

cor <- cor(methglobsd$alumed, methglobsd$linemed, use = "pairwise.complete.obs", method = "spearman")
plot(x = methglobsd$alumed, xlab = "Median Alu methylation", 
     y = methglobsd$linemed, ylab = "Median LINE-1 methylation", 
     main = paste("Spearman correlation = ", round(cor,2)),
     pch = 19, cex = 0.75)

################################################################################
# Load and summarize methylation data
#------------------------------------------------------------------------------#

meth <- readRDS(methFile.outgone) 

meth <- data.frame(id = colnames(meth),
                   globmed = apply(meth, 2, function(x) median(x, na.rm = T)),
                   globiqr = apply(meth, 2, function(x) iqr(x, na.rm = T)))
dim(meth)
# 668 3

for( e in seq_along(analyses$main_exposition)){ 
  
  phenoFile <- paste0(phenoFilePrefix, analyses$window[e], ".rds")
  
  PHENO <- readRDS( phenoFile )
  PHENO <- merge(x = PHENO[,!colnames(PHENO)%in%technical_confounders], y = methglobsd, by = "id")
  PHENO[, technical_confounders] <- apply(PHENO[, technical_confounders], 2, as.factor)
  
  ################################################################################
  # rm mixed-up samples
  #------------------------------------------------------------------------------#
  
  PHENO <- PHENO[!(as.character(PHENO$id) %in% rmSamples),]
  PHENO <- PHENO[order(as.character(PHENO$id)),]
  dim(PHENO)
  
  ################################################################################
  # Append globam methylation statistics
  #------------------------------------------------------------------------------#
  
  PHENO <- merge(PHENO, meth, by = "id", all.x = T, all.y = F)
  dim(PHENO)
  
  ################################################################################
  # Logit transform median proportions
  #------------------------------------------------------------------------------#
  
  PHENO$linemed.logit <- tapply(PHENO$linemed, 
                                seq_along(PHENO$linemed),
                                function(p) ifelse(is.na(p), NA, log(p/(100-p))))
  PHENO$alumed.logit <- tapply(PHENO$alumed, 
                               seq_along(PHENO$alumed),
                               function(p) ifelse(is.na(p), NA, log(p/(100-p))))
  
  PHENO$globmed.logit <- tapply(PHENO$globmed, 
                                seq_along(PHENO$globmed),
                                function(p) ifelse(is.na(p), NA, log(p/(100-p))))

  
  ################################################################################
  # Association tests
  #------------------------------------------------------------------------------#
  
  exposure <- analyses$exposure[e]
  
  #------------------------------------------------------------------------------#
  # LINE1
  #------------------------------------------------------------------------------#
    
  formula <- stats::as.formula(paste("linemed.logit ~", 
                                       exposure, "+",
                                       paste(c(clinical_confounders, technical_confounders, cellTypes), collapse = "+")))

  #print( formula )
    
  mod <- rlm(formula = formula, data = PHENO)
  mod <- coeftest(mod, vcov = vcovHC(mod, type = "HC0"))  
    
  analyses$Estimate.line[e] <- mod[exposure, "Estimate"]
  analyses$pval.line[e] <- mod[exposure, "Pr(>|z|)"]
    
  #------------------------------------------------------------------------------#
  # ALU
  #------------------------------------------------------------------------------#
    
  formula <- stats::as.formula(paste("alumed.logit ~", 
                                       exposure, "+",
                                       paste(c(clinical_confounders, technical_confounders, cellTypes), collapse = "+")))

  #print( formula )
    
  mod <- rlm(formula = formula, data = PHENO)
  mod <- coeftest(mod, vcov = vcovHC(mod, type = "HC0"))  
    
  analyses$Estimate.alu[e] <- mod[exposure, "Estimate"]
  analyses$pval.alu[e] <- mod[exposure, "Pr(>|z|)"]
    
  #------------------------------------------------------------------------------#
  # Global Median
  #------------------------------------------------------------------------------#
    
  formula <- stats::as.formula(paste("globmed.logit ~", 
                                       exposure, "+",
                                       paste(c(clinical_confounders, technical_confounders, cellTypes), collapse = "+")
    ))
  
    
  mod <- rlm(formula = formula, data = PHENO) 
  mod <- coeftest(mod, vcov = vcovHC(mod, type = "HC0"))  
    
  analyses$Estimate.median[e] <- mod[exposure, "Estimate"]
  analyses$pval.median[e] <- mod[exposure, "Pr(>|z|)"]
    
}

analyses

################################################################################
# Save
#------------------------------------------------------------------------------#

write.xlsx(analyses, file = outFile)

################################################################################

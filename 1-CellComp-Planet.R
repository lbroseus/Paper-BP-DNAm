#!/usr/bin/env Rscript
################################################################################
# EDEN - Estimation of reference cell type composition using planet
################################################################################
# Author: Lucile
# Date of creation: April 2021
################################################################################
# Procedure:
# 1. Import Yuan's reference as in Planet
# 2. Filter probes not in the processed methylation dataset
# 3. Use planet and the RPC method to estimate CC
# 4. Impute zeros
# 5. Apply ilr transformation
#------------------------------------------------------------------------------#
# Notes:
################################################################################
# Paths and parameters
#------------------------------------------------------------------------------#

methDataFile <- "Data/bmiq_processed_noOutliers.rds"

saveFile.prop <- "Data/EDEN_CC.planet_rpc.rds"
saveFile.ilr <- "Data/EDEN_CC.planet_rpc_ilr.rds"

################################################################################
# R packages
#------------------------------------------------------------------------------#

library(magrittr)
library(planet)
library(EpiDISH)

################################################################################
#------------------------------------------------------------------------------#
################################################################################
# 1. Processing Yuan's reference probes
#------------------------------------------------------------------------------#

data("plCellCpGsThird")

plCellCpGsThird <- plCellCpGsThird[order(rownames(plCellCpGsThird)),]
plCellCpGsThird <- plCellCpGsThird[,order(colnames(plCellCpGsThird))]

# Number of reference probes:
nbCpG <- nrow(plCellCpGsThird)
nbCpG

################################################################################
# 2. Restrict to probes available in EDEN (after QC filtering)
#------------------------------------------------------------------------------#

betaValues <- readRDS( methDataFile )
betaValues <- betaValues[order(rownames(betaValues)),]

nbCpG <- length(intersect(rownames(betaValues), rownames(plCellCpGsThird)))
nbCpG
# 506 remaining probes

################################################################################
# 3. Estimate CC using the RPC projection method
#------------------------------------------------------------------------------#

EDEN_CC.rpc <- epidish(beta.m = betaValues, 
                       ref.m = plCellCpGsThird, 
                       method = "RPC") 

EDEN_CC.rpc <- EDEN_CC.rpc$estF

# There are zero values (assumed below the LOD):
apply(EDEN_CC.rpc,2,summary)

################################################################################
# 4. Impute zero values in compositions
#------------------------------------------------------------------------------#

X <- EDEN_CC.rpc
X <- apply(X, 1:2, function(x) ifelse(x==0, NA, x))
X <- robCompositions::impCoda(x = X)$xImp
EDEN_CC.rpc <- X
rownames(EDEN_CC.rpc) <- colnames(betaValues)

apply(EDEN_CC.rpc,2,summary)

rm(X)

################################################################################
# Ilr transformation
#------------------------------------------------------------------------------#

EDEN_CC.ilr <- compositions::ilr(EDEN_CC.rpc)
EDEN_CC.ilr <- apply(EDEN_CC.ilr, 2, as.numeric)
rownames(EDEN_CC.ilr) <- colnames(betaValues)
colnames(EDEN_CC.ilr) <- paste0("CC", 1:ncol(EDEN_CC.ilr))

################################################################################
# Save reference and estimates
#------------------------------------------------------------------------------#

saveRDS(file = saveFile.prop, EDEN_CC.rpc)
saveRDS(file = saveFile.ilr, EDEN_CC.ilr)

################################################################################
rm(list = ls())
################################################################################
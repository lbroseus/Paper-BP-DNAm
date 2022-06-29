#!/usr/bin/env Rscript
cat("[BP-DNAme] Running Rscript 0-DNAme-Preprocessing...\n")
################################################################################
# EDEN - Preprocessing methylation data 
################################################################################
# Author: Lucile
# Date: Finalized in April 2021
################################################################################
# Procedure:
# 1. BMIQ normalisation
# 2. Removal of CpGs close to known SNP and Chr XY
# 3. Removal of cross-reactive CpGs 
# 4. [Removal of outlying probes (cf: Abraham2018, Jedynak2021)]
# 5. Create summary flow chart
#------------------------------------------------------------------------------#
# Notes:
################################################################################
# Paths and parameters
#------------------------------------------------------------------------------#

# BMIQ data provided by Johanna
methFile <- "~/Data/EDEN/METHYL/FDF/BMIQ/bmiq_no_duplicates.RDS"

# Final BMIQ data sets
outFile1 <- "Data/bmiq_processed.rds"             # needed for mediation analysis
outFile2 <- "Data/bmiq_processed_noOutliers.rds"  # used for EWAS analyses

################################################################################
# R packages
#------------------------------------------------------------------------------#

library(Gmisc, quietly = TRUE)
library(glue)
library(htmlTable)
library(grid)
library(magrittr)

################################################################################
# R functions
#------------------------------------------------------------------------------#

# Exact same as the function used for paper (Abraham2018)
removeOutliers <- function(probes){
  
  require(matrixStats)
  
  if(nrow(probes) < ncol(probes)) warning("expecting probes are rows (long dataset)")
  
  rowIQR <- rowIQRs(probes, na.rm = T)
  row2575 <- rowQuantiles(probes, probs = c(0.25, 0.75), na.rm = T)
  maskL <- probes < row2575[,1] - 3 * rowIQR 
  maskU <- probes > row2575[,2] + 3 * rowIQR 
  initial_NAs <- rowSums(is.na(probes))
  probes[maskL] <- NA
  removed_lower <- rowSums(is.na(probes))-initial_NAs
  probes[maskU] <- NA
  removed_upper <- rowSums(is.na(probes))-removed_lower-initial_NAs
  N_for_probe <-rowSums(!is.na(probes))
  Log <- data.frame(initial_NAs,removed_lower,removed_upper,N_for_probe)
  
  return( list(probes, Log) )
}

################################################################################
# 1. BMIQ normalisation (cf: platform)
#------------------------------------------------------------------------------#

meth <- readRDS(methFile)
dim(meth)
#379904

if(nrow(meth) < ncol(meth)) meth <- t( meth )

Nprobes <- c(norm = nrow(meth))

################################################################################
# 2. Removal of CpGs close to SNPs (dist<2 bp) and ChrXY
#------------------------------------------------------------------------------#

meth <- DMRcate::rmSNPandCH(meth, 
                            dist = 2, 
                            mafcut = 0.05,
                            rmcrosshyb = FALSE, 
                            rmXY = TRUE)
dim(meth)
# 370 752 CpGs 

Nprobes <- c(Nprobes, snps = nrow(meth))

################################################################################
# 3. Removal of known cross-reactive probes
#------------------------------------------------------------------------------#

xloci <- maxprobes::xreactive_probes(array_type = "450K")
length(xloci)
# 38941 known cross-reactive probes

x <- which(rownames(meth) %in% xloci)
# 2887 known cross-reactive probes still in the dataset

meth <- meth[-x,]
dim(meth)
# 367 865 CpGs 

Nprobes <- c(Nprobes, crossh = nrow(meth))

################################################################################
cat("[BP-DNAme] Pre-processed beta values saved in file:", outFile1, "\n")
#------------------------------------------------------------------------------#

saveRDS(meth, file = outFile1)

################################################################################
# 4. Removal of outlying probes
#------------------------------------------------------------------------------#

system.time( OutlierResults <- removeOutliers(meth) ) 

meth <- OutlierResults[[1]]
dim(meth)
# 379904    668

Log <- OutlierResults[[2]]
save(Log, file = "Data/EDEN_bmiq_Outlier_log.Rdata") #save log
dim(Log)
rm(Log, OutlierResults)

#Filtrage pour exclusion des CpGs qui ont plus de 25% de femmes manquantes
nbdmcpg <- apply(meth,1, function(x) sum(is.na(x))) #nb missing per CpG
nbdmcpg <- data.frame(nbdmcpg)
exclucpg <- nbdmcpg %>% dplyr::filter(nbdmcpg>=25*668/100) #liste des Cpg exclus car plus de 25% de données manquantes
rownames(exclucpg)

if( length(exclucpg)>0 ) meth <- meth[, which(!(rownames(meth)%in%exclucpg))]
dim(meth)

# 379 904 CpGs 
# 0.36% of the values in the data set were outliers. 
# 0 CpGs were removed as none contained >25% of missing values

Nprobes <- c(Nprobes, outliers = nrow(meth))

################################################################################
cat("[BP-DNAme] Pre-processed beta values after outlier exclusion saved in file:", outFile2, "\n")
#------------------------------------------------------------------------------#

saveRDS(meth, file = outFile2)

################################################################################
# 5. Flow Chart
#------------------------------------------------------------------------------#

meth_experiment <- boxGrob(glue("Genome-wide methylation profiling",
                                "(Illumina 450K methylation array, N={Nb})", 
                                Nb = 668, .sep = "\n"))

qc <- boxGrob(glue("Quality controls","({N} CpGs, N={Nb})",
                     N = Nprobes["norm"], Nb = 667, .sep = "\n"))

norm <- boxGrob(glue("Normalization steps","({N} CpGs, N={Nb})",
                      N = Nprobes["norm"], Nb = 667, .sep = "\n"))

snps <- boxGrob(glue("Filtering of SNPs and ChrXY probes","({N} CpGs, N={Nb})",
                     N = Nprobes["snps"], Nb = 667, .sep = "\n"))

crossh <- boxGrob(glue("Removal of cross-reactive probes","({N} CpGs, N={Nb})",
                     N = Nprobes["crossh"], Nb = 667, .sep = "\n"))

outliers <- boxGrob(glue("Removal of outlying intensities","({N} CpGs, , N={Nb})",
                         N = Nprobes["outliers"], Nb = 667, .sep = "\n"))

#______________________________________________________________________________#
grid.newpage()

vertices <- spreadVertical(meth_experiment = meth_experiment,
                           qc = qc, 
                           norm = norm,
                           snps = snps,
                           crossh = crossh,
                           outliers = outliers,
                           .from = 0.95, .to = 0.05)
vertices

for(i in 1:(length(vertices)-1)){
  connectGrob(vertices[[i]], vertices[[i+1]], type = "vertical", 
             arrow_obj = arrow(length = unit(0.15, "inches"), type = "closed")) %>%
    print()
}
#______________________________________________________________________________#

################################################################################
rm(list = ls())
################################################################################
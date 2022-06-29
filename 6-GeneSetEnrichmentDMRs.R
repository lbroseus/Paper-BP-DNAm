#!/usr/bin/env Rscript
cat("[BP-DNAme] Running Rscript 6-GeneSetEnrichmentDMRs.R...\n")
################################################################################
# BP-DNAme - Test DMRs for Gene Set Enrichment...
################################################################################
# Author: Lucile
# Date of creation: November 2021
# Note: 
# Run Rscript 5-DMRs.R before
################################################################################
# Paths and parameters
#------------------------------------------------------------------------------#

# Where DMR results are to be found
resDir <- "~/Work/BP-DNAme/Results/EWASI/Main/DMR/"
saveDir <- "~/Work/BP-DNAme/Results/EWASI/Main/GSEA/"
# Parameters (GSEA)

minNbProbes <- 3
arraytype <- "450K"
ncores <- 2

################################################################################
# R packages
#------------------------------------------------------------------------------#

library(magrittr)
library(ggplot2)
library(missMethyl)
library(GenomicRanges)

suppressPackageStartupMessages( library( IlluminaHumanMethylation450kanno.ilmn12.hg19 ) ) 
annotation_object <- IlluminaHumanMethylation450kanno.ilmn12.hg19::IlluminaHumanMethylation450kanno.ilmn12.hg19
FullAnnot <- minfi::getAnnotation(annotation_object)

################################################################################
# R function(s)
#------------------------------------------------------------------------------#

################################################################################
# RUN
#------------------------------------------------------------------------------#

if( !dir.exists(saveDir) ) dir.create( saveDir )

setwd( resDir )

inFiles <- list.files(path = resDir, pattern = "combp.csv", recursive = T)
inFiles

expositions <- tapply(inFiles, 
                      seq_along(inFiles),
                      FUN = function(x) paste(stringr::str_split(x, 
                                                                 pattern = "\\.", 
                                                                 simplify = T)[c(1,3)], collapse = "-"))
expositions

# Run GSEA on gene subset (genes associated to DMRs)
fdr.thr <- 0.15
PDE.thr <- 0.05

GSEA <- data.frame()
for(f in seq_along( inFiles )){
  
  print( inFiles[f] )
  res <- read.csv( inFiles[f] )
  res <- res %>% dplyr::filter(nprobe>=minNbProbes & !is.na(geneName))
  
  if( nrow(res>0) ){
  # regions: GRanges object of DMR coordinates to test for GO term enrichment
  regions <- GRanges(seqnames = res$chr, ranges = IRanges(start = res$start,end = res$end))
  
  # KEGG pathway enrichment
  
  # Enrichment 
  df_tmp <- missMethyl::goregion(regions = regions, 
                                 collection = "KEGG",
                                 prior.prob = TRUE,
                                 anno = FullAnnot,
                                 genomic.features = c("ALL"),
                                 array.type = arraytype)
  
  if( nrow(df_tmp) > 0 ) GSEA <- rbind.data.frame(GSEA,
                                                  data.frame(exposure = expositions[f],
                                                             geneSet = "all",
                                                             location = "ALL",
                                                            df_tmp))
  }
}

GSEA %>% dplyr::filter(FDR<0.1)
GSEA %>% dplyr::arrange(P.DE) %>% head(25)

################################################################################
# Save 
#------------------------------------------------------------------------------#

outFile <- paste0(saveDir, "/KEGG_regions.all.rds")
saveRDS(GSEA, file = outFile)

################################################################################
# Go term enrichment
#------------------------------------------------------------------------------#

GSEA <- data.frame()
for(f in seq_along( inFiles )){
  
  print( inFiles[f] )
  res <- read.csv( inFiles[f] )
  res <- res %>% dplyr::filter(nprobe>=minNbProbes & !is.na(geneName))
  
  if( nrow(res>0) ){
    # regions: GRanges object of DMR coordinates to test for GO term enrichment
    regions <- GRanges(seqnames = res$chr, ranges = IRanges(start = res$start,end = res$end))
    
    # KEGG pathway enrichment
    
    # Enrichment 
    df_tmp <- missMethyl::goregion(regions = regions, 
                                   collection = "GO",
                                   prior.prob = TRUE,
                                   anno = FullAnnot,
                                   genomic.features = c("ALL"),
                                   array.type = arraytype) 
    
    if( nrow(df_tmp) > 0 ) GSEA <- rbind.data.frame(GSEA,
                                                    data.frame(exposure = expositions[f],
                                                               geneSet = "all",
                                                               location = "ALL",
                                                               df_tmp))
  }
}

GSEA %>% dplyr::filter(FDR<0.05)
GSEA %>% dplyr::arrange(P.DE) %>% head(25)

################################################################################
# Save 
#------------------------------------------------------------------------------#

outFile <- paste0(saveDir, "/GO_regions.all.rds")
saveRDS(GSEA, file = outFile)

################################################################################
rm(list = ls())
################################################################################
#!/usr/bin/env Rscript
cat("[BP-DNAme] Running Rscript 6-Mediation-ccmm.R...\n")
################################################################################
# BP-DNAme - Mediation analysis: BP -> CellComp -> DNAme
################################################################################
# Author: Lucile
# Date: December 2021
# Note: 
# Run Rscript 5-DMRs.R before
################################################################################
################################################################################
# Paths and parameters
#------------------------------------------------------------------------------#

methFile <- "~/Work/BP-DNAme/Data/bmiq_processed.rds" # no outlier removal (mediate() incompatible with NA in Y)
CCfile <- "~/Work/BP-DNAme/Data/EDEN_CC.planet_rpc.rds"
root <- "~/Work/BP-DNAme/Results/EWASI/"

outFile <- "~/Work/BP-DNAme/Results/EWASI/Mediation/ccmmRefCellTypeMediation.rds"

categorical_confounders <- c("batch","chip", "plate", "smoke", "BMI", "education")
continuous_confounders <- c("age", "parity", "sex", "center")  #or binary

annotationFile <- "~/Work/EDEN/Data/450kanno.ilmn12.hg19.rds"
nprobes.min <- 3

exposures <- c("SBP", "DBP", "MAP", "PP")
windows <- c("pregnancy", "earlyP", "lateP")
#refCellTypes <- c("Stromal", "Syncytiotrophoblast", "Endothelial", "Trophoblasts", "Hofbauer", "nRBC")

# Parameters for ccmm
method.est.cov <- "bootstrap"
n.boot <- 2000
sig.level <- 0.05

################################################################################
# R packages
#------------------------------------------------------------------------------#

library(magrittr)
library(GenomicFeatures)
library(ccmm)

################################################################################
# R functions
#------------------------------------------------------------------------------#

logit2 <- function(p) log2(p/(1-p))

mergeDMRs <- function(fileI.dmr, fileII.dmr, nprobes.min = 3, verbose = TRUE){
  
  dmrsI <- read.csv(fileI.dmr)[,c("chr", "start", "end", "geneName", "nprobe")] %>% 
    dplyr::filter(nprobe>=nprobes.min) %>%
    dplyr::select(-nprobe) %>%
    dplyr::distinct(.keep_all = TRUE)
  dmrsI.gr <- GRanges(seqnames = paste0("chr", dmrsI$chr), 
                      ranges = IRanges(start = dmrsI$start, end = dmrsI$end),
                      strand = "*")
  
  dmrsII <- read.csv(fileII.dmr)[,c("chr", "start", "end", "geneName", "nprobe")] %>% 
    dplyr::filter(nprobe>=nprobes.min) %>%
    dplyr::select(-nprobe) %>%
    dplyr::distinct(.keep_all = TRUE)
  dmrsII.gr <- GRanges(seqnames = paste0("chr", dmrsII$chr), 
                       ranges = IRanges(start = dmrsII$start, end = dmrsII$end),
                       strand = "*")
  
  overlaps <- findOverlaps(query = dmrsI.gr, subject = dmrsII.gr, ignore.strand = TRUE)
  
  if(length(overlaps) == 0){
    
    if(verbose) cat("No overlap found. \n")
    
    dmrs.aggr <- rbind.data.frame(data.frame(dmrsI, sig.in = 'I'),
                                  data.frame(dmrsII, sig.in = 'II'))
    
  }else{
    
    if(verbose) cat("Overlap(s) found! \n")
    
    dmrs.shared <- data.frame(chr = dmrsI$chr[unique(queryHits(overlaps))],
                              start = dmrsI$start[unique(queryHits(overlaps))],
                              end = dmrsI$end[unique(queryHits(overlaps))],
                              geneName = dmrsI$geneName[unique(queryHits(overlaps))],
                              sig.in = "I+II")
    
    if(length(unique(queryHits(overlaps))) < nrow(dmrsI)){
      dmrsI <- data.frame(dmrsI[-unique(queryHits(overlaps)),], sig.in = 'I') 
    }else{
      dmrsI <- data.frame()
    }
    
    if(length(unique(subjectHits(overlaps))) < nrow(dmrsII)){
      dmrsII <- data.frame(dmrsII[-unique(subjectHits(overlaps)),], sig.in = 'II') 
    }else{
      dmrsII <- data.frame()
    }
    
    dmrs.aggr <- rbind.data.frame(dmrsI, dmrsII, dmrs.shared)
  }
  
  return( dmrs.aggr)
  
}

shapeFrame.dmp <- function(file.dmp, 
                           dmrs,
                           suffix = NULL,
                           annotFrame, 
                           nprobes.min = 3, 
                           verbose = TRUE){
  
  dmps <- readRDS(file.dmp)[, c("CpG", "MeanBeta", "Estimate", "SE", "raw_p_value")]
  
  if( !is.null(suffix) ) colnames(dmps)[-1] <- paste0(colnames(dmps)[-1], ".", suffix)
  
  dmps <- merge(dmps, annotFrame, by = "CpG")
  
  dmps <- merge(dmps, dmrs, by = "chr")
  
  # Keep only CpGs included in a DMR
  dmps <- dmps %>% dplyr::filter(pos>=start & pos<=end)
  
  if( verbose ) cat("EWAS:", nrow(dmps), "DMPs are included in a DMR \n")
  
  return( dmps )
  
} 

seekNewAnnot <- function(dmrs.unknown, 
                         genes, 
                         gene.lab = "geneName", verbose = TRUE){
  
  dmrs.unknown <- dplyr::filter(dmrs, geneName == "")
  
  dmrs.gr <- GRanges(seqnames = dmrs.unknown$chr, 
                     ranges = IRanges(start = dmrs.unknown$start, end = dmrs.unknown$end),
                     strand = "*")
  
  seqlevelsStyle(genes) <- seqlevelsStyle( dmrs.gr) <- "UCSC"
  
  overlaps <- findOverlaps(query = dmrs.gr, subject = genes, ignore.strand = TRUE)
  
  if(length(overlaps)>0){
    
    if(verbose) cat("New annotation found!! \n") 
    gene.names <- genes$gene_id[subjectHits(overlaps)]
    gene.names <- mygene::queryMany(gene.names, scopes = "ensembl.gene", fields = "symbol", species = "human")
    dmrs.unknown[queryHits(overlaps),"geneName"] <- ifelse(!is.na(gene.names$symbol), gene.names$symbol, gene.names$query)
    
    dmrs <- rbind.data.frame(dplyr::filter(dmrs, geneName != ""), dmrs.unknown)
    
  }
  
  return( dmrs )
}

getSignComp <- function(IDE.CIs){
  
  x <- which(apply(IDE.CIs, 2, function(v) prod(v)>0))
  if(length(x)>0) return(paste(colnames(IDE.CIs)[x], collapse = ";"))
  else return( "" )
}

################################################################################
# I.1 - Load cell composition data (proportions)
#------------------------------------------------------------------------------#

Comp <- readRDS(CCfile)

################################################################################
# I.2 - Loading annotations
#------------------------------------------------------------------------------#

annot <- readRDS(file = annotationFile)

# To annotate some unknown regions
#library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- makeTxDbFromGFF(file = "~/Data/Annotation/Homo_sapiens.GRCh38.104.gtf.gz")
genes <- genes(txdb)

library(rtracklayer)
path <- system.file(package = "liftOver", "extdata", "hg38ToHg19.over.chain")
chain <- import.chain(path)
seqlevelsStyle(genes) <- "UCSC"
genes <- liftOver(genes, chain) %>% unlist()

################################################################################
# RUN Mediation analysis
#------------------------------------------------------------------------------#

methC <- data.frame()
for(w in seq_along(windows)){
  
  BPfile <- paste0("~/Work/BP-DNAme/Data/N668_BP1+covariates_", windows[w], ".rds")
  
  ##############################################################################
  # 1. Load average BP estimates over the current time window + covariates
  #----------------------------------------------------------------------------#
  
  BP <- readRDS(BPfile)
  BP <- BP %>% dplyr::arrange(id)
  BP <- na.omit(BP)
  
  BP <- merge(BP, data.frame(id = rownames(Comp), Comp), by = "id", all.y = FALSE)
  
  # Dummy variables for BMI, smoke, batch, etc...for ccmm
  X <- fastDummies::dummy_cols(BP[,categorical_confounders], 
                               select_columns = categorical_confounders,
                               remove_first_dummy = TRUE)
  X <- X[, !colnames(X) %in% categorical_confounders]
  X <- cbind(X,BP[, c(continuous_confounders)])
  
  ################################################################################
  # 2. Load beta values
  #------------------------------------------------------------------------------#
  
  Bvalues <- readRDS(methFile)
  Bvalues <- Bvalues[,colnames(Bvalues) %in% BP$id]
  Bvalues <- Bvalues[,order(colnames(Bvalues))]
  
  stopifnot(identical(BP$id, colnames(Bvalues)))
  
  ################################################################################
  # 3. Pre-process covariates
  #------------------------------------------------------------------------------#
  
  BP[, categorical_confounders] <- apply(BP[, categorical_confounders], 2, as.factor)
  
  for(e in seq_along(exposures)){
    
    # DMRs detected in the main analysis (no adjustment on CC)
    resFileI.dmp <- paste0(root, "/Main/DMP/", exposures[e],".mlm2.", windows[w], ".rds")
    resFileI.dmr <- paste0(root, "/Main/DMR/", exposures[e],".mlm2.", windows[w], ".combp.csv")
    # DMRs detected after adjusting for CC
    resFileII.dmp <- paste0(root, "/CellComp/DMP/", exposures[e],".mlm2.", windows[w], ".rds")
    resFileII.dmr <- paste0(root, "/CellComp/DMR/", exposures[e],".mlm2.", windows[w], ".combp.csv")
    
    ############################################################################
    # 3. Selecting CpGs
    #--------------------------------------------------------------------------#
    
    dmrs <- mergeDMRs(fileI.dmr = resFileI.dmr, fileII.dmr = resFileII.dmr)
    
    dmrs <- seekNewAnnot(dmrs, genes = genes)
    
    # Results from EWAS and DMR detection
    dmpsI <- shapeFrame.dmp(file.dmp = resFileI.dmp, 
                            dmrs = dmrs, 
                            suffix = "I",
                            nprobes.min = 3, 
                            annotFrame = annot)
    
    dmpsII <- shapeFrame.dmp(file.dmp = resFileII.dmp, 
                             dmrs = dmrs, 
                             suffix = "II",
                             nprobes.min = 3, 
                             annotFrame = annot)
    
    dmps <- merge(dmpsI, dmpsII, by = intersect(colnames(dmpsI), colnames(dmpsII)), all = T)
    dmps <- merge(dmps, dmrs, by = intersect(colnames(dmps), colnames(dmrs)))
    
    dmps$exposure <- exposures[e]
    dmps$window <- windows[w]
    
    rm(dmpsI, dmpsII)
    
    ############################################################################
    # 4. Compositional mediator (using ccmm)
    #--------------------------------------------------------------------------#
    
    mediation.res <- data.frame()
    for(probe in dmps$CpG){
      
      CpG <- Bvalues[probe,]
      
      ccmm.res <- try( 
        ccmm(y = as.vector(logit2(CpG)), 
             M = as.matrix(BP[,colnames(Comp)]), 
             tr = as.vector(BP[,exposures[e]]), 
             x = as.matrix(X),
             method.est.cov = method.est.cov, 
             n.boot = n.boot, sig.level = sig.level, tol = 1e-04, max.iter = 2000),
        silent = TRUE
      )
      
      if( !class(ccmm.res)[1] == "try-error" ){
        
        colnames(ccmm.res$IDE.CIs) <- colnames(Comp)
        sigCellTypes <- getSignComp(ccmm.res$IDE.CIs)
        
        mediation.res <- rbind.data.frame(mediation.res,
                                           data.frame(CpG = probe,
                                                      DE.effectSize = ccmm.res$DE, 
                                                      DE.CIlow = ccmm.res$DE.CI[1],
                                                      DE.CIup = ccmm.res$DE.CI[2],
                                                      TIDE.effectSize = ccmm.res$TIDE, 
                                                      TIDE.CIlow = ccmm.res$TIDE.CI[1],
                                                      TIDE.CIup = ccmm.res$TIDE.CI[2],
                                                      sigCellTypes = sigCellTypes))
      }
    }
    
    rownames(mediation.res) <- NULL
    mediation.res <- merge(dmps, mediation.res, by = "CpG", all.x = T)
      
    methC <- rbind.data.frame(methC, data.frame(mediation.res, 
                                                exposure = exposures[e], window = windows[w]))  
  }
}


################################################################################
# Save 
################################################################################

saveRDS(methC, file = outFile)

################################################################################
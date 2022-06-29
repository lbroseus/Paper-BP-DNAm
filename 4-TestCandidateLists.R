#!/usr/bin/env Rscript
cat("[BP-DNAme] Running Rscript 4-TestCandidateLists.R...\n")
################################################################################
# BP-DNAme - Candidate lists
################################################################################
# Author: Lucile
# Date of creation: November 2021
# Update: 31st December 2021 (output as multi-sheet excel file)
################################################################################
# 
#------------------------------------------------------------------------------#
# Notes:
# Focus on main analyses : EWAS-I (not adjusted on CC)
################################################################################
# Files/Target
#------------------------------------------------------------------------------#

methFile <- "~/Work/BP-DNAme/Data/bmiq_processed_noOutliers.rds"

# Reports significant hits from Workalemahu2019 and Kazmi2019 studies
candidateFile1 <- "~/Work/Cohorts/EDEN/PJ_methylBP/Data/CandidateList.xlsx"
# Reports significant hits from studies linking DNAme to maternal HDP
candidateFile2 <- "~/Work/Cohorts/EDEN/PJ_methylBP/Data/EWAS_Atlas_association_HDP.txt"

mainOutFile <- "~/Work/BP-DNAme/Results/Tables/SupplementaryTables.xlsx"

################################################################################
# R packages
#------------------------------------------------------------------------------#

require(magrittr)
require(ggplot2)

orderGeneName <- function(name){ 
  name <- stringr::str_split(name, pattern = ";", simplify = T)
  name <- unique(c(name))
  name <- name[stringi::stri_order(name)]
  name <- paste0(name, collapse = ";")
  return( name )
}

################################################################################
# For probe annotation
#------------------------------------------------------------------------------#

suppressPackageStartupMessages( library( IlluminaHumanMethylation450kanno.ilmn12.hg19 ) ) 
annotation_object <- IlluminaHumanMethylation450kanno.ilmn12.hg19::IlluminaHumanMethylation450kanno.ilmn12.hg19

FullAnnot <- minfi::getAnnotation(annotation_object)

FullAnnot <- data.frame(chr = FullAnnot$chr,
                        start = FullAnnot$pos,
                        CpG = FullAnnot$Name,
                        geneName = FullAnnot$UCSC_RefGene_Name,
                        location = FullAnnot$UCSC_RefGene_Group)

FullAnnot$geneName <- tapply(FullAnnot$geneName,
                             INDEX = seq_along(FullAnnot$geneName),
                             FUN = orderGeneName)

################################################################################
################################################################################
# Load candidate list 1 (Workalemahu+ Kazmi) and 2 (HDP from the EWAS Atlas)
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# List from Workalemahu2019 and Kazmi2019
#------------------------------------------------------------------------------#
candidates1 <- xlsx::read.xlsx(candidateFile1, sheetIndex = 1, as.data.frame = T)
candidates1 <- merge(candidates1, FullAnnot, by = "CpG", all.x = T, all.y = F)
head(candidates1)

dplyr::filter(candidates1, Study == "Workalemahu")$CpG %>% unique() %>% length()
dplyr::filter(candidates1, Study == "Kazmi")$CpG %>% unique() %>% length()

# Reshape Phenotype into categories SBP, DBP SBP+DBP
target <- candidates1$CpG[candidates1$Phenotype == "DBP"]
candidates1$Phenotype[which(candidates1$Phenotype=="SBP" & candidates1$CpG%in%target)] <- "SBP+DBP"
target <- candidates1$CpG[candidates1$Phenotype == "SBP+DBP"]
candidates1 <- candidates1[-which(candidates1$Phenotype=="DBP" & candidates1$CpG %in% target),]
table(candidates1$Phenotype)

# CpG Position Gene DirectionOfEffect Phenotype Window
candidates1 <- candidates1 %>% 
               dplyr::filter(!is.na(start)) %>% 
               dplyr::mutate(DirectionOfEffect = ifelse(grepl(pattern = "−", x = Estimate), -1, 1),
                             Position = paste(chr, start, sep = ":")) %>%
               dplyr::select(CpG, Position, Gene = geneName, Source = Study, Phenotype, DirectionOfEffect) %>%
               dplyr::distinct(.keep_all = T)
head(candidates1)

table(candidates1$Source)

# simplify SBP, DBP SBP+DBP

#------------------------------------------------------------------------------#
# List from the EWAS Atlas database
#------------------------------------------------------------------------------#

candidates2 <- data.table::fread(candidateFile2, data.table = F)  %>% 
               dplyr::rename(CpG = probeID)
candidates2 <- merge(candidates2, FullAnnot, by = "CpG", all.x = T, all.y = F)
head(candidates2)

# CpG Position Gene DirectionOfEffect Phenotype
candidates2 <- candidates2 %>% 
  dplyr::mutate(DirectionOfEffect = ifelse(correlation == "pos", 1, -1),
                Position = paste(chr, start, sep = ":"),
                Source = "EWAS Atlas database",
                Phenotype = "HDP") %>%
  dplyr::select(CpG, Position, Gene = geneName, Source, Phenotype,  DirectionOfEffect)
head(candidates2)

#------------------------------------------------------------------------------#
# Combined List 
#------------------------------------------------------------------------------#

# CpG Position Gene Phenotype DirectionOfEffect
candidates.all <- rbind.data.frame(candidates1, candidates2)

################################################################################
# Venn diagram of probe lists
#------------------------------------------------------------------------------#

probeList <- list()
for(sc in unique(candidates.all$Source)){
  probeList[[length(probeList)+1]] <- unique(dplyr::filter(candidates.all, Source == sc)$CpG)
}
names(probeList) <- unique(candidates.all$Source)

str(probeList)

ggvenn::ggvenn(probeList, 
               fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF","#CD534CFF"),
               stroke_size = 0.5, 
               set_name_size = 4)

rm(probeList)

################################################################################
# Merge Kazmi2019 and EWAS database when shared
#------------------------------------------------------------------------------#

target <- candidates.all$CpG[candidates.all$Source == "Kazmi"]
candidates.all$Source[which(candidates.all$Source=="EWAS Atlas database" & candidates.all$CpG %in% target)] <- "EWAS Atlas + Kazmi"
candidates.all <- candidates.all %>% dplyr::filter(Source != "Kazmi")

table(candidates.all$Source)

################################################################################
# Overlap with EWAS results (Main analysis)
#------------------------------------------------------------------------------#

windows <- c("pregnancy", "earlyP", "lateP")
exposures <- c("SBP", "DBP", "MAP", "PP")

dmps.all <- data.frame()
for(w in seq_along(windows)){
  for( e in seq_along(exposures) ){

      resFile <- paste0("~/Work/BP-DNAme/Results/EWASI/Main/DMP/", exposures[e], ".mlm2.", windows[w],".rds")
        
      if( file.exists(resFile) ){  
        dmps <- readRDS(file = resFile) %>% dplyr::filter(CpG %in% candidates.all$CpG) 
        dmps <- dmps %>% dplyr::select(CpG, MeanBeta, pvalue = raw_p_value)
        dmps <- cbind.data.frame(dmps, 
                                 data.frame(window = windows[w], 
                                            exposure = exposures[e]))
        dmps <- merge(candidates.all, dmps, by = "CpG")
        dmps.all <- rbind.data.frame(dmps.all, dmps) 
      }
   }
}
rm(dmps)

################################################################################
# Main tab : summary 
#------------------------------------------------------------------------------#

p.thr <- 0.05

mainTab <- dmps.all %>% 
           dplyr::mutate(DirectionOfEffects = paste0(DirectionOfEffect, "/", sign(MeanBeta))) %>% 
           dplyr::filter(pvalue < p.thr) %>%
           dplyr::group_by(CpG, Position, Gene, Source, Phenotype, DirectionOfEffects, window) %>%
           dplyr::summarize(BPindicators = paste(as.vector(exposure), collapse = ";")) %>%
           data.frame()
head(mainTab)

mainTab <- reshape::cast(mainTab, 
                       formula = CpG+Position+Gene+Source+Phenotype+DirectionOfEffects ~ window, 
                       value = "BPindicators")
head(mainTab)

mainTab[,c("earlyP", "lateP", "pregnancy")] <- apply(mainTab[,c("earlyP", "lateP", "pregnancy")], 
                                                     MARGIN = 1:2,
                                                     FUN = function(x) ifelse(is.na(x), "-", x))

head(mainTab)

################################################################################
# Save 
#------------------------------------------------------------------------------#

outFile <- "~/Work/BP-DNAme/Results/Tables/replicatedCandidateDMPs_main.xlsx"
mainTab %>% dplyr::filter(Source == "Workalemahu") %>% 
            dplyr::arrange(desc(DirectionOfEffects)) %>%
            xlsx::write.xlsx(file = outFile, row.names = F)

################################################################################
# Supplementary Table of significant hits with figures 
#------------------------------------------------------------------------------#

p.thr <- 0.05

#dmps.all <- dmps.all %>% dplyr::filter(pvalue < p.thr & DirectionOfEffect == sign(MeanBeta)) %>%
#            dplyr::arrange(Gene, CpG)

dmps.all <- dmps.all %>% dplyr::filter(pvalue < p.thr) %>%
                         dplyr::mutate(DirectionOfEffects = paste0(DirectionOfEffect, "/", sign(MeanBeta))) %>%
                         dplyr::arrange(CpG)

pvals <- reshape::cast(dmps.all, 
                       formula = CpG+Position+Gene+Source+Phenotype ~ exposure+window, 
                       value = "pvalue")
colnames(pvals)[-c(1:5)] <- paste(colnames(pvals)[-c(1:5)], "Pvalue", sep = ".")

meanbeta <- reshape::cast(dmps.all, 
                          formula = CpG+Position+Gene+Source+Phenotype ~ exposure+window, 
                          value = "MeanBeta")
colnames(meanbeta)[-c(1:5)] <- paste(colnames(meanbeta)[-c(1:5)], "Beta", sep = ".")

DirectionOfEffects <- reshape::cast(dmps.all, 
                          formula = CpG+Position+Gene+Source+Phenotype ~ exposure+window, 
                          value = "DirectionOfEffects")
colnames(DirectionOfEffects)[-c(1:5)] <- paste(colnames(DirectionOfEffects)[-c(1:5)], "DirectionOfEffects", sep = ".")

dmps.all <- merge(pvals, meanbeta, by = intersect(colnames(pvals), colnames(meanbeta)))
dmps.all <- merge(dmps.all, DirectionOfEffects, by = intersect(colnames(dmps.all), colnames(DirectionOfEffects)))
dmps.all<- cbind.data.frame(dmps.all[,c(1:5)], dmps.all[,-c(1:5)][,order(colnames(dmps.all[,-c(1:5)]))])

head(dmps.all)

################################################################################
# Save numeric, supplementary table
#------------------------------------------------------------------------------#

dmps.all <- apply(dmps.all, 1:2, function(x) ifelse(is.na(x), "-", x)) %>% 
            data.frame() %>% 
            dplyr::arrange(dplyr::desc(Source))

xlsx::write.xlsx(dmps.all, 
                 file = mainOutFile, row.names = F,
                 append = T,
                 sheetName = "Candidate List - Main Analysis")

################################################################################
# Sensitivity Analysis: excluding PE cases
################################################################################

windows <- c("pregnancy", "earlyP", "lateP")
exposures <- c("SBP", "DBP", "MAP", "PP")

dmps.all <- data.frame()
for(w in seq_along(windows)){
  for( e in seq_along(exposures) ){
    
    resFile <- paste0("~/Work/BP-DNAme/Results/EWASI/PE//DMP/", exposures[e], ".mlm2.", windows[w],".rds")
    
    if( file.exists(resFile) ){  
      dmps <- readRDS(file = resFile) %>% dplyr::filter(CpG %in% candidates.all$CpG) 
      dmps <- dmps %>% dplyr::select(CpG, MeanBeta, pvalue = raw_p_value)
      dmps <- cbind.data.frame(dmps, 
                               data.frame(window = windows[w], 
                                          exposure = exposures[e]))
      dmps <- merge(candidates.all, dmps, by = "CpG")
      dmps.all <- rbind.data.frame(dmps.all, dmps) 
    }
  }
}
rm(dmps)

################################################################################
# Main tab : summary 
#------------------------------------------------------------------------------#

p.thr <- 0.05

mainTab <- dmps.all %>% 
  dplyr::mutate(DirectionOfEffects = paste0(DirectionOfEffect, "/", sign(MeanBeta))) %>% 
  dplyr::filter(pvalue < p.thr) %>%
  dplyr::group_by(CpG, Position, Gene, Source, Phenotype, DirectionOfEffects, window) %>%
  dplyr::summarize(BPindicators = paste(as.vector(exposure), collapse = ";")) %>%
  data.frame()
head(mainTab)

mainTab <- reshape::cast(mainTab, 
                         formula = CpG+Position+Gene+Source+Phenotype+DirectionOfEffects ~ window, 
                         value = "BPindicators")
head(mainTab)

mainTab[,c("earlyP", "lateP", "pregnancy")] <- apply(mainTab[,c("earlyP", "lateP", "pregnancy")], 
                                                     MARGIN = 1:2,
                                                     FUN = function(x) ifelse(is.na(x), "-", x))

head(mainTab)

################################################################################
# Save 
#------------------------------------------------------------------------------#

outFile <- "~/Work/BP-DNAme/Results/Tables/replicatedCandidateDMPs_main_PE.xlsx"
mainTab %>% dplyr::filter(Source == "Workalemahu") %>% 
  dplyr::arrange(desc(DirectionOfEffects)) %>%
  xlsx::write.xlsx(file = outFile, row.names = F)

################################################################################

################################################################################
# Supplementary Table of significant hits with figures 
#------------------------------------------------------------------------------#

p.thr <- 0.05

dmps.all <- dmps.all %>% dplyr::filter(pvalue < p.thr) %>%
  dplyr::mutate(DirectionOfEffects = paste0(DirectionOfEffect, "/", sign(MeanBeta))) %>%
  dplyr::arrange(CpG)

pvals <- reshape::cast(dmps.all, 
                       formula = CpG+Position+Gene+Source+Phenotype ~ exposure+window, 
                       value = "pvalue")
colnames(pvals)[-c(1:5)] <- paste(colnames(pvals)[-c(1:5)], "Pvalue", sep = ".")

meanbeta <- reshape::cast(dmps.all, 
                          formula = CpG+Position+Gene+Source+Phenotype ~ exposure+window, 
                          value = "MeanBeta")
colnames(meanbeta)[-c(1:5)] <- paste(colnames(meanbeta)[-c(1:5)], "Beta", sep = ".")

DirectionOfEffects <- reshape::cast(dmps.all, 
                                    formula = CpG+Position+Gene+Source+Phenotype ~ exposure+window, 
                                    value = "DirectionOfEffects")
colnames(DirectionOfEffects)[-c(1:5)] <- paste(colnames(DirectionOfEffects)[-c(1:5)], "DirectionOfEffects", sep = ".")

dmps.all <- merge(pvals, meanbeta, by = intersect(colnames(pvals), colnames(meanbeta)))
dmps.all <- merge(dmps.all, DirectionOfEffects, by = intersect(colnames(dmps.all), colnames(DirectionOfEffects)))
dmps.all<- cbind.data.frame(dmps.all[,c(1:5)], dmps.all[,-c(1:5)][,order(colnames(dmps.all[,-c(1:5)]))])

head(dmps.all)

################################################################################
# Save numeric, supplementary table
#------------------------------------------------------------------------------#

dmps.all <- apply(dmps.all, 1:2, function(x) ifelse(is.na(x), "-", x)) %>% 
            data.frame() %>% dplyr::arrange(dplyr::desc(Source))

xlsx::write.xlsx(dmps.all, 
                 file = mainOutFile, row.names = F,
                 append = T,
                 sheetName = "Candidate List - Sensitivity Analysis PE")

################################################################################
rm(dmps.all, meanbeta, pvals)
################################################################################
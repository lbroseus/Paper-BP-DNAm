#!/usr/bin/env Rscript
cat("[BP-DNAme] Running Rscript 5-Association-BP-CellComp.R...\n")
################################################################################
# BP-DNAme - Association test: BP <-> Cell Comp
################################################################################
# Author: Lucile
# Date: February 2022
# Note: 
# A. ANOVA
# B. Whole cell-type composition 
# C. log(Stromal/Syncytiotrophoblasts)
# Run Rscript 4-DMRs.R before
################################################################################
# Paths and parameters
#------------------------------------------------------------------------------#

window <- c("pregnancy", "earlyP", "lateP")
exposure <- c("SBP", "DBP", "MAP", "PP")

technical_confounders <- c("batch", "chip", "plate")

################################################################################
# R packages
#------------------------------------------------------------------------------#

require(compositions)

################################################################################
# Step A "Exploration" of potential general relationships between 
# composition and BP:
# ilr(x1) ~ BP + Confounders
################################################################################

test.res <- data.frame()

for(w in seq_along(window)){
  
  BPhatFile <- paste0("Data/N668_BP1+covariates_", window[w], ".rds")
  
  BPhat <- readRDS( BPhatFile )
  BPhat <- BPhat %>% 
    dplyr::filter(!(id %in% spuriousSamples)) %>%
    dplyr::arrange(id)
  
  BPhat[,technical_confounders] <- apply(BPhat[,technical_confounders],2,as.factor)
  
  Comp <- CCprop[,c("Stromal", "Syncytiotrophoblast", "Endothelial",  'Hofbauer', "nRBC", "Trophoblasts")]
  Comp <- Comp[rownames(Comp) %in% BPhat$id,]
  
  stopifnot(identical(rownames(Comp), BPhat$id))
  
  for(e in seq_along(exposure)){
    
    #Main model, all samples
    formula <- stats::as.formula(paste0("ilr(Comp) ~  ", 
                                        exposure[e], " + center + batch + chip + plate + age + sex + BMI + smoke + parity"))
    
    
    # Tests d'association Comp ~ BP
    res.lm1 <- lm(formula = formula, data = BPhat)
    res.lm1 <- anova(res.lm1)
    
    test.res <- rbind.data.frame(test.res,
                                 data.frame(exposure = exposure[e], 
                                            window = window[w], 
                                            PE = "Main model",
                                            pval.anova = res.lm1[exposure[e], "Pr(>F)"]))
    
    # Model adjusted for GA, all samples
    formula <- stats::as.formula(paste0("ilr(Comp) ~  ", 
                                        exposure[e], " + center + batch + chip + plate + age + sex + BMI + smoke + parity + gestAge"))
    
    
    # Tests d'association Comp ~ BP
    # Adjusting for gestational duration
    res.lm1 <- lm(formula = formula, data = BPhat)
    res.lm1 <- anova(res.lm1)
    
    test.res <- rbind.data.frame(test.res,
                                 data.frame(exposure = exposure[e], 
                                            window = window[w], 
                                            PE = "Adjusted for GA",
                                            pval.anova = res.lm1[exposure[e], "Pr(>F)"]))
    
    # Exclude PE cases but keep main regression model
    BPhat <- BPhat %>% dplyr::filter(HDP != 2)
    
    formula <- stats::as.formula(paste0("ilr(Comp) ~  ", 
                                        exposure[e], " + center + batch + chip + plate + age + sex + BMI + smoke + parity"))
    Comp <- Comp[rownames(Comp) %in% BPhat$id,]
    # Tests d'association
    res.lm1 <- lm(formula = formula, data = BPhat)
    res.lm1 <- anova(res.lm1)
    
    test.res <- rbind.data.frame(test.res,
                                 data.frame(exposure = exposure[e], 
                                            window = window[w], 
                                            PE = "PE excluded",
                                            pval.anova = res.lm1[exposure[e], "Pr(>F)"]))
    
  }
}

################################################################################
cat("Test P-values -Significant results (0.05) \n")
# No significant result
#------------------------------------------------------------------------------#

test.res %>% dplyr::filter(pval.anova<0.05) %>% knitr::kable() %>% print()

################################################################################
# Table of anova p-values
#------------------------------------------------------------------------------#

test.res <- test.res %>% dplyr::filter(PE == "Main model")
test.res <- reshape::cast(test.res, formula = exposure~window, fun.aggregate = mean,value = "pval.anova")

writexl::write_xlsx(test.res,path = "Results/Tables/BPCellType_ANOVA_Revised.xlsx")

################################################################################
rm(test.res)
################################################################################

################################################################################
# Step B:  BP ~ ilr(Comp)
################################################################################

Comp <- CCprop
# Choose reference cell types (as targeted mediator): 
refCellTypes <- colnames(Comp)

test.res <- data.frame()
for(refCellType in refCellTypes){
  
  #cat(refCellType, "as reference cell type :\n")
  for(w in seq_along(window)){
    
    BPhatFile <- paste0("Data/N668_BP1+covariates_", window[w], ".rds")
    
    BPhat <- readRDS( BPhatFile )
    BPhat <- BPhat %>% dplyr::filter(!(id %in% spuriousSamples))
    BPhat[,technical_confounders] <- apply(BPhat[,technical_confounders],2,as.factor)
    
    Comp <- Comp[, c(refCellType, sample(setdiff(colnames(Comp), refCellType), size = ncol(Comp)-1))]
    
    # Apply ilr transformation:
    Comp.ilr <- mixOmics::logratio.transfo(Comp, logratio = 'ILR')
    Comp.ilr <- apply(Comp.ilr, 2, as.numeric)
    rownames(Comp.ilr) <- rownames(Comp)
    colnames(Comp.ilr) <- c(refCellType, paste0("z", 1:(ncol(Comp.ilr)-1)))
    
    for(e in seq_along(exposure)){
      
      formula <- stats::as.formula(paste0(exposure[e], 
                                          " ~ center + batch + chip + plate + age + sex + BMI + smoke + parity +", 
                                          paste(colnames(Comp.ilr), collapse = "+")))
      
      data <- merge(BPhat, data.frame(id = rownames(Comp.ilr), Comp.ilr))
      
      # Tests d'association
      res.lm1 <- lm(formula = formula, data = data)
      
      x <- summary(res.lm1)$coefficients
      y <- confint(res.lm1)[refCellType,]
      z <- data.frame(x[refCellType,'Pr(>|t|)'],
                      EffectSize = x[refCellType,"Estimate"],
                      CI.low = y[1], CI.up = y[2])
      colnames(z)[1] <- "pval.z1"
      
      test.res <- rbind.data.frame(test.res,
                                   data.frame(exposure = exposure[e], 
                                              window = window[w], 
                                              refCellType = refCellType,
                                              PE = "Main model",
                                              z))
      
      formula <- stats::as.formula(paste0(exposure[e], 
                                          " ~ center + batch + chip + plate + age + sex + BMI + smoke + parity + gestDuration + ", 
                                          paste(colnames(Comp.ilr), collapse = "+")))
      
      data <- merge(BPhat, data.frame(id = rownames(Comp.ilr), Comp.ilr))
      
      # Tests d'association
      res.lm1 <- lm(formula = formula, data = data)
      
      x <- summary(res.lm1)$coefficients
      y <- confint(res.lm1)[refCellType,]
      z <- data.frame(x[refCellType,'Pr(>|t|)'],
                      EffectSize = x[refCellType,"Estimate"],
                      CI.low = y[1], CI.up = y[2])
      colnames(z)[1] <- "pval.z1"
      
      test.res <- rbind.data.frame(test.res,
                                   data.frame(exposure = exposure[e], 
                                              window = window[w], 
                                              refCellType = refCellType,
                                              PE = "Adjusted for GA",
                                              z))
      
      # Excluding PE cases
      data <- data %>% dplyr::filter(HDP != 2)
      
      #Back to the main model formula
      formula <- stats::as.formula(paste0(exposure[e], 
                                          " ~ center + batch + chip + plate + age + sex + BMI + smoke + parity +", 
                                          paste(colnames(Comp.ilr), collapse = "+")))
      
      # Tests d'association
      res.lm1 <- lm(formula = formula, data = data)
      x <- summary(res.lm1)$coefficients
      y <- confint(res.lm1)[refCellType,]
      z <- data.frame(x[refCellType,'Pr(>|t|)'],
                      EffectSize = x[refCellType,"Estimate"],
                      CI.low = y[1], CI.up = y[2])
      colnames(z)[1] <- "pval.z1"
      
      test.res <- rbind.data.frame(test.res,
                                   data.frame(exposure = exposure[e], 
                                              window = window[w], 
                                              refCellType = refCellType,
                                              PE = "PE excluded",
                                              z))
      
    }
  }
}

################################################################################
# Effect sizes
#------------------------------------------------------------------------------#

test.res$exposure <- factor(test.res$exposure, levels = c("DBP","SBP","MAP","PP"))
test.res$PE <- factor(test.res$PE, levels = c("Main model", "PE excluded", "Adjusted for GA"))

test.res <- test.res %>% 
  dplyr::mutate(window = ifelse(window == "earlyP", "Early",
                                ifelse(window == "lateP", "Late", "Whole")))

test.res %>% ggplot(aes(x = window, y = EffectSize, col = PE)) +
  geom_pointrange(aes(ymin = CI.low, ymax = CI.up,
  ), size = 0.3, position = position_dodge(width = 0.7)) +
  geom_hline(yintercept = 0) +
  scale_color_manual(values = c("black", "grey63", "darkred"), "Analysis") +
  facet_grid(refCellType ~ exposure) +
  ylab("Effect size") + xlab("Time Window") + 
  theme_minimal()

ggsave("Results/Figures/BPCellTypeEffectSize_meth1_Revised.png")

################################################################################
# Reshape and output table 
#------------------------------------------------------------------------------#

rownames(test.res) <- NULL

tmp <- dplyr::filter(test.res, PE == "Adjusted for GA")

test.res <- merge(dplyr::filter(test.res, PE == "Main model"),
                  dplyr::filter(test.res, PE == "PE excluded"),
                  by = c("exposure", "window", "refCellType"))

test.res <- merge(test.res,
                  tmp,
                  by = c("exposure", "window", "refCellType"))

writexl::write_xlsx(test.res,
                    path = "Results/Tables/BPCellType_EffectSize_Revised.xlsx")


################################################################################
rm(test.res, tmp)
################################################################################

################################################################################
# Step C "Investigation" of the ratio Stromal/Syncytiotrophoblasts
# log(Stromal/Syncytiotrophoblast) ~ BP  + Confounders
################################################################################

Comp <- CCprop
# Choose reference cell types (as targeted mediator): 
refCellTypes <- colnames(Comp)

z1 <- log(Comp[,"Stromal"]/Comp[,"Syncytiotrophoblast"])

test.res <- data.frame()

for(w in seq_along(window)){
  
  BPhatFile <- paste0("Data/N668_BP1+covariates_", window[w], ".rds")
  
  BPhat <- readRDS( BPhatFile )
  BPhat <- BPhat %>% dplyr::filter(!(id %in% spuriousSamples))
  BPhat[,technical_confounders] <- apply(BPhat[,technical_confounders],2,as.factor)
  
  for(e in seq_along(exposure)){
    
    formula <- stats::as.formula(paste0("z1 ~  ", 
                                        exposure[e], " + center + batch + chip + plate + age + sex + BMI + smoke + parity"))
    
    data <- merge(BPhat, data.frame(id = rownames(Comp), z1))
    
    # Tests d'association
    res.lm1 <- lm(formula = formula, data = data)
    
    x <- summary(res.lm1)$coefficients
    y <- confint(res.lm1)[exposure[e],]
    z <- data.frame(t(x[exposure[e],'Pr(>|t|)']),
                    EffectSize = x[exposure[e],"Estimate"],
                    CI.low = y[1], CI.up = y[2])
    colnames(z)[1] <- "pvalue"
    
    test.res <- rbind.data.frame(test.res,
                                 data.frame(exposure = exposure[e], 
                                            window = window[w], 
                                            PE = "Main model",
                                            z))
    
    # Adjusting for Gestational Duration
    
    formula <- stats::as.formula(paste0("z1 ~  ", 
                                        exposure[e], " + center + batch + chip + plate + age + sex + BMI + smoke + parity + gestDuration"))
    
    data <- merge(BPhat, data.frame(id = rownames(Comp), z1))
    
    # Tests d'association
    res.lm1 <- lm(formula = formula, data = data)
    
    x <- summary(res.lm1)$coefficients
    y <- confint(res.lm1)[exposure[e],]
    z <- data.frame(t(x[exposure[e],'Pr(>|t|)']),
                    EffectSize = x[exposure[e],"Estimate"],
                    CI.low = y[1], CI.up = y[2])
    colnames(z)[1] <- "pvalue"
    
    test.res <- rbind.data.frame(test.res,
                                 data.frame(exposure = exposure[e], 
                                            window = window[w], 
                                            PE = "Adjusted for GA",
                                            z))
    
    # Excluding PE cases, main regression model
    data <- data %>% dplyr::filter(HDP != 2)
    
    formula <- stats::as.formula(paste0("z1 ~  ", 
                                        exposure[e], " + center + batch + chip + plate + age + sex + BMI + smoke + parity"))
    
    # Tests d'association
    res.lm1 <- lm(formula = formula, data = data)
    
    x <- summary(res.lm1)$coefficients
    y <- confint(res.lm1)[exposure[e],]
    z <- data.frame(t(x[exposure[e],'Pr(>|t|)']),
                    EffectSize = x[exposure[e],"Estimate"],
                    CI.low = y[1], CI.up = y[2])
    colnames(z)[1] <- "pvalue"
    
    test.res <- rbind.data.frame(test.res,
                                 data.frame(exposure = exposure[e], 
                                            window = window[w], 
                                            PE = "PE excluded",
                                            z))
    
  }
}

################################################################################
# Display effect sizes
#------------------------------------------------------------------------------#

test.res$exposure <- factor(test.res$exposure, levels = c("DBP","SBP","MAP","PP"))
test.res$PE <- factor(test.res$PE, levels = c("Main model", "PE excluded", "Adjusted for GA"))

test.res <- test.res %>% 
  dplyr::mutate(window = ifelse(window == "earlyP", "Early",
                                ifelse(window == "lateP", "Late", "Whole")))

test.res %>% ggplot(aes(x = window, y = EffectSize, col = PE)) +
  geom_pointrange(aes(ymin = CI.low, ymax = CI.up,
  ), size = 0.25, position = position_dodge(width = 0.7)) +
  geom_hline(yintercept = 0) +
  scale_color_manual(values = c("black", "grey63", "darkred"), "Analysis") +
  facet_grid( ~ exposure) +
  ggtitle("") + 
  ylab("Effect size") + xlab("Time Window") + 
  theme_minimal()

################################################################################
# Save graph of effect sizes
#------------------------------------------------------------------------------#

ggsave("Results/Figures/BPCellTypeEffectSize_ratio.pdf", width = 9, height = 3)

################################################################################
# Reshape and output table / Add as sheet to Supplementary Tables
#------------------------------------------------------------------------------#

rownames(test.res) <- NULL

tmp <- dplyr::filter(test.res, PE == "Adjusted for GA")

test.res <- merge(dplyr::filter(test.res, PE == "Main model"),
                  dplyr::filter(test.res, PE == "PE excluded"),
                  by = c("exposure", "window"))

test.res <- merge(test.res,
                  tmp,
                  by = c("exposure", "window"))

writexl::write_xlsx(test.res,
                    path = "Results/Tables/BPCellType_Ratio_Revised.xlsx")


################################################################################
rm(test.res, tmp)
################################################################################
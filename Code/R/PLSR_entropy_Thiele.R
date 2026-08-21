# This scripts computes Partial Least Square Regressions for entropy data 
# calculated according to Thiele et al. (2022) [1].
# Retrieved from: https://github.com/IS-UMK/complexity/tree/master/PLSR
# The original script was used in Dreszer et al. (2020) [2] and was adapted by
# Christoph Fruehlinger.
#
# [1] Thiele, J. A., Richter, A. & Hilger, K. (2023). Multimodal brain signal
# complexity predicts human intelligence. eNeuro, 10(2), ENEURO.0345-22.2022.
# https://doi.org/10.1523/eneuro.0345-22.2022
# 
# [2] Dreszer, J., Grochowski, M., Lewandowska, M., Nikadon, J., Gorgol, J., 
# Bałaj, B., Finc, K., Duch, W., Kałamała, P., Chuderski, A. & Piotrowski, T. 
# (2020). Spatiotemporal complexity patterns of resting‐state bioelectrical 
# activity explain fluid intelligence: Sex matters. Human Brain Mapping, 41(17),
# 4846–4865. https://doi.org/10.1002/hbm.25162
#
# Last edit: July 2025

if (!require("pacman")) install.packages("pacman")
pacman::p_load(pls, mvdalab, tidyverse, rlist, plspm, purrr)


options(scipen = 999)

rm(list = ls())

# Create output folder
savepath <- "Results/"

# load gf and demographic data
gf <- read.csv("Data/IST_table.csv", sep = ";")
gf$gf_score <- scale(rowSums(gf[,2:21]))

age <- read.csv("Data/SocioDemographics.txt")
age <- age %>%
  select(ID, Gender, Age) %>%
  mutate(Gender = case_when(Gender == 1 ~ "female",
                            Gender == 2 ~ "male"))

gf_dem <- inner_join(gf[,c("ID", "gf_score")], age, by = "ID")

# Load entropy data and define clusters
entropy_all <- read.csv("Data/Entropy_Thiele/df_entropy_first_run_eyes_open.csv")

# F cluster
F_cluster <- entropy_all %>%
  select(ID, Condition)

F_cluster$Set <- "F"

for(i in 0:20) {
  F_cluster[[paste0("MSE_", i)]] <- rowMeans(
    entropy_all[, c(
      paste0("MSE_F7_", i),
      paste0("MSE_F8_", i),
      paste0("MSE_F3_", i),
      paste0("MSE_F4_", i)
    )],
    na.rm = TRUE
  )
}

## AUC
x <- as.matrix(F_cluster[, paste0("MSE_", 0:11)])
F_cluster$AUC <- rowSums((x[, -1] + x[, -ncol(x)]) / 2)

## maxSlope
x <- as.matrix(F_cluster[, paste0("MSE_", 0:3)])
min_v  <- apply(x, 1, min)
max_v  <- apply(x, 1, max)
min_id <- apply(x, 1, which.min)
max_id <- apply(x, 1, which.max)
F_cluster$max_slope <- (max_v - min_v) / (max_id - min_id)

## avg_ent
F_cluster$avg_ent <- rowMeans(
  F_cluster[, paste0("MSE_", 8:11)],
  na.rm = TRUE
)

F_cluster$Shannon <- rowMeans(
  entropy_all[, c("shannon_ch_F7", "shannon_ch_F8", 
                  "shannon_ch_F3", "shannon_ch_F4")],
  na.rm = TRUE
)

F_cluster$Fuzzy <- rowMeans(
  entropy_all[, c("fuzzy_ch_F7", "fuzzy_ch_F8", 
                  "fuzzy_ch_F3", "fuzzy_ch_F4")],
  na.rm = TRUE
)


# FL cluster
FL_cluster <- entropy_all %>%
  select(ID, Condition)

FL_cluster$Set <- "FL"

for(i in 0:20) {
  FL_cluster[[paste0("MSE_", i)]] <- rowMeans(
    entropy_all[, c(
      paste0("MSE_FP1_", i),
      paste0("MSE_F7_", i),
      paste0("MSE_F3_", i),
      paste0("MSE_FC3_", i)
    )],
    na.rm = TRUE
  )
}

## AUC
x <- as.matrix(FL_cluster[, paste0("MSE_", 0:11)])
FL_cluster$AUC <- rowSums((x[, -1] + x[, -ncol(x)]) / 2)

## maxSlope
x <- as.matrix(FL_cluster[, paste0("MSE_", 0:3)])
min_v  <- apply(x, 1, min)
max_v  <- apply(x, 1, max)
min_id <- apply(x, 1, which.min)
max_id <- apply(x, 1, which.max)
FL_cluster$max_slope <- (max_v - min_v) / (max_id - min_id)

## avg_ent
FL_cluster$avg_ent <- rowMeans(
  FL_cluster[, paste0("MSE_", 8:11)],
  na.rm = TRUE
)

FL_cluster$Shannon <- rowMeans(
  entropy_all[, c("shannon_ch_FP1", "shannon_ch_F7", 
                  "shannon_ch_F3", "shannon_ch_FC3")],
  na.rm = TRUE
)

FL_cluster$Fuzzy <- rowMeans(
  entropy_all[, c("fuzzy_ch_FP1", "fuzzy_ch_F7", 
                  "fuzzy_ch_F3", "fuzzy_ch_FC3")],
  na.rm = TRUE
)


# FR cluster
FR_cluster  <- entropy_all %>%
  select(ID, Condition)

FR_cluster$Set <- "FR"

for(i in 0:20) {
  FR_cluster[[paste0("MSE_", i)]] <- rowMeans(
    entropy_all[, c(
      paste0("MSE_FP2_", i),
      paste0("MSE_F8_", i),
      paste0("MSE_F4_", i),
      paste0("MSE_FC4_", i)
    )],
    na.rm = TRUE
  )
}

## AUC
x <- as.matrix(FR_cluster[, paste0("MSE_", 0:11)])
FR_cluster$AUC <- rowSums((x[, -1] + x[, -ncol(x)]) / 2)

## maxSlope
x <- as.matrix(FR_cluster[, paste0("MSE_", 0:3)])
min_v  <- apply(x, 1, min)
max_v  <- apply(x, 1, max)
min_id <- apply(x, 1, which.min)
max_id <- apply(x, 1, which.max)
FR_cluster$max_slope <- (max_v - min_v) / (max_id - min_id)

## avg_ent
FR_cluster$avg_ent <- rowMeans(
  FR_cluster[, paste0("MSE_", 8:11)],
  na.rm = TRUE
)

FR_cluster$Shannon <- rowMeans(
  entropy_all[, c("shannon_ch_FP2", "shannon_ch_F8", 
                  "shannon_ch_F4", "shannon_ch_FC4")],
  na.rm = TRUE
)

FR_cluster$Fuzzy <- rowMeans(
  entropy_all[, c("fuzzy_ch_FP2", "fuzzy_ch_F2", 
                  "fuzzy_ch_F4", "fuzzy_ch_FC4")],
  na.rm = TRUE
)


# C cluster
C_cluster  <- entropy_all %>%
  select(ID, Condition)

C_cluster$Set <- "C"

for(i in 0:20) {
  C_cluster[[paste0("MSE_", i)]] <- rowMeans(
    entropy_all[, c(
      paste0("MSE_FZ_", i),
      paste0("MSE_CZ_", i),
      paste0("MSE_PZ_", i),
      paste0("MSE_OZ_", i)
    )],
    na.rm = TRUE
  )
}

## AUC
x <- as.matrix(C_cluster[, paste0("MSE_", 0:11)])
C_cluster$AUC <- rowSums((x[, -1] + x[, -ncol(x)]) / 2)

## maxSlope
x <- as.matrix(C_cluster[, paste0("MSE_", 0:3)])
min_v  <- apply(x, 1, min)
max_v  <- apply(x, 1, max)
min_id <- apply(x, 1, which.min)
max_id <- apply(x, 1, which.max)
C_cluster$max_slope <- (max_v - min_v) / (max_id - min_id)

## avg_ent
C_cluster$avg_ent <- rowMeans(
  C_cluster[, paste0("MSE_", 8:11)],
  na.rm = TRUE
)

C_cluster$Shannon <- rowMeans(
  entropy_all[, c("shannon_ch_FZ", "shannon_ch_CZ", 
                  "shannon_ch_PZ", "shannon_ch_OZ")],
  na.rm = TRUE
)

C_cluster$Fuzzy <- rowMeans(
  entropy_all[, c("fuzzy_ch_FZ", "fuzzy_ch_CZ", 
                  "fuzzy_ch_PZ", "fuzzy_ch_OZ")],
  na.rm = TRUE
)


# P Cluster
P_cluster  <- entropy_all %>%
  select(ID, Condition)

P_cluster$Set <- "P"

for(i in 0:20) {
  P_cluster[[paste0("MSE_", i)]] <- rowMeans(
    entropy_all[, c(
      paste0("MSE_P3_", i),
      paste0("MSE_P4_", i),
      paste0("MSE_P7_", i),
      paste0("MSE_P8_", i)
    )],
    na.rm = TRUE
  )
}

## AUC
x <- as.matrix(P_cluster[, paste0("MSE_", 0:11)])
P_cluster$AUC <- rowSums((x[, -1] + x[, -ncol(x)]) / 2)

## maxSlope
x <- as.matrix(P_cluster[, paste0("MSE_", 0:3)])
min_v  <- apply(x, 1, min)
max_v  <- apply(x, 1, max)
min_id <- apply(x, 1, which.min)
max_id <- apply(x, 1, which.max)
P_cluster$max_slope <- (max_v - min_v) / (max_id - min_id)

## avg_ent
P_cluster$avg_ent <- rowMeans(
  P_cluster[, paste0("MSE_", 8:11)],
  na.rm = TRUE
)

P_cluster$Shannon <- rowMeans(
  entropy_all[, c("shannon_ch_P3", "shannon_ch_P4", 
                  "shannon_ch_P7", "shannon_ch_P8")],
  na.rm = TRUE
)

P_cluster$Fuzzy <- rowMeans(
  entropy_all[, c("fuzzy_ch_P3", "fuzzy_ch_P4", 
                  "fuzzy_ch_P7", "fuzzy_ch_P8")],
  na.rm = TRUE
)


# PL Cluster
PL_cluster  <- entropy_all %>%
  select(ID, Condition)

PL_cluster$Set <- "PL"

for(i in 0:20) {
  PL_cluster[[paste0("MSE_", i)]] <- rowMeans(
    entropy_all[, c(
      paste0("MSE_P7_", i),
      paste0("MSE_P3_", i),
      paste0("MSE_O1_", i),
      paste0("MSE_PO3_", i)
    )],
    na.rm = TRUE
  )
}

## AUC
x <- as.matrix(PL_cluster[, paste0("MSE_", 0:11)])
PL_cluster$AUC <- rowSums((x[, -1] + x[, -ncol(x)]) / 2)

## maxSlope
x <- as.matrix(PL_cluster[, paste0("MSE_", 0:3)])
min_v  <- apply(x, 1, min)
max_v  <- apply(x, 1, max)
min_id <- apply(x, 1, which.min)
max_id <- apply(x, 1, which.max)
PL_cluster$max_slope <- (max_v - min_v) / (max_id - min_id)

## avg_ent
PL_cluster$avg_ent <- rowMeans(
  PL_cluster[, paste0("MSE_", 8:11)],
  na.rm = TRUE
)

PL_cluster$Shannon <- rowMeans(
  entropy_all[, c("shannon_ch_P7", "shannon_ch_P3", 
                  "shannon_ch_O1", "shannon_ch_PO3")],
  na.rm = TRUE
)

PL_cluster$Fuzzy <- rowMeans(
  entropy_all[, c("fuzzy_ch_P7", "fuzzy_ch_P3", 
                  "fuzzy_ch_O1", "fuzzy_ch_PO3")],
  na.rm = TRUE
)


# PR Cluster
PR_cluster  <- entropy_all %>%
  select(ID, Condition)

PR_cluster$Set <- "PR"

for(i in 0:20) {
  PR_cluster[[paste0("MSE_", i)]] <- rowMeans(
    entropy_all[, c(
      paste0("MSE_P8_", i),
      paste0("MSE_P4_", i),
      paste0("MSE_O2_", i),
      paste0("MSE_PO4_", i)
    )],
    na.rm = TRUE
  )
}

## AUC
x <- as.matrix(PR_cluster[, paste0("MSE_", 0:11)])
PR_cluster$AUC <- rowSums((x[, -1] + x[, -ncol(x)]) / 2)

## maxSlope
x <- as.matrix(PR_cluster[, paste0("MSE_", 0:3)])
min_v  <- apply(x, 1, min)
max_v  <- apply(x, 1, max)
min_id <- apply(x, 1, which.min)
max_id <- apply(x, 1, which.max)
PR_cluster$max_slope <- (max_v - min_v) / (max_id - min_id)

## avg_ent
PR_cluster$avg_ent <- rowMeans(
  PR_cluster[, paste0("MSE_", 8:11)],
  na.rm = TRUE
)

PR_cluster$Shannon <- rowMeans(
  entropy_all[, c("shannon_ch_P8", "shannon_ch_P4", 
                  "shannon_ch_O2", "shannon_ch_PO4")],
  na.rm = TRUE
)

PR_cluster$Fuzzy <- rowMeans(
  entropy_all[, c("fuzzy_ch_P8", "fuzzy_ch_P4", 
                  "fuzzy_ch_O2", "fuzzy_ch_PO4")],
  na.rm = TRUE
)


# ML Cluster
ML_cluster  <- entropy_all %>%
  select(ID, Condition)

ML_cluster$Set <- "ML"

for(i in 0:20) {
  ML_cluster[[paste0("MSE_", i)]] <- rowMeans(
    entropy_all[, c(
      paste0("MSE_T7_", i),
      paste0("MSE_C3_", i),
      paste0("MSE_CP5_", i),
      paste0("MSE_CP1_", i)
    )],
    na.rm = TRUE
  )
}

## AUC
x <- as.matrix(ML_cluster[, paste0("MSE_", 0:11)])
ML_cluster$AUC <- rowSums((x[, -1] + x[, -ncol(x)]) / 2)

## maxSlope
x <- as.matrix(ML_cluster[, paste0("MSE_", 0:3)])
min_v  <- apply(x, 1, min)
max_v  <- apply(x, 1, max)
min_id <- apply(x, 1, which.min)
max_id <- apply(x, 1, which.max)
ML_cluster$max_slope <- (max_v - min_v) / (max_id - min_id)

## avg_ent
ML_cluster$avg_ent <- rowMeans(
  ML_cluster[, paste0("MSE_", 8:11)],
  na.rm = TRUE
)

ML_cluster$Shannon <- rowMeans(
  entropy_all[, c("shannon_ch_T7", "shannon_ch_C3", 
                  "shannon_ch_CP5", "shannon_ch_CP1")],
  na.rm = TRUE
)

ML_cluster$Fuzzy <- rowMeans(
  entropy_all[, c("fuzzy_ch_T7", "fuzzy_ch_C3", 
                  "fuzzy_ch_CP5", "fuzzy_ch_CP1")],
  na.rm = TRUE
)


# MR Cluster
MR_cluster  <- entropy_all %>%
  select(ID, Condition)

MR_cluster$Set <- "MR"

for(i in 0:20) {
  MR_cluster[[paste0("MSE_", i)]] <- rowMeans(
    entropy_all[, c(
      paste0("MSE_T8_", i),
      paste0("MSE_C4_", i),
      paste0("MSE_CP6_", i),
      paste0("MSE_CP2_", i)
    )],
    na.rm = TRUE
  )
}

## AUC
x <- as.matrix(MR_cluster[, paste0("MSE_", 0:11)])
MR_cluster$AUC <- rowSums((x[, -1] + x[, -ncol(x)]) / 2)

## maxSlope
x <- as.matrix(MR_cluster[, paste0("MSE_", 0:3)])
min_v  <- apply(x, 1, min)
max_v  <- apply(x, 1, max)
min_id <- apply(x, 1, which.min)
max_id <- apply(x, 1, which.max)
MR_cluster$max_slope <- (max_v - min_v) / (max_id - min_id)

## avg_ent
MR_cluster$avg_ent <- rowMeans(
  MR_cluster[, paste0("MSE_", 8:11)],
  na.rm = TRUE
)

MR_cluster$Shannon <- rowMeans(
  entropy_all[, c("shannon_ch_T8", "shannon_ch_C4", 
                  "shannon_ch_CP6", "shannon_ch_CP2")],
  na.rm = TRUE
)

MR_cluster$Fuzzy <- rowMeans(
  entropy_all[, c("fuzzy_ch_T8", "fuzzy_ch_C4", 
                  "fuzzy_ch_CP6", "fuzzy_ch_CP2")],
  na.rm = TRUE
)

# Merge dataframes
mmse_data <- dplyr::bind_rows(
  F_cluster,
  FL_cluster,
  FR_cluster,
  C_cluster,
  P_cluster,
  PL_cluster,
  PR_cluster,
  ML_cluster,
  MR_cluster
)

entropy_data_all <- inner_join(gf_dem, mmse_data, by = "ID")

# Define main condition to analyse 
# first_run_eyes_open was the preregistered condition. We analyzed the remaining
# exploratively. Just copy-paste the condition you are interested in.

main_cond <- "first_run_eyes_open"

if (!dir.exists(savepath)) {
  dir.create(savepath)
}

# List of analysis conditions
entropies <- c("mMSE", "Shannon", "Fuzzy")
subsamples <- c("Full", "Male", "Female")

# loop through data types and subsamples

cat("***Starting PLSR***\n")

for (entropy in entropies) {
  for (sample in subsamples) {
    
    # get rid of unnecessary information
    entropy_data <- entropy_data_all %>%
      filter(
        # Length == 40,
        # Condition == main_cond,
        !is.na(Condition) & Condition != "<undefined>",
        !is.na(Set) & Set != "<undefined>",
        rowSums(!is.na(across(-c(ID, gf_score, Gender, Age,
                                 Condition, Set)))) > 0
      )
    
    # transform to wide format
    entropy_data <- entropy_data %>%
      mutate(SetCond = paste(Set, Condition, sep = "_"))
    
    if (entropy == "mMSE") {
      value_cols <- c("AUC", "max_slope", "avg_ent",
                      grep("^mmse_", names(entropy_data), value = TRUE))
      
      entropy_data_wide <- entropy_data %>%
        pivot_wider(
          id_cols = c(ID, gf_score, Gender, Age),
          names_from = SetCond,
          values_from = all_of(value_cols)
        )
    } else if (entropy == "Shannon") {
      
      entropy_data_wide <- entropy_data %>%
        pivot_wider(
          id_cols = c(ID, gf_score, Gender, Age),
          names_from = SetCond,
          values_from = all_of("Shannon")
        )
      
    } else if (entropy == "Fuzzy") {
      
      entropy_data_wide <- entropy_data %>%
        pivot_wider(
          id_cols = c(ID, gf_score, Gender, Age),
          names_from = SetCond,
          values_from = all_of("Fuzzy")
        )
      
    }
    
    # Select sample data
    if (entropy == "mMSE") {
      if (sample != "Full") {
        
        entropy_data_wide_subsample <- entropy_data_wide %>%
          filter(Gender == tolower(sample))
        
        X_AUC <- as.matrix(entropy_data_wide_subsample %>% 
                             select(contains("auc")))
        
        X_AVG <- as.matrix(entropy_data_wide_subsample %>% 
                             select(contains("avg_entropy")))
        
        X_MaxSlope <- as.matrix(entropy_data_wide_subsample %>% 
                                  select(contains("max_slope")))
        
        X <- cbind(X_AUC,X_AVG,X_MaxSlope)
        
        y <- entropy_data_wide_subsample$gf_score
        
      } else {
        
        entropy_data_wide_subsample <- entropy_data_wide
        
        X_AUC <- as.matrix(entropy_data_wide_subsample %>% 
                             select(contains("AUC")))
        
        X_AVG <- as.matrix(entropy_data_wide_subsample %>% 
                             select(contains("avg_ent")))
        
        X_MaxSlope <- as.matrix(entropy_data_wide_subsample %>% 
                                  select(contains("max_slope")))
        
        X <- cbind(X_AUC,X_AVG,X_MaxSlope)
        
        y <- entropy_data_wide_subsample$gf_score
        
      } # end for if sample != "Full"
    } else if (entropy == "Shannon") {
      if (sample != "Full") {
        
        entropy_data_wide_subsample <- entropy_data_wide %>%
          filter(Gender == tolower(sample))
        
        X <- as.matrix(entropy_data_wide_subsample %>% 
                             select(-c("ID", "gf_score", "Gender", "Age")))
        
        y <- entropy_data_wide_subsample$gf_score
        
      } else {
        
        entropy_data_wide_subsample <- entropy_data_wide
        
        X <- as.matrix(entropy_data_wide_subsample %>% 
                                 select(-c("ID", "gf_score", "Gender", "Age")))
        
        y <- entropy_data_wide_subsample$gf_score
        
      } # end for if sample != "Full"
      
    } else if (entropy == "Fuzzy") {
      
      if (sample != "Full") {
        
        entropy_data_wide_subsample <- entropy_data_wide %>%
          filter(Gender == tolower(sample))
        
        X <- as.matrix(entropy_data_wide_subsample %>% 
                                 select(-c("ID", "gf_score", "Gender", "Age")))
        
        y <- entropy_data_wide_subsample$gf_score
        
      } else {
        
        entropy_data_wide_subsample <- entropy_data_wide
        
        X <- as.matrix(entropy_data_wide_subsample %>% 
                                 select(-c("ID", "gf_score", "Gender", "Age")))
        
        y <- entropy_data_wide_subsample$gf_score
        
      } # end for if sample != "Full"
      
    }
    
    # cases used for training
    PLSmodel <- list(y,X)
    names(PLSmodel) <- list('dep','indep')
    
    cat("***************************\n")
    cat(paste0("N = ", nrow(X)), "\n")
    
    # maximum number of components for PLS models
    max_comp <- ncol(X)
    
    # PLSR model with CV to determine number of components
    PLSR_CV <- plsr(dep ~ indep, ncomp = max_comp, data = PLSmodel, 
                    na.action = na.omit, method = "simpls", 
                    validation = "LOO", jackknife = TRUE, scale = TRUE)
    
    numberSignifComp <- selectNcomp(PLSR_CV, method = "randomization", 
                                    alpha=0.05, nperm=10000)
    
    # Check for absolute minimum
    abs_min <- which.min(RMSEP(PLSR_CV)$val["adjCV", 1, ]) - 1
    
    if (numberSignifComp == 0 && abs_min == 0) {
      
      cat(" Entropy: ", entropy, "\n",
          "Sample: ", sample, "\n",
          r"(No significant components selected ¯\_(ツ)_/¯)", "\n")
      
      cat("Abs. Min: ", which.min(RMSEP(PLSR_CV)$val["adjCV", 1, ]) - 1, '\n')
      
      # tiff(paste(savepath, "PLSR_selected_comps_", main_cond, "_", data, "_", 
      #            sample, ".tiff", sep = ""), width = 2400, height = 1800, 
      #      res = 600)
      # selectNcomp(PLSR_CV, method = "randomization", alpha=0.05, nperm=10000, 
      #             plot = TRUE, ylim = c(0.9, 1.5))
      # dev.off()
      
    } else {
      
      cat("Entropy: ", entropy, "\n",
          "Sample: ", sample, "\n",
          "Abs. Min: ", abs_min, "\n",
          "Sign. Component: ", numberSignifComp, "\n", sep = "")
      
      # tiff(paste(savepath, "PLSR_selected_comps_", main_cond, "_", data, "_", 
      #            sample, ".tiff", sep = ""), width = 2400, height = 1800, 
      #      res = 600)
      # selectNcomp(PLSR_CV, method = "randomization", alpha=0.05, nperm=10000, 
      #             plot = TRUE, ylim = c(0.9, 1.5))
      # dev.off()
      
      # Select component
      selected_comp <- max(numberSignifComp, abs_min)
      
      # Jack-knife estimator of regression coefficients
      jack.test(PLSR_CV,ncomp=selected_comp)
      
      # R-squared and RMSEP for CV model
      cat(r"(R2: PLSR_CV)", "\n")
      print(R2(PLSR_CV))
      cat(r"(RMSEP: PLSR_CV)", "\n")
      print(RMSEP(PLSR_CV))
      
      # fitted (fixed-effect) PLS model with preselected number of components
      PLSR_NO_CV <- plsr(dep ~ indep, ncomp = max_comp, data = PLSmodel, 
                         method = "simpls", scale = TRUE)
      
      # R-squared and RMSEP for fitted model
      R2(PLSR_NO_CV)
      RMSEP(PLSR_NO_CV)
      
      # fitted (fixed-effect) PLS model using MVDALAB function allowing for 
      # bootstrapped confidence interval estimation 
      PLSR_BOOT <- plsFit(dep ~ indep, scale = TRUE, data = PLSmodel, 
                          ncomp = selected_comp, validation = "oob", 
                          boots = 10000)
      
      # R-squared for fitted model
      cat(r"(R2: PLSR_BOOT)", "\n")
      print(R2s(PLSR_BOOT))
      
      # bootstrapped confidence interval estimation
      cfs_int <- coefficients.boots(PLSR_BOOT, ncomp = 1, 
                                    conf = (1 - (0.05 / 3)))
      
      # all intervals
      x <- cfs_int[[1]]
      
      # significant regression coefficients obtained from bootstrapped 
      # confidence interval estimation with alpha=0.05
      print(filter(x, `0.8333333%` > 0 | `99.16667%` < 0))
      
    } # end for if numberSignifComp == 0
    
    # Clear variables
    rm(X_AUC, X_AVG, X_MaxSlope, X, y, numberSignifComp, PLSmodel, PLSR_CV)
    
  } # end for subsample
}


# # PLSPM Model
# 
# # MMSE
# 
# # Define inner model
# auc <- c(0,0,0,0)
# max_slope <- c(0,0,0,0)
# avg_ent <- c(0,0,0,0)
# gf <- c(1,1,1,0)
# 
# mmse_path <- rbind(auc,max_slope, avg_ent, gf)
# colnames(mmse_path) <- rownames(mmse_path)
# 
# innerplot(mmse_path)
# 
# # Define outer model
# mmse_blocks <- list(
#   colnames(mmse_data_wide_subsample %>% select(contains("AUC"))), 
#   colnames(mmse_data_wide_subsample %>% select(contains("max_slope"))), 
#   colnames(mmse_data_wide_subsample %>% select(contains("avg_ent"))), 
#   "gf_score")
# 
# # Define mode
# mmse_mode <- rep("A",4)
# 
# # Filter for Gender and run model separately
# mmse_data_wide_plspm_female <-  mmse_data_wide %>%
#   filter(Gender == "female") %>%
#   drop_na(all_of(unlist(mmse_blocks)))
# 
# mmse_pls_female <- plspm(mmse_data_wide_plspm_female, mmse_path, mmse_blocks, 
#                          modes = mmse_mode)
# 
# mmse_data_wide_plspm_male <-  mmse_data_wide %>%
#   filter(Gender == "male") %>%
#   drop_na(all_of(unlist(mmse_blocks)))
# 
# mmse_pls_male <- plspm(mmse_data_wide_plspm_male, mmse_path, mmse_blocks, 
#                        modes = mmse_mode)
# 
# # Run Group Model
# mmse_data_wide_plspm <-  mmse_data_wide %>%
#   drop_na(all_of(unlist(mmse_blocks))) %>%
#   filter(Gender != "<undefined>")
# 
# mmse_data_wide_plspm$Gender <- as.factor(mmse_data_wide_plspm$Gender)
# 
# mmse_pls <- plspm(mmse_data_wide_plspm, mmse_path, mmse_blocks, 
#                   modes = mmse_mode)
# 
# mmse_plspm_boot <- plspm.groups(mmse_pls, mmse_data_wide_plspm$Gender, 
#                                 method = "bootstrap")
# 
# mmse_plspm_boot
# 
# bar_colors <- c("female" = "#FEB24C", "male" = "#74A9CF")
# 
# bar_mmse <- mmse_plspm_boot$test[, 2:3] %>%
#   rownames_to_column("Path") %>%
#   pivot_longer(-Path, names_to = "Group", values_to = "Coefficient") %>%
#   mutate(Group = case_when(
#     Group == "group.female" ~ "female",
#     Group == "group.male" ~ "male",
#     TRUE ~ Group
#   ))
# 
# bar_mmse$Path <- sub("->gf", "", bar_mmse$Path)
# bar_mmse$Path <- sub("_", ". ", bar_mmse$Path)
# bar_mmse$Path <- sub("ent", "entropy", bar_mmse$Path)
# bar_mmse$Path[1:2] <- toupper(bar_mmse$Path[1:2])
# bar_mmse$Path[3:6] <- tools::toTitleCase(bar_mmse$Path[3:6])
# 
# plspm_plot_mmse <- ggplot(bar_mmse, aes(x = Path, y = Coefficient, 
#                                         fill = Group)) +
#   geom_bar(stat = "identity", position = position_dodge(width = 0.8), 
#            width = 0.7) +
#   geom_hline(yintercept = 0, color = "gray50") +
#   scale_fill_manual(values = bar_colors, name = NULL) +
#   scale_y_continuous() +
#   labs(x = NULL, y = "Path Coefficient") +
#   theme_classic(base_size = 14) +
#   theme(
#     panel.grid.major = element_line(color = "grey80"),
#     panel.grid.minor = element_line(color = "grey90"),
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     axis.text = element_text(color = "gray30"),
#     legend.position = "right",
#     legend.text = element_text(color = "gray40")
#   )
# 
# plspm_plot_mmse
# 
# ggsave(filename = paste0(savepath, "PLSPM_MMSE_" , main_cond, ".tiff"), 
#        plot = plspm_plot_mmse, width = 8, height = 5, dpi = 600)
# 
# 
# # Microstates
# 
# # Define inner model
# n_peaks <- c(0,0,0,0,0,0,0,0)
# coverage <- c(0,0,0,0,0,0,0,0)
# lifespan <- c(0,0,0,0,0,0,0,0)
# lifespan_peaks <- c(0,0,0,0,0,0,0,0)
# frequency <- c(0,0,0,0,0,0,0,0)
# transition_probability <- c(0,0,0,0,0,0,0,0)
# transition_probability_peaks <- c(0,0,0,0,0,0,0,0)
# gf <- c(1,1,1,1,1,1,1,0)
# 
# microstate_path <- rbind(n_peaks, coverage, lifespan, lifespan_peaks, frequency,
#                          transition_probability, transition_probability_peaks, gf)
# colnames(microstate_path) <- rownames(microstate_path)
# 
# innerplot(microstate_path)
# 
# # Define outer model
# microstate_blocks <- list(
#   colnames(microstate_data %>% select(contains("n_gfp_peaks"))), 
#   colnames(microstate_data %>% select(contains("coverage"))), 
#   colnames(microstate_data %>% select(contains("lifespan"))),
#   colnames(microstate_data %>% select(contains("lifespan_peaks"))),
#   colnames(microstate_data %>% select(contains("frequence"))),
#   colnames(microstate_data %>% select(contains("transition_probability"))),
#   colnames(microstate_data %>% select(contains("transition_probability_peaks"))), 
#   "gf_score")
# 
# # Define mode
# microstate_mode <- rep("A",8)
# 
# # Filter for Gender and run model separately
# microstate_data_plspm_female <-  microstate_data %>%
#   filter(Gender == "female") %>%
#   drop_na(all_of(unlist(microstate_blocks)))
# 
# microstate_pls_female <- plspm(microstate_data_plspm_female, microstate_path, 
#                                microstate_blocks, modes = microstate_mode)
# 
# microstate_data_plspm_male <-  microstate_data %>%
#   filter(Gender == "male") %>%
#   drop_na(all_of(unlist(microstate_blocks)))
# 
# microstate_pls_male <- plspm(microstate_data_plspm_male, microstate_path, 
#                              microstate_blocks, modes = microstate_mode)
# 
# # Run Group Model
# microstate_data_plspm <-  microstate_data %>%
#   drop_na(all_of(unlist(microstate_blocks))) %>%
#   filter(Gender != "<undefined>")
# 
# microstate_data_plspm$Gender <- as.factor(microstate_data_plspm$Gender)
# 
# microstate_pls <- plspm(microstate_data_plspm, microstate_path, 
#                         microstate_blocks, modes = microstate_mode)
# 
# microstate_plspm_boot <- plspm.groups(microstate_pls, microstate_data_plspm$Gender, 
#                                       method = "bootstrap")
# 
# microstate_plspm_boot
# 
# bar_colors <- c("female" = "#FEB24C", "male" = "#74A9CF")
# 
# bar_ms <- microstate_plspm_boot$test[, 2:3] %>%
#   rownames_to_column("Path") %>%
#   pivot_longer(-Path, names_to = "Group", values_to = "Coefficient") %>%
#   mutate(Group = case_when(
#     Group == "group.female" ~ "female",
#     Group == "group.male" ~ "male",
#     TRUE ~ Group
#   ))
# 
# bar_ms$Path <- sub("->gf", "", bar_ms$Path)
# bar_ms$Path <- sub("_", " ", bar_ms$Path)
# bar_ms$Path <- sub("_", " ", bar_ms$Path)
# bar_ms$Path <- tools::toTitleCase(bar_ms$Path)
# bar_ms$Path[1:2] <- sub("n", "Number of GFP", bar_ms$Path[1:2])
# 
# plspm_plot_microstate <- ggplot(bar_ms, aes(x = Path, y = Coefficient, 
#                                             fill = Group)) +
#   geom_bar(stat = "identity", position = position_dodge(width = 0.8), 
#            width = 0.7) +
#   geom_hline(yintercept = 0, color = "gray50") +
#   scale_fill_manual(values = bar_colors, name = NULL) +
#   scale_y_continuous() +
#   labs(x = NULL, y = "Path Coefficient") +
#   theme_classic(base_size = 14) +
#   theme(
#     panel.grid.major = element_line(color = "grey80"),
#     panel.grid.minor = element_line(color = "grey90"),
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     axis.text = element_text(color = "gray30"),
#     legend.position = "right",
#     legend.text = element_text(color = "gray40")
#   )
# 
# plspm_plot_microstate
# 
# ggsave(filename = paste0(savepath, "PLSPM_Microstates_" , main_cond, ".tiff"), 
#        plot = plspm_plot_microstate, width = 8, height = 5, dpi = 600)

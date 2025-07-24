# This scripts computes Partial Least Square Regressions for MMSE and Microstate
# data. Retrieved from: https://github.com/IS-UMK/complexity/tree/master/PLSR
# The original script was used in Dreszer et al. (2020) [1] and was adapted by
# Christoph Fruehlinger.
#
# [1] Dreszer, J., Grochowski, M., Lewandowska, M., Nikadon, J., Gorgol, J., 
# Bałaj, B., Finc, K., Duch, W., Kałamała, P., Chuderski, A. & Piotrowski, T. 
# (2020). Spatiotemporal complexity patterns of resting‐state bioelectrical 
# activity explain fluid intelligence: Sex matters. Human Brain Mapping, 41(17),
# 4846–4865. https://doi.org/10.1002/hbm.25162
#
# Last edit: July 2025

if (!require("pacman")) install.packages("pacman")
pacman::p_load(pls, mvdalab, tidyverse, rlist)


options(scipen = 999)

rm(list = ls())

# Create output folder
savepath = "Results/"

# Define main condition to analyse
main_cond = "first_run_eyes_open"

if (!dir.exists(savepath)) {
  dir.create(savepath)
}

# List of analysis conditions
data_types = c("MMSE", "Microstates")
subsamples = c("Full", "Male", "Female")

# loop through data types and subsamples

cat("***Starting PLSR***\n")

for (data in data_types) {
  
  for (sample in subsamples) {
    
    cat(paste0("\n", "Data-Type: ", data, "; Subsample: ", sample, "\n"))
    
    if (data == "MMSE") {
      ## MMSE data
      # load data set
      mmse_datapath = "Data/MMSEData/MMSE_data_full.csv"
      
      # Check if file exists else continue with Microstates data
      if (!file.exists(mmse_datapath)) {
        
        cat("No MMSE data found...\n")
        next
        
      }
      
      mmse_data <- read_csv(mmse_datapath, show_col_types = FALSE)
      
      # get rid of unnecessary information
      mmse_data <- mmse_data %>%
        filter(
          Length == 40,
          Condition == main_cond,
          !is.na(Condition) & Condition != "<undefined>",
          !is.na(Set) & Set != "<undefined>",
          rowSums(!is.na(across(-c(ID, gf_score, Gender, Age,
                                   Condition, Set)))) > 0
        )
      
      # transform to wide format
      mmse_data <- mmse_data %>%
        mutate(SetCond = paste(Set, Condition, sep = "_"))
      
      
      value_cols <- c("auc", "max_slope", "avg_entropy",
                      grep("^mmse_", names(mmse_data), value = TRUE))
      
      mmse_data_wide <- mmse_data %>%
        pivot_wider(
          id_cols = c(ID, gf_score, Gender, Age),
          names_from = SetCond,
          values_from = all_of(value_cols)
        )
      
      # Select sample data
      if (sample != "Full") {
        
        mmse_data_wide_subsample <- mmse_data_wide %>%
          filter(Gender == tolower(sample))
        
        X_AUC <- as.matrix(mmse_data_wide_subsample %>% 
                             select(contains("auc")))
        
        X_AVG <- as.matrix(mmse_data_wide_subsample %>% 
                             select(contains("avg_entropy")))
        
        X_MaxSlope <- as.matrix(mmse_data_wide_subsample %>% 
                                  select(contains("max_slope")))
        
        X <- cbind(X_AUC,X_AVG,X_MaxSlope)
        
        y <- mmse_data_wide_subsample$gf_score
        
      } else {
        
        mmse_data_wide_subsample <- mmse_data_wide
        
        X_AUC <- as.matrix(mmse_data_wide_subsample %>% 
                             select(contains("auc")))
        
        X_AVG <- as.matrix(mmse_data_wide_subsample %>% 
                             select(contains("avg_entropy")))
        
        X_MaxSlope <- as.matrix(mmse_data_wide_subsample %>% 
                                  select(contains("max_slope")))
        
        X <- cbind(X_AUC,X_AVG,X_MaxSlope)
        
        y <- mmse_data_wide_subsample$gf_score
        
      } # end for if sample != "Full"
      
      # cases used for training
      PLSmodel <- list(y,X)
      names(PLSmodel) <- list('dep','indep')
      
      cat(paste0("N = ", nrow(X)), "\n")
      
      # maximum number of components for PLS models
      max_comp=30
      
      # PLSR model with CV to determine number of components
      PLSR_CV <- plsr(dep ~ indep, ncomp = max_comp, data = PLSmodel, 
                      na.action = na.omit, method = "simpls", 
                      validation = "LOO", jackknife = TRUE, scale = TRUE)
      
      numberSignifComp <- selectNcomp(PLSR_CV, method = "randomization", 
                                      alpha=0.05, nperm=10000)
      
      # Check for absolute minimum
      abs_min <- which.min(RMSEP(PLSR_CV)$val["adjCV", 1, ]) - 1
      
      if (numberSignifComp == 0 && abs_min == 0) {
        
        cat(r"(No significant components selected ¯\_(ツ)_/¯)", "\n")
        
        cat("Abs. Min: ", which.min(RMSEP(PLSR_CV)$val["adjCV", 1, ]) - 1, '\n')
        
        tiff(paste(savepath, "PLSR_selected_comps_", data, "_", sample, ".tiff", 
                   sep = ""), width = 2400, height = 1800, res = 600)
        selectNcomp(PLSR_CV, method = "randomization", alpha=0.05, nperm=10000, 
                    plot = TRUE, ylim = c(0.9, 1.5))
        dev.off()
        
      } else {
        
        # Select component
        selected_comp = max(numberSignifComp, abs_min)
        
        # Jack-knife estimator of regression coefficients
        jack.test(PLSR_CV,ncomp=numberSignifComp)
        
        # R-squared and RMSEP for CV model
        R2(PLSR_CV)
        RMSEP(PLSR_CV)
        
        # fitted (fixed-effect) PLS model with preselected number of components
        PLSR_NO_CV <- plsr(dep ~ indep, ncomp = max_comp, data = PLSmodel, 
                           method = "simpls", scale = TRUE)
        
        # R-squared and RMSEP for fitted model
        R2(PLSR_NO_CV)
        RMSEP(PLSR_NO_CV)
        
        # fitted (fixed-effect) PLS model using MVDALAB function allowing for 
        # bootstrapped confidence interval estimation 
        PLSR_BOOT <- plsFit(dep ~ indep, scale = TRUE, data = PLSmodel, 
                            ncomp = numberSignifComp, validation = "oob", 
                            boots = 10000)
        
        # R-squared for fitted model
        R2s(PLSR_BOOT)
        
        # bootstrapped confidence interval estimation
        cfs_int <- coefficients.boots(PLSR_BOOT, ncomp = 1, conf = .95)
        
        # all intervals
        x <- cfs_int[[1]]
        x
        
        # significant regression coefficients obtained from bootstrapped 
        # confidence interval estimation with alpha=0.05
        filter(x,x$"2.5%">0)
        filter(x,x$"97.5%"<0)
        
        filter(x,x$"2.5%">0 & x$"97.5" < 0)
        
      } # end for if numberSignifComp == 0
      
      # Clear variables
      rm(X_AUC, X_AVG, X_MaxSlope, X, y, numberSignifComp, PLSmodel, PLSR_CV)
      
    } else {
      ## Microstate data
      microstate_datapath = "Data/Microstates/Microstates_data_full.csv"
      
      # Check if file exists
      if (!file.exists(microstate_datapath)) {
        
        cat("No microstate data found...\n")
        break
        
      }
      
      microstate_data <- read_csv(microstate_datapath, show_col_types = FALSE)
      
      microstate_names <- c("0" = "F", "1" = "C", "2" = "D", "3" = "B", "4" = "A")
      
      colnames(microstate_data) <- str_replace_all(colnames(microstate_data), 
                                                   microstate_names)
      
      # get rid of unnecessary information
      microstate_data <- microstate_data %>%
        filter(
          Length == 40,
          Condition == main_cond,
          !is.na(Condition) & Condition != "<undefined>",
          rowSums(!is.na(across(-c(ID, gf_score, Gender, Age, Condition)))) > 0
        )
      
      # Select sample data
      if (sample != "Full") {
        
        microstate_data_subsample <- microstate_data %>%
          filter(Gender == tolower(sample))
        
        X_n_gfp_peaks <- as.matrix(microstate_data_subsample %>%
                                     select(contains("n_gfp_peaks")))
        
        X_coverage <- as.matrix(microstate_data_subsample %>% 
                             select(contains("coverage")))
        
        X_lifespan <- as.matrix(microstate_data_subsample %>%
                                  select(contains("lifespan") 
                                         & !contains("lifespan_peaks")))
        
        X_lifespan_peaks <- as.matrix(microstate_data_subsample %>% 
                                  select(contains("lifespan_peaks")))
        
        X_frequency <- as.matrix(microstate_data_subsample %>% 
                                  select(contains("frequence")))
        
        X_trans_prob <- as.matrix(microstate_data_subsample %>%
                                    select(contains("transition_probability") 
                                           & !contains("transition_probability_peaks")))
        
        X_trans_prob_peaks <- as.matrix(microstate_data_subsample %>% 
                                        select(contains("transition_probability_peaks")))
        
        X <- cbind(X_n_gfp_peaks, X_coverage, X_lifespan, X_lifespan_peaks,
                   X_frequency, X_trans_prob, X_trans_prob_peaks)
        
        y <- microstate_data_subsample$gf_score
        
      } else {
        
        microstate_data_subsample <- microstate_data
        
        X_n_gfp_peaks <- as.matrix(microstate_data_subsample %>%
                                     select(contains("n_gfp_peaks")))
        
        X_coverage <- as.matrix(microstate_data_subsample %>% 
                                  select(contains("coverage")))
        
        X_lifespan <- as.matrix(microstate_data_subsample %>%
                                  select(contains("lifespan") 
                                         & !contains("lifespan_peaks")))
        
        X_lifespan_peaks <- as.matrix(microstate_data_subsample %>% 
                                        select(contains("lifespan_peaks")))
        
        X_frequency <- as.matrix(microstate_data_subsample %>% 
                                   select(contains("frequence")))
        
        X_trans_prob <- as.matrix(microstate_data_subsample %>%
                                    select(contains("transition_probability") 
                                           & !contains("transition_probability_peaks")))
        
        X_trans_prob_peaks <- as.matrix(microstate_data_subsample %>% 
                                          select(contains("transition_probability_peaks")))
        
        X <- cbind(X_n_gfp_peaks, X_coverage, X_lifespan, X_lifespan_peaks,
                   X_frequency, X_trans_prob, X_trans_prob_peaks)
        
        y <- microstate_data_subsample$gf_score
        
      } # end for if sample != "Full"
      
      # cases used for training
      PLSmodel <- list(y,X)
      names(PLSmodel) <- list('dep','indep')
      
      cat(paste0("N = ", nrow(X)), "\n")
      
      # maximum number of components for PLS models
      max_comp=30
      
      # PLSR model with CV to determine number of components
      PLSR_CV <- plsr(dep ~ indep, ncomp = max_comp, data = PLSmodel, 
                      na.action = na.omit, method = "simpls", 
                      validation = "LOO", jackknife = TRUE, scale = TRUE)
      
      numberSignifComp <- selectNcomp(PLSR_CV, method = "randomization", 
                                      alpha=0.05, nperm=10000)
      
      # Check for absolute minimum
      abs_min <- which.min(RMSEP(PLSR_CV)$val["adjCV", 1, ]) - 1
      
      if (numberSignifComp == 0 && abs_min == 0) {
        
        cat(r"(No significant components selected ¯\_(ツ)_/¯)", "\n")
        
        cat("Abs. Min:", which.min(RMSEP(PLSR_CV)$val["adjCV", 1, ]) - 1, '\n')
        
        tiff(paste(savepath, "PLSR_selected_comps_", data, "_", sample, ".tiff", 
                   sep = ""), width = 2400, height = 1800, res = 600)
        selectNcomp(PLSR_CV, method = "randomization", alpha=0.05, nperm=10000, 
                    plot = TRUE, ylim = c(0.9, 1.5))
        dev.off()
        
      } else {
        
        # Select component
        selected_comp = max(numberSignifComp, abs_min)
        
        # Jack-knife estimator of regression coefficients
        jack.test(PLSR_CV,ncomp=selected_comp)
        
        # R-squared and RMSEP for CV model
        R2(PLSR_CV)
        RMSEP(PLSR_CV)
        
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
        R2s(PLSR_BOOT)
        
        # bootstrapped confidence interval estimation
        cfs_int <- coefficients.boots(PLSR_BOOT, ncomp = 1, conf = .95)
        
        # all intervals
        x <- cfs_int[[1]]
        x
        
        # significant regression coefficients obtained from bootstrapped 
        # confidence interval estimation with alpha=0.05
        filter(x,x$"2.5%">0)
        filter(x,x$"97.5%"<0)
        
        filter(x,x$"2.5%">0 | x$"97.5%" < 0)
        
      } # end for if numberSignifComp == 0
      
    } # end for if data == "MMSE"
    
  } # end for subsample
  
} # end for data type
  
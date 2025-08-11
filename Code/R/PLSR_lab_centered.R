# This scripts computes Partial Least Square Regressions for lab-mean-centered 
# MMSE and Microstate data. 
# Retrieved from: https://github.com/IS-UMK/complexity/tree/master/PLSR
# The original script was used in Dreszer et al. (2020) [1] and was adapted by
# Christoph Fruehlinger.
#
# [1] Dreszer, J., Grochowski, M., Lewandowska, M., Nikadon, J., Gorgol, J., 
# Bałaj, B., Finc, K., Duch, W., Kałamała, P., Chuderski, A. & Piotrowski, T. 
# (2020). Spatiotemporal complexity patterns of resting‐state bioelectrical 
# activity explain fluid intelligence: Sex matters. Human Brain Mapping, 41(17),
# 4846–4865. https://doi.org/10.1002/hbm.25162
#
# Last edit: August 2025

if (!require("pacman")) install.packages("pacman")
pacman::p_load(pls, mvdalab, tidyverse, rlist, plspm)


options(scipen = 999)

rm(list = ls())

# Create output folder
savepath <- "Results/"

# Define main condition to analyse 
# first_run_eyes_open was the preregistered condition. We analyzed the remaining
# exploratively. Just copy-paste the condition you are interested in.

main_cond <- "first_run_eyes_open"

if (!dir.exists(savepath)) {
  dir.create(savepath)
}

lab_data <- read_tsv("Data/Labs.txt", col_select = c(1, 2))

# List of analysis conditions
data_types <- c("MMSE", "Microstates")
subsamples <- c("Full", "Male", "Female")

# loop through data types and subsamples

cat("***Starting PLSR***\n")

for (data in data_types) {
  
  for (sample in subsamples) {
    
    cat(paste0("\n", "Data-Type: ", data, "; Subsample: ", sample, 
               "; Main Condition: ", main_cond, "\n"))
    
    if (data == "MMSE") {
      ## MMSE data
      # load data set
      mmse_datapath <- "Data/MMSEData/MMSE_data_full.csv"
      
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
      
      mmse_data_wide <- merge(mmse_data_wide, lab_data, by = "ID", all.x = TRUE)
      
      mmse_data_centered <- mmse_data_wide %>%
        group_by(Lab) %>%
        mutate(across(
          .cols = where(is.numeric) & !all_of("Age"),
          .fns = ~ . - mean(., na.rm = TRUE)
        )) %>%
        ungroup()
      
      # Hypothesis 3 - Aggregate Score
      aggregate_score <- rowMeans(cbind(
        mmse_data_centered$`avg_entropy_FL-FR_first_run_eyes_open` * 0.133,
        mmse_data_centered$avg_entropy_FL_first_run_eyes_open * 0.12,
        mmse_data_centered$avg_entropy_P_first_run_eyes_open * 0.113,
        mmse_data_centered$`max_slope_PL-PR_first_run_eyes_open` * 0.101,
        mmse_data_centered$`max_slope_FR-PR_first_run_eyes_open` * 0.083,
        mmse_data_centered$`avg_entropy_FL-PL_first_run_eyes_open` * 0.082,
        mmse_data_centered$auc_FL_first_run_eyes_open * 0.076,
        mmse_data_centered$`auc_FL-PL_first_run_eyes_open` * 0.07,
        mmse_data_centered$`auc_FR-PR_first_run_eyes_open` * 0.053,
        mmse_data_centered$max_slope_PR_first_run_eyes_open * -0.139,
        mmse_data_centered$`max_slope_ML-MR_first_run_eyes_open` * -0.134,
        mmse_data_centered$`avg_entropy_F-P_first_run_eyes_open` * -0.134,
        mmse_data_centered$max_slope_ML_first_run_eyes_open * -0.102,
        mmse_data_centered$`auc_F-P_first_run_eyes_open` * -0.084))
      
      hyp3 <- cor.test(aggregate_score, mmse_data_centered$gf_score, 
                       method = "pearson")
      
      print(hyp3)
      
      # Select sample data
      if (sample != "Full") {
        
        mmse_data_centered_subsample <- mmse_data_centered %>%
          filter(Gender == tolower(sample))
        
        X_AUC <- as.matrix(mmse_data_centered_subsample %>% 
                             select(contains("auc")))
        
        X_AVG <- as.matrix(mmse_data_centered_subsample %>% 
                             select(contains("avg_entropy")))
        
        X_MaxSlope <- as.matrix(mmse_data_centered_subsample %>% 
                                  select(contains("max_slope")))
        
        X <- cbind(X_AUC,X_AVG,X_MaxSlope)
        
        y <- mmse_data_centered_subsample$gf_score
        
      } else {
        
        mmse_data_centered_subsample <- mmse_data_centered
        
        X_AUC <- as.matrix(mmse_data_centered_subsample %>% 
                             select(contains("auc")))
        
        X_AVG <- as.matrix(mmse_data_centered_subsample %>% 
                             select(contains("avg_entropy")))
        
        X_MaxSlope <- as.matrix(mmse_data_centered_subsample %>% 
                                  select(contains("max_slope")))
        
        X <- cbind(X_AUC,X_AVG,X_MaxSlope)
        
        y <- mmse_data_centered_subsample$gf_score
        
      } # end for if sample != "Full"
      
      # cases used for training
      PLSmodel <- list(y,X)
      names(PLSmodel) <- list('dep','indep')
      
      cat(paste0("N = ", nrow(X)), "\n")
      
      # maximum number of components for PLS models
      max_comp <- 30
      
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
        
        tiff(paste(savepath, "PLSR_centered_selected_comps_", main_cond, "_", 
                   data, "_", sample, ".tiff", sep = ""), width = 2400, 
             height = 1800, res = 600)
        selectNcomp(PLSR_CV, method = "randomization", alpha=0.05, nperm=10000, 
                    plot = TRUE, ylim = c(0.9, 1.5))
        dev.off()
        
      } else {
        
        cat("Abs. Min: ", abs_min, "\n",
            "Sign. Component: ", numberSignifComp, sep = "")
        
        tiff(paste(savepath, "PLSR_centered_selected_comps_", main_cond, "_", 
                   data, "_", sample, ".tiff", sep = ""), width = 2400, 
             height = 1800, res = 600)
        selectNcomp(PLSR_CV, method = "randomization", alpha=0.05, nperm=10000, 
                    plot = TRUE, ylim = c(0.9, 1.5))
        dev.off()
        
        # Select component
        selected_comp <- max(numberSignifComp, abs_min)
        
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
        
        # significant regression coefficients obtained from bootstrapped 
        # confidence interval estimation with alpha=0.05
        filter(x,x$"2.5%">0)
        filter(x,x$"97.5%"<0)
        
        print(filter(x, `0.8333333%` > 0 | `99.16667%` < 0))
        
      } # end for if numberSignifComp == 0
      
      # Clear variables
      rm(X_AUC, X_AVG, X_MaxSlope, X, y, numberSignifComp, PLSmodel, PLSR_CV)
      
    } else {
      ## Microstate data
      microstate_datapath <- "Data/Microstates/Microstates_data_full.csv"
      
      # Check if file exists
      if (!file.exists(microstate_datapath)) {
        
        cat("No microstate data found...\n")
        break
        
      }
      
      microstate_data <- read_csv(microstate_datapath, show_col_types = FALSE)
      
      # get rid of unnecessary information
      microstate_data <- microstate_data %>%
        filter(
          Length == 40,
          Condition == main_cond,
          !is.na(Condition) & Condition != "<undefined>",
          rowSums(!is.na(across(-c(ID, gf_score, Gender, Age, Condition)))) > 0
        )
      
      microstate_data <- merge(microstate_data, lab_data, by = "ID", all.x = TRUE)
      
      microstate_data_centered <- microstate_data %>%
        group_by(Lab) %>%
        mutate(across(
          .cols = where(is.numeric) & !all_of("Age"),
          .fns = ~ . - mean(., na.rm = TRUE)
        )) %>%
        ungroup()
      
      # Select sample data
      if (sample != "Full") {
        
        microstate_data_centered_subsample <- microstate_data_centered %>%
          filter(Gender == tolower(sample))
        
        X_n_gfp_peaks <- as.matrix(microstate_data_centered_subsample %>%
                                     select(contains("n_gfp_peaks")))
        
        X_coverage <- as.matrix(microstate_data_centered_subsample %>% 
                                  select(contains("coverage")))
        
        X_lifespan <- as.matrix(microstate_data_centered_subsample %>%
                                  select(contains("lifespan") 
                                         & !contains("lifespan_peaks")))
        
        X_lifespan_peaks <- as.matrix(microstate_data_centered_subsample %>% 
                                        select(contains("lifespan_peaks")))
        
        X_frequency <- as.matrix(microstate_data_centered_subsample %>% 
                                   select(contains("frequence")))
        
        X_trans_prob <- as.matrix(microstate_data_centered_subsample %>%
                                    select(contains("transition_probability") 
                                           & !contains("transition_probability_peaks")))
        
        X_trans_prob_peaks <- as.matrix(microstate_data_centered_subsample %>% 
                                          select(contains("transition_probability_peaks")))
        
        X <- cbind(X_n_gfp_peaks, X_coverage, X_lifespan, X_lifespan_peaks,
                   X_frequency, X_trans_prob, X_trans_prob_peaks)
        
        y <- microstate_data_centered_subsample$gf_score
        
      } else {
        
        microstate_data_centered_subsample <- microstate_data_centered
        
        X_n_gfp_peaks <- as.matrix(microstate_data_centered_subsample %>%
                                     select(contains("n_gfp_peaks")))
        
        X_coverage <- as.matrix(microstate_data_centered_subsample %>% 
                                  select(contains("coverage")))
        
        X_lifespan <- as.matrix(microstate_data_centered_subsample %>%
                                  select(contains("lifespan") 
                                         & !contains("lifespan_peaks")))
        
        X_lifespan_peaks <- as.matrix(microstate_data_centered_subsample %>% 
                                        select(contains("lifespan_peaks")))
        
        X_frequency <- as.matrix(microstate_data_centered_subsample %>% 
                                   select(contains("frequence")))
        
        X_trans_prob <- as.matrix(microstate_data_centered_subsample %>%
                                    select(contains("transition_probability") 
                                           & !contains("transition_probability_peaks")))
        
        X_trans_prob_peaks <- as.matrix(microstate_data_centered_subsample %>% 
                                          select(contains("transition_probability_peaks")))
        
        X <- cbind(X_n_gfp_peaks, X_coverage, X_lifespan, X_lifespan_peaks,
                   X_frequency, X_trans_prob, X_trans_prob_peaks)
        
        y <- microstate_data_centered_subsample$gf_score
        
      } # end for if sample != "Full"
      
      # cases used for training
      PLSmodel <- list(y,X)
      names(PLSmodel) <- list('dep','indep')
      
      cat(paste0("N = ", nrow(X)), "\n")
      
      # maximum number of components for PLS models
      max_comp<-30
      
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
        
        tiff(paste(savepath, "PLSR_Centered_selected_comps_", main_cond, "_", 
                   data, "_", sample, ".tiff", sep = ""), width = 2400,
             height = 1800, res = 600)
        selectNcomp(PLSR_CV, method = "randomization", alpha=0.05, nperm=10000, 
                    plot = TRUE, ylim = c(0.9, 1.5))
        dev.off()
        
      } else {
        
        cat("Abs. Min: ", abs_min, "\n",
            "Sign. Component: ", numberSignifComp, sep = "")
        
        tiff(paste(savepath, "PLSR_Centered_selected_comps_", main_cond, "_", 
                   data, "_", sample, ".tiff", sep = ""), width = 2400, 
             height = 1800, res = 600)
        selectNcomp(PLSR_CV, method = "randomization", alpha=0.05, nperm=10000, 
                    plot = TRUE, ylim = c(0.9, 1.5))
        dev.off()
        
        # Select component
        selected_comp <- max(numberSignifComp, abs_min)
        
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
        cfs_int <- coefficients.boots(PLSR_BOOT, ncomp = 1, conf = (1 - (0.05 / 3)))
        
        # all intervals
        x <- cfs_int[[1]]
        
        # significant regression coefficients obtained from bootstrapped 
        # confidence interval estimation with alpha=0.05
        print(filter(x, `0.8333333%` > 0 | `99.16667%` < 0))
        
      } # end for if numberSignifComp == 0
      
    } # end for if data == "MMSE"
    
  } # end for subsample
  
} # end for data type


# PLSPM Model

# MMSE

# Define inner model
auc <- c(0,0,0,0)
max_slope <- c(0,0,0,0)
avg_ent <- c(0,0,0,0)
gf <- c(1,1,1,0)

mmse_path <- rbind(auc,max_slope, avg_ent, gf)
colnames(mmse_path) <- rownames(mmse_path)

innerplot(mmse_path)

# Define outer model
mmse_blocks <- list(
  colnames(mmse_data_centered_subsample %>% select(contains("auc"))), 
  colnames(mmse_data_centered_subsample %>% select(contains("max_slope"))), 
  colnames(mmse_data_centered_subsample %>% select(contains("avg_entropy"))), 
  "gf_score")

# Define mode
mmse_mode <- rep("A",4)

# Filter for Gender and run model separately
mmse_data_centered_plspm_female <-  mmse_data_centered %>%
  filter(Gender == "female") %>%
  drop_na(all_of(unlist(mmse_blocks)))

mmse_pls_female <- plspm(mmse_data_centered_plspm_female, mmse_path, mmse_blocks, modes = mmse_mode)

mmse_data_centered_plspm_male <-  mmse_data_centered %>%
  filter(Gender == "male") %>%
  drop_na(all_of(unlist(mmse_blocks)))

mmse_pls_male <- plspm(mmse_data_centered_plspm_male, mmse_path, mmse_blocks, modes = mmse_mode)

# Run Group Model
mmse_data_centered_plspm <-  mmse_data_centered %>%
  drop_na(all_of(unlist(mmse_blocks))) %>%
  filter(Gender != "<undefined>")

mmse_data_centered_plspm$Gender <- as.factor(mmse_data_centered_plspm$Gender)

mmse_pls <- plspm(mmse_data_centered_plspm, mmse_path, mmse_blocks, modes = mmse_mode)

mmse_plspm_boot <- plspm.groups(mmse_pls, mmse_data_centered_plspm$Gender, 
                                method = "bootstrap")

mmse_plspm_boot

bar_colors <- c("female" = "#FEB24C", "male" = "#74A9CF")

bar_mmse <- mmse_plspm_boot$test[, 2:3] %>%
  rownames_to_column("Path") %>%
  pivot_longer(-Path, names_to = "Group", values_to = "Coefficient") %>%
  mutate(Group = case_when(
    Group == "group.female" ~ "female",
    Group == "group.male" ~ "male",
    TRUE ~ Group
  ))

bar_mmse$Path <- sub("->gf", "", bar_mmse$Path)
bar_mmse$Path <- sub("_", ". ", bar_mmse$Path)
bar_mmse$Path <- sub("ent", "entropy", bar_mmse$Path)
bar_mmse$Path[1:2] <- toupper(bar_mmse$Path[1:2])
bar_mmse$Path[3:6] <- tools::toTitleCase(bar_mmse$Path[3:6])

plspm_plot_mmse <- ggplot(bar_mmse, aes(x = Path, y = Coefficient, 
                                        fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), 
           width = 0.7) +
  geom_hline(yintercept = 0, color = "gray50") +
  scale_fill_manual(values = bar_colors, name = NULL) +
  scale_y_continuous() +
  labs(x = NULL, y = "Path Coefficient") +
  theme_classic(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey90"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(color = "gray30"),
    legend.position = "right",
    legend.text = element_text(color = "gray40")
  )

plspm_plot_mmse

ggsave(filename = paste0(savepath, "PLSPM_Centered_MMSE_" , main_cond, ".tiff"), 
       plot = plspm_plot_mmse, width = 8, height = 5, dpi = 600)


# Microstates

# Define inner model
n_peaks <- c(0,0,0,0,0,0,0,0)
coverage <- c(0,0,0,0,0,0,0,0)
lifespan <- c(0,0,0,0,0,0,0,0)
lifespan_peaks <- c(0,0,0,0,0,0,0,0)
frequency <- c(0,0,0,0,0,0,0,0)
transition_probability <- c(0,0,0,0,0,0,0,0)
transition_probability_peaks <- c(0,0,0,0,0,0,0,0)
gf <- c(1,1,1,1,1,1,1,0)

microstate_path <- rbind(n_peaks, coverage, lifespan, lifespan_peaks, frequency,
                         transition_probability, transition_probability_peaks ,gf)
colnames(microstate_path) <- rownames(microstate_path)

innerplot(microstate_path)

# Define outer model
microstate_blocks <- list(
  colnames(microstate_data_centered %>% select(contains("n_gfp_peaks"))), 
  colnames(microstate_data_centered %>% select(contains("coverage"))), 
  colnames(microstate_data_centered %>% select(contains("lifespan"))),
  colnames(microstate_data_centered %>% select(contains("lifespan_peaks"))),
  colnames(microstate_data_centered %>% select(contains("frequence"))),
  colnames(microstate_data_centered %>% select(contains("transition_probability"))),
  colnames(microstate_data_centered %>% select(contains("transition_probability_peaks"))), 
  "gf_score")

# Define mode
microstate_mode <- rep("A",8)

# Filter for Gender and run model separately
microstate_data_centered_plspm_female <-  microstate_data_centered %>%
  filter(Gender == "female") %>%
  drop_na(all_of(unlist(microstate_blocks)))

microstate_pls_female <- plspm(microstate_data_centered_plspm_female, microstate_path, 
                               microstate_blocks, modes = microstate_mode)

microstate_data_centered_plspm_male <-  microstate_data_centered %>%
  filter(Gender == "male") %>%
  drop_na(all_of(unlist(microstate_blocks)))

mmse_pls_male <- plspm(microstate_data_centered_plspm_male, microstate_path, 
                       microstate_blocks, modes = microstate_mode)

# Run Group Model
microstate_data_centered_plspm <-  microstate_data_centered %>%
  drop_na(all_of(unlist(microstate_blocks))) %>%
  filter(Gender != "<undefined>")

microstate_data_centered_plspm$Gender <- as.factor(microstate_data_centered_plspm$Gender)

microstate_pls <- plspm(microstate_data_centered_plspm, microstate_path, 
                        microstate_blocks, modes = microstate_mode)

microstate_plspm_boot <- plspm.groups(microstate_pls, microstate_data_centered_plspm$Gender, 
                                      method = "bootstrap")

microstate_plspm_boot

bar_colors <- c("female" = "#FEB24C", "male" = "#74A9CF")

bar_ms <- microstate_plspm_boot$test[, 2:3] %>%
  rownames_to_column("Path") %>%
  pivot_longer(-Path, names_to = "Group", values_to = "Coefficient") %>%
  mutate(Group = case_when(
    Group == "group.female" ~ "female",
    Group == "group.male" ~ "male",
    TRUE ~ Group
  ))

bar_ms$Path <- sub("->gf", "", bar_ms$Path)
bar_ms$Path <- sub("_", " ", bar_ms$Path)
bar_ms$Path <- sub("_", " ", bar_ms$Path)
bar_ms$Path <- tools::toTitleCase(bar_ms$Path)
bar_ms$Path[1:2] <- sub("n", "Number of GFP", bar_ms$Path[1:2])

plspm_plot_microstate <- ggplot(bar_ms, aes(x = Path, y = Coefficient, 
                                            fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), 
           width = 0.7) +
  geom_hline(yintercept = 0, color = "gray50") +
  scale_fill_manual(values = bar_colors, name = NULL) +
  scale_y_continuous() +
  labs(x = NULL, y = "Path Coefficient") +
  theme_classic(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey90"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(color = "gray30"),
    legend.position = "right",
    legend.text = element_text(color = "gray40")
  )

plspm_plot_microstate

ggsave(filename = paste0(savepath, "PLSPM_Centered_Microstates_" , main_cond, 
                         ".tiff"), plot = plspm_plot_microstate, width = 8, 
       height = 5, dpi = 600)

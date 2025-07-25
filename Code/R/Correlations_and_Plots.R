# This script calculates re-test correlations for the MMSE and Microstate 
# features and plots them.
#
# This script was created by Christoph Fruehlinger (June 2025)
# Last edit: July 2025

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, rlist, cowplot, ggtext)

options(scipen = 999)

rm(list = ls())

source("./Code/R/filter_plots.R")

# Create output folder
savepath = "Results/"

if (!dir.exists(savepath)) {
  dir.create(savepath)
}

# List of analysis conditions
data_types = c("MMSE", "Microstates")

# loop through data types

cat("***Starting Re-Test Correlations***\n")

for (data in data_types) {
  
  cat(paste0("\n", "Data-Type: ", data, "\n"))
  
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

    mmse_data <- mmse_data %>%
      filter(
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
    
    
    # Set conditions and sets to compare
    cond1 <- "first_run_eyes_open"
    otherconds <- c("first_run_eyes_closed", "second_run_eyes_open", 
                    "third_run_eyes_open")
    condnames <- c("run 1 EC", "run 2 EO", "run 3 EO")
    sets <- c("C", "F", "FL", "FR", "ML", "MR", "P", "PL", "PR")
    feat_vars <- c("auc", "max_slope", "avg_entropy")
    featurenames = c("AUC", "MaxSlope", "AvgEnt")
    
    # Collect correlations
    p_vals <- list()
    plot_list <- list()
    correlations_mmse <- tibble()
    
    # loop over conditions and sets
    for (cond2 in otherconds) {
      
      for (target_set in sets) {
        
        # filter
        T1 <- mmse_data %>%
          filter(Condition == cond1, Set == target_set)
        
        T2 <- mmse_data %>%
          filter(Condition == cond2, Set == target_set)
        
        # sort IDs
        common_ids <- intersect(T1$ID, T2$ID)
        T1 <- T1 %>% 
          filter(ID %in% common_ids) %>% 
          arrange(ID)
        T2 <- T2 %>% 
          filter(ID %in% common_ids) %>% 
          arrange(ID)
        
        if (nrow(T1) == 0 || nrow(T2) == 0) next
        
        stopifnot(identical(T1$ID, T2$ID))
        
        # calculate correlations
        rho_vec <- numeric(length(feat_vars))
        p_vec <- numeric(length(feat_vars))
        
        for (i in seq_along(feat_vars)) {
          x <- T1[[feat_vars[i]]]
          y <- T2[[feat_vars[i]]]
          
          test <- cor.test(x, y, method = "spearman", exact = FALSE)
          rho_vec[i] <- test$estimate
          p_vec[i] <- test$p.value
        }
        
        p_vals[[paste(cond2, target_set, sep = "_")]] <- p_vec
        
        # save results
        corr_tbl <- tibble(
          Set = target_set,
          Condition1 = cond1,
          Condition2 = cond2,
          Feature = feat_vars,
          SpearmanRho = rho_vec,
          pValue = p_vec
        )
        
        correlations_mmse <- bind_rows(
          correlations_mmse,
          corr_tbl %>% 
            mutate(Comparison = paste(cond2, target_set, sep = "_"))
        )
        
        correlations_mmse <- correlations_mmse %>%
          mutate(pValue_adj = p.adjust(pValue, method = "holm"))
        
        # Scatter plots
        
        for (i in seq_along(feat_vars)) {
          
          feature <- feat_vars[i]
          featurename <- featurenames[i]
          plot_id <- paste(cond2, target_set, feature, sep = "_")
          
          if (p_vec[i] < 0.001) {
            
            p_disp <- '< .001'
            
          } else {
            
            p_disp <- paste0('= ', round(p_vec[i], 3))
            
          }
          
          cond_label <- condnames[which(otherconds == cond2)]
          
          plot_list[[plot_id]] <- ggplot(data = data.frame(x = T1[[feature]],
                                                           y = T2[[feature]]), 
                                         aes(x = x, y = y)) +
            geom_point(color = "#2c3e50", alpha = 0.7) +
            geom_smooth(method = "lm", formula = y ~ x, 
                        se = TRUE, color = "#e74c3c") +
            labs(
              title = paste0("Set: ", target_set, "<br>",
                             "*r* = ", round(rho_vec[i], 2), ", *p* ", p_disp),
              x = "run 1 EO",
              y = cond_label
            ) +
            theme_classic() +
            theme(
              panel.grid.major = element_line(color = "grey80"),
              panel.grid.minor = element_line(color = "grey90"),
              panel.grid.major.x = element_line(),
              panel.grid.major.y = element_line(),
              plot.title = element_markdown(size = 14))
        }
      }
    }
    
    # Descriptive
    correlations_mmse$Set <- as.factor(correlations_mmse$Set)
    correlations_mmse$Condition1 <- as.factor(correlations_mmse$Condition1)
    correlations_mmse$Condition2 <- as.factor(correlations_mmse$Condition2)
    correlations_mmse$Feature <- as.factor(correlations_mmse$Feature)
    
    write_csv(correlations_mmse, 
              paste0(savepath, "retest_correlations_", data, '.csv'))
    
    for (cond2 in otherconds) {
    
      heatmap_data <- correlations_mmse %>%
        filter(Condition1 == cond1,
               Condition2 == cond2) %>%
        mutate(Feature = recode(Feature,
                                auc = 'AUC',
                                max_slope = 'Max. Slope',
                                avg_entropy =  'Avg. Entropy'))
      
      heatmap_data <- heatmap_data %>%
        mutate(signif = case_when(
          pValue < 0.001 ~ "***",
          pValue < 0.01  ~ "**",
          pValue < 0.05  ~ "*",
          TRUE           ~ ""
        ))
      

      heatmap_plot <- ggplot(heatmap_data, aes(x = Feature, y = Set, 
                                               fill = SpearmanRho)) +
        geom_tile(color = "white") +
        geom_text(aes(label = signif), color = "white", size = 6) +
        theme_minimal(base_size = 14) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.title = element_blank(),
              panel.grid = element_blank(),
              plot.background = element_rect(fill = "white", color = NA)) +
        scale_fill_viridis_c(name = "Correlation")
      
      heatmap_filename <- paste0(savepath, cond2, '_', '_hm_mmse_retest.tiff')
      
      ggsave(filename = heatmap_filename, plot = heatmap_plot, width = 8, 
             height = 4, dpi = 600)
      
    }
    
    # Plot vector for indexing
    plot_vector = c(1,4,7,10,13,16,19,22,25)
    
    # Create plot grids
    auc_plot_1EC <- plot_grid(plotlist = plot_list[plot_vector], 
                              labels = "AUTO", ncol = 3)
    max_slope_plot_1EC <- plot_grid(plotlist = plot_list[plot_vector + 1], 
                                    labels = "AUTO", ncol = 3)
    avg_ent_plot_1EC <- plot_grid(plotlist = plot_list[plot_vector + 2], 
                                  labels = "AUTO", ncol = 3)
    
    auc_plot_2EO <- plot_grid(plotlist = plot_list[plot_vector + 27], 
                              labels = "AUTO", ncol = 3)
    max_slope_plot_2EO <- plot_grid(plotlist = plot_list[plot_vector + 28], 
                                    labels = "AUTO", ncol = 3)
    avg_ent_plot_2EO <- plot_grid(plotlist = plot_list[plot_vector + 29], 
                                  labels = "AUTO", ncol = 3)
    
    auc_plot_3EO <- plot_grid(plotlist = plot_list[plot_vector + 54], 
                              labels = "AUTO", ncol = 3)
    max_slope_plot_3EO <- plot_grid(plotlist = plot_list[plot_vector + 55], 
                                    labels = "AUTO", ncol = 3)
    avg_ent_plot_3EO <- plot_grid(plotlist = plot_list[plot_vector + 56], 
                                  labels = "AUTO", ncol = 3)
    
    # save plots
    ggsave(filename = paste0(savepath, "auc_1EC.tiff"), plot = auc_plot_1EC, 
           width = 7, height = 7, dpi = 600)
    ggsave(filename = paste0(savepath, "auc_2EO.tiff"), plot = auc_plot_2EO, 
           width = 7, height = 7, dpi = 600)
    ggsave(filename = paste0(savepath, "auc_3EO.tiff"), plot = auc_plot_3EO, 
           width = 7, height = 7, dpi = 600)
    
    ggsave(filename = paste0(savepath, "max_slope_1EC.tiff"), 
           plot = max_slope_plot_1EC, width = 7, height = 7, dpi = 600)
    ggsave(filename = paste0(savepath, "max_slope_2EO.tiff"), 
           plot = max_slope_plot_2EO, width = 7, height = 7, dpi = 600)
    ggsave(filename = paste0(savepath, "max_slope_3EO.tiff"), 
           plot = max_slope_plot_3EO, width = 7, height = 7, dpi = 600)
    
    ggsave(filename = paste0(savepath, "avg_ent_1EC.tiff"), 
           plot = avg_ent_plot_1EC, width = 7, height = 7, dpi = 600)
    ggsave(filename = paste0(savepath, "avg_ent_2EO.tiff"), 
           plot = avg_ent_plot_2EO, width = 7, height = 7, dpi = 600)
    ggsave(filename = paste0(savepath, "avg_ent_3EO.tiff"), 
           plot = avg_ent_plot_3EO, width = 7, height = 7, dpi = 600)
    
    ## MMSE Vector Plots
    
    # Convert to long format
    mmse_long <- mmse_data %>%
      pivot_longer(cols = starts_with("mmse_"), names_to = "mmse", 
                   values_to = "value")
    
    # Summarize data
    summary_mmse_gender <- mmse_long %>%
      filter(Gender != "<undefined>") %>%
      group_by(Gender, Set, mmse, Condition) %>%
      summarise(
        mean = mean(value, na.rm = TRUE),
        sd = sd(value, na.rm = TRUE),
        n = sum(!is.na(value)),
        se = sd / sqrt(n),
        ci_lower = mean - 1.96 * se,
        ci_upper = mean + 1.96 * se,
        .groups = "drop"
      ) %>%
      drop_na()
    
    # Rename columns for plotting
    summary_mmse_gender$mmse <- factor(summary_mmse_gender$mmse, 
                                       levels = paste0("mmse_", 1:12))
    
    # Repeat for full sample
    summary_mmse <- mmse_long %>%
      group_by(Set, mmse, Condition) %>%
      summarise(
        mean = mean(value, na.rm = TRUE),
        sd = sd(value, na.rm = TRUE),
        n = sum(!is.na(value)),
        se = sd / sqrt(n),
        ci_lower = mean - 1.96 * se,
        ci_upper = mean + 1.96 * se,
        .groups = "drop"
      )
    
    summary_mmse$mmse <- factor(summary_mmse$mmse, 
                                levels = paste0("mmse_", 1:12))
    
    # Set gender to "full" in full sample
    summary_mmse_full <- summary_mmse %>%
      mutate(Gender = "full") %>%
      drop_na()
    
    # Combine and plot
    summary_all <- bind_rows(summary_mmse_gender, summary_mmse_full)
    
    summary_all$Gender <- factor(summary_all$Gender, 
                                 levels = c("full", "female", "male"))
    
    # Plot for Main Condition
    vector_plot <- summary_all %>%
      filter(Condition == "first_run_eyes_open") %>%
      ggplot(aes(x = mmse, y = mean, color = Set, group = Set)) +
      geom_line(size = 0.5) +
      geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), 
                    width = 0.5, size = 0.4)  +
      geom_point(shape = 21, colour = "black", stroke = 0.25, size = 0.5) +
      scale_x_discrete(labels = as.character(1:12)) +
      labs(x = "Timescale Factor", y = "mMSE") +
      theme_bw() +
      theme(
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14, face = "bold"),
        legend.position = "right",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 9)) +
      scale_colour_viridis_d() +
      facet_grid(. ~ Gender, labeller = labeller(Gender = c(
        "full" = "Full Sample",
        "female" = "Female",
        "male" = "Male"
      )))
    
    ggsave(filename = paste0(savepath, "MMSE_vectors.tiff"), 
           plot = vector_plot, width = 8, height = 3, dpi = 600)
    
    # Plot for All Condition
    vector_plot_all_conds <- summary_all %>%
      ggplot(aes(x = mmse, y = mean, color = Set, group = Set)) +
      geom_line(size = 0.5) +
      geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), 
                    width = 0.5, size = 0.4)  +
      geom_point(shape = 21, colour = "black", stroke = 0.25, size = 0.5) +
      scale_x_discrete(labels = as.character(1:12)) +
      labs(x = "Timescale Factor", y = "mMSE") +
      theme_bw() +
      theme(
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14, face = "bold"),
        legend.position = "right",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 9)) +
      scale_colour_viridis_d() +
      facet_grid(Condition ~ Gender, labeller = labeller(Gender = c(
        "full" = "Full Sample", "female" = "Female", "male" = "Male"), 
        Condition = c("first_run_eyes_closed"  = "Run 1 EC",
                      "first_run_eyes_open"    = "Run 1 EO",
                      "second_run_eyes_closed" = "Run 2 EC",
                      "second_run_eyes_open"   = "Run 2 EO",
                      "third_run_eyes_closed"  = "Run 3 EC",
                      "third_run_eyes_open"  = "Run 3 EO")))
    
    ggsave(filename = paste0(savepath, "MMSE_vectors_all_conds.tiff"), 
           plot = vector_plot_all_conds, width = 8, height = 8, dpi = 600)
    
  } else {
    
    ## Microstate data
    microstate_datapath = "Data/Microstates/Microstates_data_full.csv"
    
    # Check if file exists
    if (!file.exists(microstate_datapath)) {
      
      cat("No microstate data found...\n")
      break
      
    }
    
    microstate_data <- read_csv(microstate_datapath, show_col_types = FALSE)
    
    # Microstates = [A: 4, B: 3, C: 1, D: 2, F: 0]
    # Rename columns 
    microstate_names <- c("0" = "F", "1" = "C", "2" = "D", "3" = "B", "4" = "A")
    
    colnames(microstate_data) <- str_replace_all(colnames(microstate_data), 
                                                 microstate_names)
    
    microstate_data <- microstate_data %>%
      filter(
        !is.na(Condition),
        rowSums(!is.na(across(-c(ID, gf_score, Gender, Age, Condition)))) > 0
      )
    
    # Microstate Hypotheses 2 and 3
    # Full sample
    microstate_data_sim_gev <- microstate_data %>%
      select(ID, Gender, gf_score, Condition, 
             similarity_subject_micosates, gev_group) %>%
      filter(Condition == cond1)
    
    x <- microstate_data_sim_gev$similarity_subject_micosates
    y <- microstate_data_sim_gev$gf_score
    full_sim <- cor.test(x, y, method = "pearson")
    
    x <- microstate_data_sim_gev$gev_group
    full_gev <- cor.test(x, y, method = "pearson")
    
    # Female
    microstate_data_sim_gev_female <- microstate_data_sim_gev %>%
      filter(Gender == "female")
    
    x <- microstate_data_sim_gev_female$similarity_subject_micosates
    y <- microstate_data_sim_gev_female$gf_score
    female_sim <- cor.test(x, y, method = "pearson")
    
    x <- microstate_data_sim_gev_female$gev_group
    female_gev <- cor.test(x, y, method = "pearson")
    
    # Male
    microstate_data_sim_gev_male <- microstate_data_sim_gev %>%
      filter(Gender == "male")
    
    x <- microstate_data_sim_gev_male$similarity_subject_micosates
    y <- microstate_data_sim_gev_male$gf_score
    male_sim <- cor.test(x, y, method = "pearson")
    
    x <- microstate_data_sim_gev_male$gev_group
    male_gev <- cor.test(x, y, method = "pearson")
    
    hypothesis_2_3_correlations <- tibble(
      Sample = c("Full", "Full", "Female", "Female", "Male", "Male"),
      X = c("similarity_subject_micosates", "gev_group",
            "similarity_subject_micosates", "gev_group",
            "similarity_subject_micosates", "gev_group"),
      Y = rep("gf_score", 6),
      Pearson_R = c(full_sim$estimate, full_gev$estimate,
                    female_sim$estimate, female_gev$estimate,
                    male_sim$estimate, male_gev$estimate),
      p_value = c(full_sim$p.value, full_gev$p.value,
                  female_sim$p.value, female_gev$p.value,
                  male_sim$p.value, male_gev$p.value)
    )
    
    hypothesis_2_3_correlations <- hypothesis_2_3_correlations %>%
      mutate(p_value_adj = p.adjust(p_value, method = "holm"))
    
    # Set conditions and microstates to compare
    cond1 <- "first_run_eyes_open"
    otherconds <- c("first_run_eyes_closed", "second_run_eyes_open", 
                    "third_run_eyes_open")
    condnames <- c("run 1 EC", "run 2 EO", "run 3 EO")
    ms <- c("A", "B", "C", "D", "F")
    feat_vars_no_trans_prob <- c("coverage", "lifespan", 
                                 "lifespan_peaks", "frequence")
    feat_vars <- c("n_gfp_peaks", "coverage", "lifespan", "lifespan_peaks", 
                   "frequence","transition_probability", 
                   "transition_probability_peaks")
    featurenames = c("Number of GFP Peaks", "Coverage", "Lifespan", 
                     "Lifespan at GFP Peaks", "Frequency", "Transition Probability",
                     "Transition Probability at GFP Peaks")
    
    features <- unlist(lapply(feat_vars_no_trans_prob, function(f) paste0(f, "_", ms)))
    trans_prob <- names(microstate_data)[startsWith(names(microstate_data), 
                                                    "transition_probability_")]
    
    all_features <- c("n_gfp_peaks", features, trans_prob)
    
    # Collect correlations
    plot_list_ms <- list()
    correlations_microstates <- tibble()
    
    # loop over conditions and microstates
    for (cond2 in otherconds) {
      
      # filter
      T1 <- microstate_data %>%
        filter(Condition == cond1)
      
      T2 <- microstate_data %>%
        filter(Condition == cond2)
      
      # sort IDs
      common_ids <- intersect(T1$ID, T2$ID)
      T1 <- T1 %>% 
        filter(ID %in% common_ids) %>% 
        arrange(ID)
      T2 <- T2 %>% 
        filter(ID %in% common_ids) %>% 
        arrange(ID)
      
      if (nrow(T1) == 0 || nrow(T2) == 0) next
      
      stopifnot(identical(T1$ID, T2$ID))
      
      # Loop through microstate features
      for (feat in all_features) {
        
        if (!feat %in% names(T1) || !feat %in% names(T2)) next
        
        x <- T1[[feat]]
        y <- T2[[feat]]
        
        # calculate correlation
        test <- cor.test(x, y, method = "spearman", exact = FALSE)

        # Store data
        correlations_microstates <- bind_rows(
          correlations_microstates,
          tibble(
            Condition1 = cond1,
            Condition2 = cond2,
            Feature = if (startsWith(feat, "transition_probability")) {
              gsub("_[A-Z]_[A-Z]$", "", feat)
            } else {
              gsub("_[A-Z]$", "", feat)
            },
            Microstate = if (startsWith(feat, "transition_probability")) {
              sub(".*_([A-Z]_[A-Z])$", "\\1", feat)
            } else {
              sub(".*_([A-Z])$", "\\1", feat) 
            },
            SpearmanRho = test$estimate,
            pValue = test$p.value))
      }
    }
        
    # Correct for multiple comparisons
    correlations_microstates <- correlations_microstates %>%
      mutate(pValue_adj = p.adjust(pValue, method = "holm"))
    
    for (cond2 in otherconds) {
      
      T1 <- microstate_data %>%
        filter(Condition == cond1)
      
      T2 <- microstate_data %>%
        filter(Condition == cond2)
      
      # sort IDs
      common_ids <- intersect(T1$ID, T2$ID)
      T1 <- T1 %>% 
        filter(ID %in% common_ids) %>% 
        arrange(ID)
      T2 <- T2 %>% 
        filter(ID %in% common_ids) %>% 
        arrange(ID)
      
      if (nrow(T1) == 0 || nrow(T2) == 0) next
      
      stopifnot(identical(T1$ID, T2$ID))

      # Scatter plots
      for (feat in all_features) {
        
        # Define for display
        cond_label <- condnames[which(otherconds == cond2)]
        feature <- feat
        
        feat_base = if (startsWith(feat, "transition_probability")) {
          gsub("_[A-Z]_[A-Z]$", "", feat)
        } else {
          gsub("_[A-Z]$", "", feat)
        }
        
        microstate = if (startsWith(feat, "transition_probability")) {
          sub(".*_([A-Z]_[A-Z])$", "\\1", feat)
        } else {
          sub(".*_([A-Z])$", "\\1", feat) 
        }
        
        microstate_to_plot = if (startsWith(feat, "transition_probability")) {
          r <- sub(".*_([A-Z]_[A-Z])$", "\\1", feat)
          gsub("_", " to ", x = r)
        } else if (feat_base == "n_gfp_peaks") { 
          ""
        } else {
          sub(".*_([A-Z])$", "\\1", feat) 
        }
        
        feat_name <- featurenames[which(feat_vars == feat_base)]
        plot_id <- paste(cond2, feat, sep = "_")
        
        # Get Rho and p-Values
        corr_row <- correlations_microstates %>%
          filter(Condition1 == cond1,
                 Condition2 == cond2,
                 Feature == feat_base,
                 Microstate == microstate)
        
        rho_curr <- corr_row$SpearmanRho
        p_val_curr <- corr_row$pValue_adj
        
        if (p_val_curr < 0.001) {
          
          p_disp <- '< .001'
          
        } else {
          
          p_disp <- paste0('= ', round(p_val_curr, 3))
          
        }
        
        # Plot
        plot_list_ms[[plot_id]] <- ggplot(data = data.frame(x = T1[[feat]],
                                                         y = T2[[feat]]), 
                                       aes(x = x, y = y)) +
          geom_point(color = "#2c3e50", alpha = 0.7) +
          geom_smooth(method = "lm", formula = y ~ x, 
                      se = TRUE, color = "#e74c3c") +
          labs(
            title = paste0(feat_name, " ", microstate_to_plot, "<br>",
                           "*r* = ", round(rho_curr, 2), ", *p* ", p_disp),
            x = "run 1 EO",
            y = cond_label
          ) +
          theme_classic() +
          theme(
            panel.grid.major = element_line(color = "grey80"),
            panel.grid.minor = element_line(color = "grey90"),
            panel.grid.major.x = element_line(),
            panel.grid.major.y = element_line(),
            plot.title = element_markdown(size = 14))
      }
    }
    
    # Descriptive
    correlations_microstates$Condition1 <- as.factor(correlations_microstates$Condition1)
    correlations_microstates$Condition2 <- as.factor(correlations_microstates$Condition2)
    correlations_microstates$Feature <- as.factor(correlations_microstates$Feature)
    correlations_microstates$Microstate <- as.factor(correlations_microstates$Microstate)
    
    write_csv(correlations_microstates, 
              paste0(savepath, "retest_correlations_", data, '.csv'))
    
    # Create Plots for Re-Test for each Microstate and Condition
    plot_names <- names(plot_list_ms)
    
    for (cond2 in otherconds) {

      for (target_ms in ms) {
        
        selected_plots <- filter_plots(plot_names, cond2, target_ms)
        
        subplots <- plot_grid(plotlist = plot_list_ms[selected_plots], 
                                    labels = "AUTO", ncol = 2)
        
        subplot_filename <- paste0(savepath, target_ms, '_', 
                                   cond2, '_retest.tiff')
          
        ggsave(filename = subplot_filename, plot = subplots, width = 10,
               height = 13, dpi = 600)
        
      }
      
      heatmap_data <- correlations_microstates %>%
        filter(Condition1 == cond1,
               Condition2 == cond2,
               !Feature %in% c("n_gfp_peaks", "transition_probability", 
                               "transition_probability_peaks")) %>%
        mutate(Feature = recode(Feature, 
                                coverage = "Coverage",
                                lifespan = "Lifespan",
                                lifespan_peaks  =  "Lifespan (GFP Peaks)",
                                frequence = "Frequency"))
      
      heatmap_data <- heatmap_data %>%
        mutate(signif = case_when(
          pValue_adj < 0.001 ~ "***",
          pValue_adj < 0.01  ~ "**",
          pValue_adj < 0.05  ~ "*",
          TRUE               ~ ""
        ))
      
      heatmap_data$Microstate <- factor(heatmap_data$Microstate, 
                                   levels = c("F", "D", "C", "B", "A"))
      
      heatmap_plot <- ggplot(heatmap_data, aes(x = Feature, y = Microstate, 
                                               fill = SpearmanRho)) +
        geom_tile(color = "white") +
        geom_text(aes(label = signif), color = "white", size = 6) +
        theme_minimal(base_size = 14) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.title = element_blank(),
              panel.grid = element_blank(),
              plot.background = element_rect(fill = "white", color = NA)) +
        scale_fill_viridis_c(name = "Correlation")
      
      heatmap_filename <- paste0(savepath, cond2, '_', '_hm_microstate_retest.tiff')
      
      ggsave(filename = heatmap_filename, plot = heatmap_plot, width = 8, 
             height = 4, dpi = 600)
      
    }
    
    # Last plot for Number of GFP Peaks
    n_gfp_peaks_plot <- plot_grid(plotlist = plot_list_ms[c(1, 72, 143)], 
                                  labels = "AUTO", ncol = 3)
    
    subplot_filename <- paste0(savepath, 'n_gfp_peaks_retest.tiff')
    
    ggsave(filename = subplot_filename, plot = n_gfp_peaks_plot, width = 8,
           height = 3, dpi = 600)
    
  } # end for if data == "MMSE"
  
} # end for data type

cat("Done")

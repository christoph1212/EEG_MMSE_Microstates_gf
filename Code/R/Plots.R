# This script calculates re-test correlations for the mMSE features and plots
# them.
#
# This script was created by Christoph Fruehlinger (June 2025)
# Last edit: June 2025

if (!require("pacman")) install.packages("pacman")
pacman::p_load(pls, mvdalab, tidyverse, haven, rlist, cowplot, ggtext)


options(max.print=1000000)

rm(list = ls())

# load data set
mmse_datapath = "Data/MMSEData/MMSE_data_full.csv"
mmse_data <- read_csv(mmse_datapath)

mmse_data <- mmse_data %>%
  filter(
    !is.na(Condition) & Condition != "<undefined>",
    !is.na(Set) & Set != "<undefined>",
    rowSums(!is.na(across(-c(ID, gf_score, Gender, Age, Condition, Set)))) > 0
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
otherconds <- c("first_run_eyes_closed", "second_run_eyes_open", "third_run_eyes_open")
condnames <- c("run 1 EC", "run 2 EO", "run 3 EO")
sets <- c("C", "F", "FL", "FR", "ML", "MR", "P", "PL", "PR")
feat_vars <- c("auc", "max_slope", "avg_entropy")
featurenames = c("AUC", "MaxSlope", "AvgEnt")

# Collect correlations
p_vals <- list()
plot_list <- list()
correlations <- tibble()

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
    T1 <- T1 %>% filter(ID %in% common_ids) %>% arrange(ID)
    T2 <- T2 %>% filter(ID %in% common_ids) %>% arrange(ID)
    
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
    
    correlations <- bind_rows(
      correlations,
      corr_tbl %>% 
        mutate(Comparison = paste(cond2, target_set, sep = "_"))
    )
    
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
      
      plot_list[[plot_id]] <- ggplot(data = data.frame(x = T1[[feature]], y = T2[[feature]]), aes(x = x, y = y)) +
        geom_point(color = "#2c3e50", alpha = 0.7) +
        geom_smooth(method = "lm", se = TRUE, color = "#e74c3c") +
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
correlations$Set <- as.factor(correlations$Set)
correlations$Condition1 <- as.factor(correlations$Condition1)
correlations$Condition2 <- as.factor(correlations$Condition2)
correlations$Feature <- as.factor(correlations$Feature)

psych::describeBy(correlations, group = "Feature")

# Plot vector for indexing
plot_vector = c(1,4,7,10,13,16,19,22,25)

# Create plot grids
auc_plot_1EC <- plot_grid(plotlist = plot_list[plot_vector], labels = "AUTO", ncol = 3)
max_slope_plot_1EC <- plot_grid(plotlist = plot_list[plot_vector + 1], labels = "AUTO", ncol = 3)
avg_ent_plot_1EC <- plot_grid(plotlist = plot_list[plot_vector + 2], labels = "AUTO", ncol = 3)

auc_plot_2EO <- plot_grid(plotlist = plot_list[plot_vector + 27], labels = "AUTO", ncol = 3)
max_slope_plot_2EO <- plot_grid(plotlist = plot_list[plot_vector + 28], labels = "AUTO", ncol = 3)
avg_ent_plot_2EO <- plot_grid(plotlist = plot_list[plot_vector + 29], labels = "AUTO", ncol = 3)

auc_plot_3EO <- plot_grid(plotlist = plot_list[plot_vector + 54], labels = "AUTO", ncol = 3)
max_slope_plot_3EO <- plot_grid(plotlist = plot_list[plot_vector + 55], labels = "AUTO", ncol = 3)
avg_ent_plot_3EO <- plot_grid(plotlist = plot_list[plot_vector + 56], labels = "AUTO", ncol = 3)

# save plots
savepath = "Results/"

if (!dir.exists(savepath)) {
  dir.create(savepath)
}

ggsave(filename = paste0(savepath, "auc_1EC.jpeg"), plot = auc_plot_1EC, dpi = 600)
ggsave(filename = paste0(savepath, "auc_2EO.jpeg"), plot = auc_plot_2EO, dpi = 600)
ggsave(filename = paste0(savepath, "auc_3EO.jpeg"), plot = auc_plot_3EO, dpi = 600)

ggsave(filename = paste0(savepath, "max_slope_1EC.jpeg"), plot = max_slope_plot_1EC, dpi = 600)
ggsave(filename = paste0(savepath, "max_slope_2EO.jpeg"), plot = max_slope_plot_2EO, dpi = 600)
ggsave(filename = paste0(savepath, "max_slope_3EO.jpeg"), plot = max_slope_plot_3EO, dpi = 600)

ggsave(filename = paste0(savepath, "avg_ent_1EC.jpeg"), plot = avg_ent_plot_1EC, dpi = 600)
ggsave(filename = paste0(savepath, "avg_ent_2EO.jpeg"), plot = avg_ent_plot_2EO, dpi = 600)
ggsave(filename = paste0(savepath, "avg_ent_3EO.jpeg"), plot = avg_ent_plot_3EO, dpi = 600)
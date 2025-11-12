# This script loads and analyses re-test correlations for  MMSE and Microstate 
# features created in Correlations_and_Plots.R.
#
# This script was created by Christoph Fruehlinger (November 2025)
# Last edit: November 2025

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse)

rm(list = ls())

## mMSE

mmse_retest <- read.csv("Results/retest_correlations_MMSE.csv", header = TRUE)

# No Difference Scores
mmse_retest_no_diff <- mmse_retest %>%
  filter(nchar(Set) <= 2)

mmse_retest_no_diff %>%
  group_by(Condition2, Feature) %>%
  summarise(
    mean_rho = mean(SpearmanRho, na.rm = TRUE),
    sd_rho = sd(SpearmanRho, na.rm = TRUE),
    min_rho = min(SpearmanRho),
    max_rho = max(SpearmanRho)
  )

mmse_retest_no_diff %>%
  group_by(Feature) %>%
  summarise(
    mean_rho = mean(SpearmanRho, na.rm = TRUE),
    sd_rho = sd(SpearmanRho, na.rm = TRUE),
    min_rho = min(SpearmanRho),
    max_rho = max(SpearmanRho)
  )

mmse_retest_no_diff %>%
  group_by(Set) %>%
  summarise(
    mean_rho = mean(SpearmanRho, na.rm = TRUE),
    sd_rho = sd(SpearmanRho, na.rm = TRUE),
    min_rho = min(SpearmanRho),
    max_rho = max(SpearmanRho)
  ) %>%
  arrange(desc(mean_rho)) %>% print(n = 50)

# Only Difference Scores

mmse_retest_diff <- mmse_retest %>%
  filter(nchar(Set) > 2)

mmse_retest_diff %>%
  group_by(Condition2, Feature) %>%
  summarise(
    mean_rho = round(mean(SpearmanRho, na.rm = TRUE), 2),
    sd_rho = round(sd(SpearmanRho, na.rm = TRUE), 2),
    min_rho = round(min(SpearmanRho), 2),
    max_rho = round(max(SpearmanRho), 2)
  )

mmse_retest_diff %>%
  group_by(Feature) %>%
  summarise(
    mean_rho = round(mean(SpearmanRho, na.rm = TRUE), 2),
    sd_rho = round(sd(SpearmanRho, na.rm = TRUE), 2),
    min_rho = round(min(SpearmanRho), 2),
    max_rho = round(max(SpearmanRho), 2)
  )

mmse_retest_diff %>%
  group_by(Set) %>%
  summarise(
    mean_rho = round(mean(SpearmanRho, na.rm = TRUE), 2),
    sd_rho = round(sd(SpearmanRho, na.rm = TRUE), 2),
    min_rho = round(min(SpearmanRho), 2),
    max_rho = round(max(SpearmanRho), 2)
  ) %>%
  arrange(desc(mean_rho)) %>% 
  print(n = 50)

# All

mmse_retest %>%
  group_by(Condition2, Feature) %>%
  summarise(
    mean_rho = round(mean(SpearmanRho, na.rm = TRUE), 2),
    sd_rho = round(sd(SpearmanRho, na.rm = TRUE), 2),
    min_rho = round(min(SpearmanRho), 2),
    max_rho = round(max(SpearmanRho), 2)
  )

mmse_retest %>%
  group_by(Feature) %>%
  summarise(
    mean_rho = round(mean(SpearmanRho, na.rm = TRUE), 2),
    sd_rho = round(sd(SpearmanRho, na.rm = TRUE), 2),
    min_rho = round(min(SpearmanRho), 2),
    max_rho = round(max(SpearmanRho), 2)
  )

mmse_retest %>%
  group_by(Set) %>%
  summarise(
    mean_rho = round(mean(SpearmanRho, na.rm = TRUE), 2),
    sd_rho = round(sd(SpearmanRho, na.rm = TRUE), 2),
    min_rho = round(min(SpearmanRho), 2),
    max_rho = round(max(SpearmanRho), 2)
  ) %>%
  arrange(desc(mean_rho)) %>% print(n = 50)


# Microstates
microstate_retest <- read.csv("Results/retest_correlations_Microstates.csv", header = TRUE)

microstate_retest %>%
  group_by(Condition2, Feature) %>%
  summarise(
    mean_rho = round(mean(SpearmanRho, na.rm = TRUE), 2),
    sd_rho = round(sd(SpearmanRho, na.rm = TRUE), 2),
    min_rho = round(min(SpearmanRho), 2),
    max_rho = round(max(SpearmanRho), 2),
  ) %>% 
  print(n = 50)

microstate_retest %>%
  group_by(Feature) %>%
  summarise(
    mean_rho = round(mean(SpearmanRho, na.rm = TRUE), 2),
    sd_rho = round(sd(SpearmanRho, na.rm = TRUE), 2),
    min_rho = round(min(SpearmanRho), 2),
    max_rho = round(max(SpearmanRho),  2),
  ) %>%
  arrange(desc(mean_rho)) %>% 
  print(n = 50)

microstate_retest %>%
  group_by(Microstate) %>%
  summarise(
    mean_rho = round(mean(SpearmanRho, na.rm = TRUE), 2),
    sd_rho = round(sd(SpearmanRho, na.rm = TRUE), 2),
    min_rho = round(min(SpearmanRho), 2),
    max_rho = round(max(SpearmanRho), 2)
  ) %>%
  # arrange(desc(mean_rho)) %>% 
  print(n = 50)

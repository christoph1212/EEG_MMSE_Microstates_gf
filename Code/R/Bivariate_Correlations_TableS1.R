# This script calculates bivariate correlations for the relevant features found
# in the PLSR.R script
#
# This script was created by Christoph Fruehlinger (November 2025)
# Last edit: November 2025

# Load data and get female subsample
microstate_datapath <- "Data/Microstates/Microstates_data_full.csv"

microstate_data <- read_csv(microstate_datapath, show_col_types = FALSE)

sample <- "Female"
microstate_data_subsample <- microstate_data %>%
  filter(Gender == tolower(sample))

# Create dataframes for each run
microstate_data_run1EO <- microstate_data_subsample %>%
  filter(
    Length == 40,
    Condition == "first_run_eyes_open",
    !is.na(Condition) & Condition != "<undefined>",
    rowSums(!is.na(across(-c(ID, gf_score, Gender, Age, Condition)))) > 0
    )

microstate_data_run2EO <- microstate_data_subsample %>%
  filter(
    Length == 40,
    Condition == "second_run_eyes_open",
    !is.na(Condition) & Condition != "<undefined>",
    rowSums(!is.na(across(-c(ID, gf_score, Gender, Age, Condition)))) > 0
  )

microstate_data_run3EO <- microstate_data_subsample %>%
  filter(
    Length == 40,
    Condition == "third_run_eyes_open",
    !is.na(Condition) & Condition != "<undefined>",
    rowSums(!is.na(across(-c(ID, gf_score, Gender, Age, Condition)))) > 0
  )

# Calculate correlations (see Table S1)
run1_coverageD <- cor.test(microstate_data_run1EO$coverage_D, 
                           microstate_data_run1EO$gf_score)
run2_coverageD <- cor.test(microstate_data_run2EO$coverage_D, 
                           microstate_data_run2EO$gf_score)
run3_coverageD <- cor.test(microstate_data_run3EO$coverage_D, 
                           microstate_data_run3EO$gf_score)

run1_transprobAD <- cor.test(microstate_data_run1EO$transition_probability_A_D, 
                             microstate_data_run1EO$gf_score)
run2_transprobAD <- cor.test(microstate_data_run2EO$transition_probability_A_D, 
                             microstate_data_run2EO$gf_score)
run3_transprobAD <- cor.test(microstate_data_run3EO$transition_probability_A_D, 
                             microstate_data_run3EO$gf_score)

run1_lifespanA <- cor.test(microstate_data_run1EO$lifespan_A, 
                           microstate_data_run1EO$gf_score)

run1_transprobAA <- cor.test(microstate_data_run1EO$transition_probability_A_A, 
                             microstate_data_run1EO$gf_score)

run2_lifespan_GFP_D <- cor.test(microstate_data_run2EO$lifespan_peaks_D, 
                                microstate_data_run2EO$gf_score)
run3_lifespan_GFP_D <- cor.test(microstate_data_run3EO$lifespan_peaks_D, 
                                microstate_data_run3EO$gf_score)

run2_transprob_GFP_CD <- cor.test(microstate_data_run2EO$transition_probability_peaks_C_D, 
                                  microstate_data_run2EO$gf_score)
run3_transprob_GFP_CD <- cor.test(microstate_data_run3EO$transition_probability_peaks_C_D, 
                                  microstate_data_run3EO$gf_score)

run2_transprob_GFP_DD <- cor.test(microstate_data_run2EO$transition_probability_peaks_D_D, 
                                  microstate_data_run2EO$gf_score)
run3_transprob_GFP_DD <- cor.test(microstate_data_run3EO$transition_probability_peaks_D_D, 
                                  microstate_data_run3EO$gf_score)

run2_coverageB <- cor.test(microstate_data_run2EO$coverage_B, 
                           microstate_data_run2EO$gf_score)

run2_transprob_GFP_DB <- cor.test(microstate_data_run2EO$transition_probability_peaks_D_B, 
                                  microstate_data_run2EO$gf_score)

run2_transprob_GFP_BD <- cor.test(microstate_data_run2EO$transition_probability_peaks_B_D, 
                                  microstate_data_run2EO$gf_score)

run2_transprob_GFP_FD <- cor.test(microstate_data_run2EO$transition_probability_peaks_F_D, 
                                  microstate_data_run2EO$gf_score)

run3_transprob_GFP_AB <- cor.test(microstate_data_run3EO$transition_probability_peaks_A_B, 
                                  microstate_data_run3EO$gf_score)

run3_transprob_GFP_AD <- cor.test(microstate_data_run3EO$transition_probability_peaks_A_D, 
                                  microstate_data_run3EO$gf_score)

run3_transprob_GFP_FA <- cor.test(microstate_data_run3EO$transition_probability_peaks_F_A, 
                                  microstate_data_run3EO$gf_score)

# Adjust p-values for run 1
run1_padj <- p.adjust(c(run1_coverageD$p.value, run1_transprobAD$p.value, 
                        run1_lifespanA$p.value, run1_transprobAA$p.value),
                      method = "holm")
round(run1_padj, 3)

# Adjust p-values for run 2
run2_padj <- p.adjust(c(run2_coverageD$p.value, run2_transprobAD$p.value, 
                        run2_lifespan_GFP_D$p.value, run2_transprob_GFP_CD$p.value, 
                        run2_transprob_GFP_DD$p.value, run2_coverageB$p.value,
                        run2_transprob_GFP_DB$p.value, run2_transprob_GFP_BD$p.value, 
                        run2_transprob_GFP_FD$p.value),
                      method = "holm")
round(run2_padj, 3)

# Adjust p-values for run 3
run3_padj <- p.adjust(c(run3_coverageD$p.value, run3_transprobAD$p.value, 
                        run3_lifespan_GFP_D$p.value, run3_transprob_GFP_CD$p.value, 
                        run3_transprob_GFP_DD$p.value, run3_transprob_GFP_AB$p.value,
                        run3_transprob_GFP_AD$p.value, run3_transprob_GFP_FA$p.value),
                      method = "holm")
round(run3_padj, 3)

# Adjust p-values for all runs (not reported)
padj <- p.adjust(c(run1_coverageD$p.value, run2_coverageD$p.value, run3_coverageD$p.value, 
                   run1_transprobAD$p.value, run2_transprobAD$p.value, run3_transprobAD$p.value, 
                   run1_lifespanA$p.value, run1_transprobAA$p.value, run2_lifespan_GFP_D$p.value, 
                   run3_lifespan_GFP_D$p.value, run2_transprob_GFP_CD$p.value, run3_transprob_GFP_CD$p.value,
                   run2_transprob_GFP_DD$p.value, run3_transprob_GFP_DD$p.value, run2_coverageB$p.value, 
                   run2_transprob_GFP_DB$p.value,run2_transprob_GFP_BD$p.value, run2_transprob_GFP_FD$p.value, 
                   run3_transprob_GFP_AB$p.value, run3_transprob_GFP_AD$p.value, run3_transprob_GFP_FA$p.value), 
                 method = "holm")
round(padj, 3)

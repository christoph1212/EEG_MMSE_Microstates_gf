####################################
####### Descriptive Analysis ####### 
####################################

## Load relevant packages for further analysis
if (!require("pacman")) install.packages("pacman")
pacman::p_load(psych, apaTables, tidyverse, effsize, ggsignif, ggpubr, psychometric, 
               stringr, ggstance, ggh4x)

options(scipen=999)

# ----------------------------------
# ---------- Demographics ----------
# ----------------------------------

## Load Participant IDs for Pre- and Post-Subs
data_mmse <- read.csv("Data/MMSEData/MMSE_data_full.csv", header = TRUE)
subs_mmse_main <- data_mmse %>% 
  filter(Condition == "first_run_eyes_open" & Length == 40) %>% 
  pull("ID") %>% 
  unique() %>%
  sort()

data_microstate <- read.csv("Data/Microstates/Microstates_data_full.csv", header = TRUE)
subs_microstate_main <- data_microstate %>% 
  filter(Condition == "first_run_eyes_open" & Length == 40) %>% 
  pull("ID") %>% 
  unique()%>%
  sort()

if (sum(subs_mmse_main == subs_microstate_main) == length(subs_mmse_main)) {
  subs = subs_microstate_main
}

## Descriptive Analysis

# Add Sociodemographic Data

socio_table = read.csv("Data/SocioDemographics.txt", header = T, na.strings = c("", "NA"))

## Convert columns to factors and rename levels
socio_table$Gender <- as.factor(socio_table$Gender)
socio_table$HormonalContraceptives <- as.factor(socio_table$HormonalContraceptives)
socio_table$Ethnicity <- as.factor(socio_table$Ethnicity)
socio_table$MaritalStatus <- as.factor(socio_table$MaritalStatus)
socio_table$Occupancy <- as.factor(socio_table$Occupancy)
socio_table$HighestDegree <- as.factor(socio_table$HighestDegree)
socio_table$monthlyBruttoIncome <- as.factor(socio_table$monthlyBruttoIncome)

levels(socio_table$Gender)[levels(socio_table$Gender) == 1] <- "Female"
levels(socio_table$Gender)[levels(socio_table$Gender) == 2] <- "Male"
levels(socio_table$HormonalContraceptives)[levels(socio_table$HormonalContraceptives) == 1] <- "Yes"
levels(socio_table$HormonalContraceptives)[levels(socio_table$HormonalContraceptives) == 2] <- "No"
levels(socio_table$Ethnicity)[levels(socio_table$Ethnicity) == 1] <- "European"
levels(socio_table$Ethnicity)[levels(socio_table$Ethnicity) == 2] <- "Arabic"
levels(socio_table$Ethnicity)[levels(socio_table$Ethnicity) == 3] <- "African"
levels(socio_table$Ethnicity)[levels(socio_table$Ethnicity) == 4] <- "Asian"
levels(socio_table$Ethnicity)[levels(socio_table$Ethnicity) == 5] <- "Hisp"
levels(socio_table$Ethnicity)[levels(socio_table$Ethnicity) == 6] <- "other"
levels(socio_table$MaritalStatus)[levels(socio_table$MaritalStatus) == 1] <- "Single"
levels(socio_table$MaritalStatus)[levels(socio_table$MaritalStatus) == 2] <- "Relationship"
levels(socio_table$MaritalStatus)[levels(socio_table$MaritalStatus) == 3] <- "legal Partnership"
levels(socio_table$MaritalStatus)[levels(socio_table$MaritalStatus) == 4] <- "Divorced"
levels(socio_table$MaritalStatus)[levels(socio_table$MaritalStatus) == 5] <- "Widowed"
levels(socio_table$Occupancy)[levels(socio_table$Occupancy) == 1] <- "Working"
levels(socio_table$Occupancy)[levels(socio_table$Occupancy) == 2] <- "Student"
levels(socio_table$Occupancy)[levels(socio_table$Occupancy) == 3] <- "unemployed"
levels(socio_table$HighestDegree)[levels(socio_table$HighestDegree) == 1] <- "Hauptschule"
levels(socio_table$HighestDegree)[levels(socio_table$HighestDegree) == 2] <- "Realschule"
levels(socio_table$HighestDegree)[levels(socio_table$HighestDegree) == 3] <- "Abitur"
levels(socio_table$HighestDegree)[levels(socio_table$HighestDegree) == 4] <- "Hochschulabschluss"
levels(socio_table$HighestDegree)[levels(socio_table$HighestDegree) == 5] <- "Promotion"
levels(socio_table$monthlyBruttoIncome)[levels(socio_table$monthlyBruttoIncome) == 1] <- "<500"
levels(socio_table$monthlyBruttoIncome)[levels(socio_table$monthlyBruttoIncome) == 2] <- "500-1000"
levels(socio_table$monthlyBruttoIncome)[levels(socio_table$monthlyBruttoIncome) == 3] <- "1000-2000"
levels(socio_table$monthlyBruttoIncome)[levels(socio_table$monthlyBruttoIncome) == 4] <- "2000-3000"
levels(socio_table$monthlyBruttoIncome)[levels(socio_table$monthlyBruttoIncome) == 5] <- ">3000"
levels(socio_table$monthlyBruttoIncome)[levels(socio_table$monthlyBruttoIncome) == 6] <- "no info"

# Filter included participants
socio_fullSample <- socio_table %>% 
  filter(ID %in% subs)

# Add missing sub
socio_fullSample <- socio_fullSample %>%
  add_row(ID = "sub-NN04AN26")

# Calculate descriptive statistics of the sample demographics
# absolute Values
table(socio_fullSample$HighestDegree, useNA = "always")
table(socio_fullSample$Occupancy, useNA = "always")
table(socio_fullSample$Ethnicity, useNA = "always")

# relative Values
prop.table(table(socio_fullSample$HighestDegree, useNA = "always"))
prop.table(table(socio_fullSample$Occupancy, useNA = "always"))
prop.table(table(socio_fullSample$Ethnicity, useNA = "always"))


# Calculate McDonald's Omega and Cronbach's Alpha 
AUC_omega <- omega(data_mmse[,11:22]) # mmse vectors - AUC
MaxSlope_omega <- omega(data_mmse[,11:14]) # MaxSlope
AvgEnt_omega <- omega(data_mmse[,19:22]) # AvgEnt

cat(paste0("\n", "MMSE McDonald's Omega:", "\n", "AUC: ", round(AUC_omega$omega.tot, 3),
           "\n", "MaxSlope: ", round(MaxSlope_omega$omega.tot, 3), "\n", 
           "AvgEnt: ", round(AvgEnt_omega$omega.tot, 3), "\n"))

cov_omega <- omega(data_microstate[, c("coverage_A", "coverage_B", "coverage_C", 
                                       "coverage_D", "coverage_F")])
freq_omega <- omega(data_microstate[, c("frequence_A", "frequence_B", "frequence_C", 
                                       "frequence_D", "frequence_F")])
ls_omega <- omega(data_microstate[, c("lifespan_A", "lifespan_B", "lifespan_C",
                                       "lifespan_D", "lifespan_F")])
ls_peak_omega <- omega(data_microstate[, c("lifespan_peaks_A", "lifespan_peaks_B", 
                                       "lifespan_peaks_C", "lifespan_peaks_D", 
                                       "lifespan_peaks_F")])
trans_omega <- omega(data_microstate[, grep("^transition_probability_[^_]+_[^_]+$",
                                          names(data_microstate), value = TRUE)])
trans_peak_omega <- omega(data_microstate[, grep("^transition_probability_peaks_",
                                          names(data_microstate), value = TRUE)])

cat(paste0("\n", "Microstate McDonald's Omega:", "\n", 
           "Coverage: ", round(cov_omega$omega.tot, 3), "\n", 
           "Frequency: ", round(freq_omega$omega.tot, 3), "\n", 
           "Lifespan: ", round(ls_omega$omega.tot, 3), "\n",
           "Lifespan (GFP): ", round(ls_peak_omega$omega.tot, 3), "\n", 
           "Trans. Prob.: ", round(trans_omega$omega.tot, 3), "\n", 
           "Trans. Prob. (GFP): ", round(trans_peak_omega$omega.tot, 3), "\n"))


IST_full <- read.csv("Data/IST_table.csv", sep = ";", header = T)

IST_full$fluid <- rowSums(IST_full[,2:21])

IST_filtered <- IST_full %>%
  filter(ID %in% subs)

gf_omega <- omega(IST_filtered[,2:21])

cat(paste0("\n", "Fluid Intelligence McDonald's Omega:", "\n",
           "gf: ", round(gf_omega$omega.tot, 3)))

if (!require("pacman")) install.packages("pacman")
pacman::p_load(pls, mvdalab, tidyverse, haven, rlist, cowplot, ggtext)


options(max.print=1000000)

rm(list = ls())

# load data set
mmse_datapath = "Data/MMSEData/MMSE_data_full.csv"
mmse_data <- read_csv(mmse_datapath)

mmse_data <- mmse_data %>%
  filter(
    Length == 40,
    Condition == "first_run_eyes_open",
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

mmse_data_wide_male <- mmse_data_wide %>%
  filter(Gender == "male")

X_AUC <- as.matrix(mmse_data_wide_male %>% select(contains("auc")))
X_AVG <- as.matrix(mmse_data_wide_male %>% select(contains("avg_entropy")))
X_MaxSlope <- as.matrix(mmse_data_wide_male %>% select(contains("max_slope")))
X <- cbind(X_AUC,X_AVG,X_MaxSlope)

y <- mmse_data_wide_male$gf_score

# cases used for training
PLSmodel <- list(y,X)
names(PLSmodel) <- list('dep','indep')

# maximum number of components for PLS models
max_comp=30

# PLSR model with CV to determine number of components
PLSR_CV <- plsr(dep ~ indep, ncomp = max_comp, data = PLSmodel, na.action = na.omit, method = "simpls", validation = "LOO", jackknife = TRUE, scale = TRUE)

# selection criterion for number of components
# IN STUDY, SELECTED NCOMP WERE SET AS FOLLOWS:
# OVERALL SAMPLE -> NCOMP = 3 (LOCAL MINIMUM CRITERION)
# WOMEN -> NCOMP = 1 (RANDOMIZATION TEST CRITERION)
# MEN -> NCOMP = 8 (LOCAL MINIMUM CRITERION)
numberSignifComp <- selectNcomp(PLSR_CV, method = "randomization", alpha=0.05, nperm=10000, plot = TRUE)

# Jack-knife estimator of regression coefficients
# jack.test(PLSR_CV,ncomp=numberSignifComp)
# jack.test(PLSR_CV,ncomp=3)

# R-squared and RMSEP for CV model
R2(PLSR_CV)
RMSEP(PLSR_CV)

# fitted (fixed-effect) PLS model with preselected number of components
PLSR_NO_CV <- plsr(dep ~ indep, ncomp = max_comp, data = PLSmodel, method = "simpls", scale = TRUE)

# R-squared and RMSEP for fitted model
R2(PLSR_NO_CV)
RMSEP(PLSR_NO_CV)

# fitted (fixed-effect) PLS model using MVDALAB function allowing for bootstrapped confidence interval estimation 
PLSR_BOOT <- plsFit(dep ~ indep, scale = TRUE, data = PLSmodel, ncomp = 1, validation = "oob", boots = 10000)

# R-squared for fitted model - the same (up to 1e^{-5}) as R2(UJ_PLSR_NO_CV) above
R2s(PLSR_BOOT)


# bootstrapped confidence interval estimation
cfs_int <- coefficients.boots(PLSR_BOOT, ncomp = 1, conf = .95)

# all intervals
x <- cfs_int[[1]]
x

# significant regression coefficients obtained from bootstrapped confidence interval estimation with alpha=0.05
filter(x,x$"2.5%">0)
filter(x,x$"97.5%"<0)

filter(x,x$"2.5%">0 & x$"97.5" < 0)

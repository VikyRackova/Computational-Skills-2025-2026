rm(list = ls())
setwd("C:/one drive/Počítač/UNI/Master studies/Computational skills")
source("EOR functions.R")
library(dplyr)
library(ggplot2)
library(knitr)
library(kableExtra)
records <- read.csv("ScanRecords.csv")
set.seed(1)
alpha <- 0.05
B <- 999
S<- 1000
#####################################################################################################################
################################################## Data Preparation #############################################
#####################################################################################################################
patient1<- records%>%
  dplyr::filter(PatientType=="Type 1")%>%  # pick the patient type 
  mutate(Duration= Duration*60)%>%        # transform the duration to minutes
  group_by(Date)%>%
  mutate(Demand = n())                    # count the number of calls in the day - proxy for the daily demand 

Patient1_Demand <- patient1%>%            # store the number of calls in a day as a unique data set to avoid repetition 
  select(Date,Demand)%>%
  distinct(Date, .keep_all = TRUE)

patient2<- records%>%
  dplyr::filter(PatientType=="Type 2")%>%  # pick the patient type 
  mutate(Duration= Duration*60)%>%        # transform the duration to minutes
  group_by(Date)%>%
  mutate(Demand = n())                    # count the number of calls in the day - proxy for the daily demand 


Patient2_Demand <- patient2%>%            # store the number of calls in a day as a unique data set to avoid repetition
  select(Date,Demand)%>%
  distinct(Date, .keep_all = TRUE)
# descriptive statistics
summary(patient1)
kable(summary(patient1),
      format = "latex",
      booktabs = TRUE,
      longtable = FALSE,
      escape = TRUE) %>%
  kable_styling(latex_options = c("hold_position"))

summary(patient2)
kable(summary(patient2),
      format = "latex",
      booktabs = TRUE,
      longtable = FALSE,
      escape = TRUE) %>%
  kable_styling(latex_options = c("hold_position"))


#####################################################################################################################
################################ Patient 1 - Testing for the distribution parameters ###############################
#####################################################################################################################
# Scan Duration
results_P1_scans <- parametric.bootstrap.normal(patient1$Duration, alpha=alpha,B=B,capacities = seq(30,35,40),percentiles = c(0.90, 0.95,0.99))
results_P1_scans$results # print the bootstrap results for the main metrics with their standard errors and confidence intervals 
#make an overleaf-ready table
kable(results_P1_scans$results,
      format = "latex",
      booktabs = TRUE,
      longtable = FALSE,
      escape = TRUE) %>%
  kable_styling(latex_options = c("hold_position"))
##############
# Expected number of patients
results_P1_demand <- parametric.bootstrap.poisson(Patient1_Demand$Demand, alpha=alpha,B=B, capacities = c(18,19,20),percentiles = c(0.90, 0.95,0.99))
results_P1_demand$results # print the bootstrap results for the main metrics with their standard errors and confidence intervals 
#make an overleaf-ready table
kable(results_P1_demand$results,
      format = "latex",
      booktabs = TRUE,
      longtable = FALSE,
      escape = TRUE) %>%
  kable_styling(latex_options = c("hold_position"))

#####################################################################################################################
################################ Patient 2 - Testing for the distribution parameters ###############################
#####################################################################################################################
# Scan Duration
results_P2_scans <-nonparametric.bootstrap(patient2$Duration,alpha=alpha,B=B,capacities = c(45,50,55,60),percentiles = c(0.90, 0.95))
results_P2_scans$results # print the bootstrap results for the main metrics with their standard errors and confidence intervals
#make an overleaf-ready table
kable(results_P2_scans$results,
      format = "latex",
      booktabs = TRUE,
      longtable = FALSE,
      escape = TRUE) %>%
  kable_styling(latex_options = c("hold_position"))

##############
# Expected number of patients
results_P2_demand <-nonparametric.bootstrap(Patient2_Demand$Demand,alpha=alpha,B=B,capacities = c(10,11),percentiles = c(0.90, 0.95))
results_P2_demand$results# print the bootstrap results for the main metrics with their standard errors and confidence intervals
#make an overleaf-ready table
kable(results_P2_demand$results,
      format = "latex",
      booktabs = TRUE,
      longtable = FALSE,
      escape = TRUE) %>%
  kable_styling(latex_options = c("hold_position"))

#####################################################################################################################
################################ Validation tests  ##################################################################
#####################################################################################################################
validation_duration_P1 <- validate.parametric.normal(patient1$Duration,S=S,B=B,alpha=alpha,capacities = c(30,35,40),percentiles = c(0.90, 0.95)) 
print(validation_duration_P1) # results 
#make an overleaf-ready table
kable(validation_duration_P1$results,
      format = "latex",
      booktabs = TRUE,
      longtable = FALSE,
      escape = TRUE) %>%
  kable_styling(latex_options = c("hold_position"))

validation_demand_P1 <- validate.parametric.poisson(Patient1_Demand$Demand,S=S,B=B,alpha=alpha,capacities = c(10, 11),percentiles = c(0.90, 0.95)) 
print(validation_demand_P1)
kable(validation_demand_P1$results,
      format = "latex",
      booktabs = TRUE,
      longtable = FALSE,
      escape = TRUE) %>%
  kable_styling(latex_options = c("hold_position"))

validation_duration_P2 <- validate.nonparametric(patient2$Duration,S=S,B=B,alpha=alpha,capacities = c(45,50,55,60),percentiles = c(0.90, 0.95)) 
print(validation_duration_P2)
kable(validation_duration_P2$results,
      format = "latex",
      booktabs = TRUE,
      longtable = FALSE,
      escape = TRUE) %>%
  kable_styling(latex_options = c("hold_position"))

validation_demand_P2 <- validate.nonparametric(Patient2_Demand$Demand,S=S,B=B,alpha=alpha,capacities = c(10),percentiles = c(0.90, 0.95)) 
print(validation_demand_P2)
kable(validation_demand_P2$results,
      format = "latex",
      booktabs = TRUE,
      longtable = FALSE,
      escape = TRUE) %>%
  kable_styling(latex_options = c("hold_position"))

# Robustness check for the nonparametric bootstrap - might be unnecessary
robustness_P2_duration<- nonparametric.bootstrap.robustness(patient2$Duration,alpha=alpha,S=S,B=B,binom_size=max(patient2$Duration, na.rm = TRUE))
print(robustness_P2_duration)
kable(robustness_P2_duration,
      format = "latex",
      booktabs = TRUE,
      longtable = FALSE,
      escape = TRUE) %>%
  kable_styling(latex_options = c("hold_position"))

robustness_P2_demand<- nonparametric.bootstrap.robustness(Patient2_Demand$Demand,alpha=alpha,S=S,B=B,binom_size=max(Patient2_Demand$Demand, na.rm = TRUE))
print(robustness_P2_demand)
kable(robustness_P2_demand,
      format = "latex",
      booktabs = TRUE,
      longtable = FALSE,
      escape = TRUE) %>%
  kable_styling(latex_options = c("hold_position"))

#####################################################################################################################
################################################## Visualization #############################################
#####################################################################################################################
# Histograms - might be unnecessary, maybe just for the appendix
############################
# Patient 1
############################
histograms(patient1$Duration,
           main = "Patient 1: Scan Duration",
           xlab = "Scan Duration (minutes)",
           bins = 60,
           x_break_by = 4)

histograms(Patient1_Demand$Demand,
           main = "Patient 1: Daily Demand",
           xlab = "Number of patients",
           bins = 20,
           x_break_by = 2)

############################
# Patient 2
############################
histograms(
  x = patient2$Duration,
  main = "Patient 2: Scan Duration",
  xlab = "Scan Duration (minutes)",
  bins = 60,
  x_break_by = 4)

histograms(
  x = Patient2_Demand$Demand,
  main = "Patient 2: Daily Demand",
  xlab = "Number of patients",
  bins = 10,
  x_break_by = 1)

# Capacity exceedance probabilities - risk that the capacity is exceeded
###########################
# Patient 1
############################
plot.exceedance.probability(
  x = Patient1_Demand$Demand,
  boot_fun = parametric.bootstrap.poisson,
  alpha = alpha,
  B = B,
  grid_n = 200,
  main = "Patient 1: Exceedance Probability of Daily Demand",
  xlab = "Number of patients threshold",
  ylab = "Exceedance Probability",
)

plot.exceedance.probability(
  x = patient1$Duration,
  boot_fun = parametric.bootstrap.normal,
  alpha = alpha,
  B = B,
  grid_n = 200,
  main = " Patient 1: Exceedance Probability of Scan Duration",
  xlab = "Scan Duration threshold (minutes)",
  ylab = "Exceedance Probability",
)
############################
# Patient 2
############################
plot.exceedance.probability(
  x = Patient2_Demand$Demand,
  boot_fun = nonparametric.bootstrap,
  alpha = alpha,
  B = B,
  grid_n = 200,
  main = "",
  xlab = "Number of patients threshold",
  ylab = "Exceedance Probability"
)

plot.exceedance.probability(
  x = patient2$Duration,
  boot_fun = nonparametric.bootstrap,
  alpha = alpha,
  B = B,
  grid_n = 200,
  main = " Patient 2: Exceedance Probability of Scan Duration",
  xlab = "Scan Duration threshold (minutes)",
  ylab = "Exceedance Probability"
)

# Quantiles - in X% of cases the capacity does not exceed a certain value
############################
# Patient 1
############################

plot.quantiles(
  x = patient1$Duration,
  boot_fun = parametric.bootstrap.normal,
  main = "Patient 1: Quantiles of Scan Duration",
  ylab = "Scan Duration (minutes)",
  B = B,
  alpha = alpha
)

plot.quantiles(
  x = Patient1_Demand$Demand,
  boot_fun = parametric.bootstrap.poisson,
  main = "Patient 1: Quantiles of Daily Demand ",
  ylab = "Number of patients",
  B = B,
  alpha = alpha
)

############################
# Patient 2
############################

plot.quantiles(
  x = patient2$Duration,
  boot_fun = nonparametric.bootstrap,
  main = "Patient 2: Quantiles of Scan Duration",
  ylab = "Duration (minutes)",
  B = B,
  alpha = alpha
)

plot.quantiles(
  x = Patient2_Demand$Demand,
  boot_fun = nonparametric.bootstrap,
  main = "Patient 2: Quantiles of Daily Demand",
  ylab = "Number of patients",
  B = B,
  alpha = alpha
)



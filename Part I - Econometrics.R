rm(list = ls())
setwd("C:/Computational skills")
source("EOR functions.R")
library(dplyr)
library(ggplot2)
records <- read.csv("ScanRecords.csv")
set.seed(1)
alpha <- 0.05
B <- 999
S<- 1000

#####################################################################################################################
################################################## Data Preparation #############################################
#####################################################################################################################
patient1<- records%>%
  dplyr::filter(PatientType=="Type 1")%>%
  mutate(Duration= Duration*60)%>%
  group_by(Date)%>%
  mutate(Demand = n())

# here I assume that the number of calls represents the expected number of patients (demand) for the day
Patient1_Demand <- patient1%>%
  select(Date,Demand)%>%
  distinct(Date, .keep_all = TRUE)

patient2<- records%>%
  dplyr::filter(PatientType=="Type 2")%>%
  mutate(Duration= Duration*60)%>%
  group_by(Date)%>%
  mutate(Demand = n())

Patient2_Demand <- patient2%>%
  select(Date,Demand)%>%
  distinct(Date, .keep_all = TRUE)

#####################################################################################################################
################################ Patient 1 - Testing for the distribution parameters ###############################
#####################################################################################################################
# Scan Duration
results_P1_scans <- parametric.bootstrap.normal(patient1$Duration, alpha,B,capacities = seq(30,35,40),percentiles = c(0.90, 0.95,0.99))
##############
# Expected number of patients
results_P1_demand <- parametric.bootstrap.poisson(Patient1_Demand$Demand, alpha,B, capacities = c(18,19,20),percentiles = c(0.90, 0.95,0.99))
#####################################################################################################################
################################ Patient 2 - Testing for the distribution parameters ###############################
#####################################################################################################################
# Scan Duration
results_P2_scans <-nonparametric.bootstrap(patient2$Duration,alpha,B,capacities = c(45,50,55,60),percentiles = c(0.90, 0.95))
##############
# Expected number of patients
results_P2_demand <-nonparametric.bootstrap(Patient2_Demand$Demand,alpha,B,capacities = c(10,11),percentiles = c(0.90, 0.95))


#####################################################################################################################
################################ Validation tests  ##################################################################
#####################################################################################################################
validation_duration_P1 <- validate.parametric.normal(patient1$Duration,S=S,B=B,alpha=alpha,capacities = c(30,35,40),percentiles = c(0.90, 0.95)) 
validation_demand_P1 <- validate.parametric.poisson(Patient1_Demand$Demand,S=S,B=B,alpha=alpha,capacities = c(10, 11),percentiles = c(0.90, 0.95)) 

validation_duration_P2 <- validate.nonparametric(patient2$Duration,S=S,B=B,alpha=alpha,capacities = c(45,50,55,60),percentiles = c(0.90, 0.95)) 
validation_demand_P2 <- validate.nonparametric(Patient2_Demand$Demand,S=S,B=B,alpha=alpha,capacities = c(10, 11),percentiles = c(0.90, 0.95)) 

# Robustness check for the nonparametric bootstrap - might be unnecessary
robustness_P2_duration<- nonparametric.bootstrap.robustness(patient2$Duration,alpha=alpha,S=S,B=B,binom_size=max(patient2$Duration))
robustness_P2_demand<- nonparametric.bootstrap.robustness(Patient2_Demand$Demand,alpha=alpha,S=S,B=B,binom_size=max(Patient2_Demand$Demand))

#####################################################################################################################
################################################## Visualization #############################################
#####################################################################################################################
# Histograms - might be unnecessary, maybe just for the appendix
############################
# Patient 1
############################
histograms(patient1$Duration,
           main = "Patient 1: Scan Duration",
           xlab = "Minutes",
           bins = 60,
           x_break_by = 4)

histograms(Patient1_Demand$Demand,
           main = "Patient 1: Daily Demand",
           xlab = "Demand (patients/day)",
           bins = 20,
           x_break_by = 2)

############################
# Patient 2
############################
histograms(
  x = patient2$Duration,
  main = "Patient 2: Scan Duration",
  xlab = "Minutes",
  bins = 60,
  x_break_by = 4)

histograms(
  x = Patient2_Demand$Demand,
  main = "Patient 2: Daily Demand",
  xlab = "Demand (patients/day)",
  bins = 10,
  x_break_by = 1)

# Capacity exceedance probabilities - risk that the capacity is exceeded
###########################
# Patient 1
############################
plot.exceedance.probability(
  x = Patient1_Demand$Demand,
  boot_fun = parametric.bootstrap.poisson,
  alpha = 0.05,
  B = 999,
  grid_n = 200,
  main = "Patient 1: Exceedance Probability of Daily Demand",
  xlab = "Daily demand threshold",
  ylab = "P(Demand > C)",
)

plot.exceedance.probability(
  x = patient1$Duration,
  boot_fun = parametric.bootstrap.normal,
  alpha = 0.05,
  B = 999,
  grid_n = 200,
  main = "Patient 1: Exceedance Probability of Scan Duration ",
  xlab = "Duration threshold",
  ylab = "P(Duration > C)",
)
############################
# Patient 2
############################
plot.exceedance.probability(
  x = Patient2_Demand$Demand,
  boot_fun = nonparametric.bootstrap,
  alpha = 0.05,
  B = 999,
  grid_n = 200,
  main = "Patient 2: Exceedance Probability of Daily Demand",
  xlab = "Demand threshold",
  ylab = "P(Duration > C)"
)

plot.exceedance.probability(
  x = patient2$Duration,
  boot_fun = nonparametric.bootstrap,
  alpha = 0.05,
  B = 999,
  grid_n = 200,
  main = "Patient 2: Exceedance Probability of Scan Duration ",
  xlab = "Duration threshold",
  ylab = "P(Duration > C)"
)
# Quantiles - in X% of cases the capacity does not exceed a certain value
############################
# Patient 1
############################

plot.quantiles(
  x = patient1$Duration,
  boot_fun = parametric.bootstrap.normal,
  main = "Patient 1: Quantiles of Scan Duration",
  ylab = "Duration (minutes)",
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




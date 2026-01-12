## ===============================
## MRI Case — Part I (R)
## Data: ScanRecords.csv
## Base R only (no packages)
## ===============================

set.seed(515)
B <- 999          # bootstrap reps
alpha <- 0.05     # for 95% CIs

## 1) Load & parse ------------------------------------------------------
df <- read.csv("/Users/klaasquak/Desktop/ComputationalResearchSkills/ScanRecords.csv", stringsAsFactors = FALSE)

# Expecting columns: Date, Time, Duration, PatientType
df$Date <- as.Date(df$Date)
df$PatientType <- factor(df$PatientType)

if (!all(c("Date","Time","Duration","PatientType") %in% names(df))) {
  stop("Unexpected columns in ScanRecords.csv")
}

## 2) Helpers -----------------------------------------------------------
perc_ci <- function(x, level = 0.95) {
  probs <- c((1 - level)/2, 1 - (1 - level)/2)
  quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 6)
}

## 3) Split -------------------------------------------------------------
d1 <- subset(df, PatientType=="Type 1")
d2 <- subset(df, PatientType=="Type 2")

## ---- TYPE 1: durations ~ Normal; arrivals/day ~ Poisson --------------
# Durations
x1 <- d1$Duration
n1 <- length(x1)
mu1_hat <- mean(x1); sd1_hat <- sd(x1)

boot_mu1   <- numeric(B)
boot_p90_1 <- numeric(B)
boot_p95_1 <- numeric(B)

for (b in 1:B) {
  xb <- rnorm(n1, mean = mu1_hat, sd = sd1_hat)
  boot_mu1[b]   <- mean(xb)
  boot_p90_1[b] <- quantile(xb, 0.90, names = FALSE, type = 6)
  boot_p95_1[b] <- quantile(xb, 0.95, names = FALSE, type = 6)
}

ci_mu1   <- perc_ci(boot_mu1, 0.95)
ci_p90_1 <- perc_ci(boot_p90_1, 0.95)
ci_p95_1 <- perc_ci(boot_p95_1, 0.95)

# Arrivals/day
counts1 <- aggregate(list(Count = d1$PatientType), list(Date = d1$Date), length)$Count
lambda1_hat <- mean(counts1)

K <- c(8,10,12)  # thresholds to report
boot_lambda1 <- numeric(B)
boot_prob_ge_k1 <- matrix(NA_real_, nrow=B, ncol=length(K),
                          dimnames=list(NULL, paste0("P(N>=",K,")")))

for (b in 1:B) {
  nb_days <- length(counts1)
  Ns <- rpois(nb_days, lambda = lambda1_hat)
  boot_lambda1[b] <- mean(Ns)
  boot_prob_ge_k1[b, ] <- colMeans(matrix(as.numeric(outer(Ns, K, `>=`)),
                                          nrow = nb_days))
}
ci_lambda1     <- perc_ci(boot_lambda1, 0.95)
ci_prob_ge_k1  <- t(apply(boot_prob_ge_k1, 2, perc_ci, level=0.95))

## ---- TYPE 2: unknown dists — nonparametric bootstrap -----------------
# Durations
x2 <- d2$Duration
n2 <- length(x2)
mu2_hat   <- mean(x2)
p90_2_hat <- quantile(x2, 0.90, names = FALSE, type = 6)
p95_2_hat <- quantile(x2, 0.95, names = FALSE, type = 6)

boot_mu2   <- numeric(B)
boot_p90_2 <- numeric(B)
boot_p95_2 <- numeric(B)

for (b in 1:B) {
  xb <- sample(x2, size = n2, replace = TRUE)
  boot_mu2[b]   <- mean(xb)
  boot_p90_2[b] <- quantile(xb, 0.90, names = FALSE, type = 6)
  boot_p95_2[b] <- quantile(xb, 0.95, names = FALSE, type = 6)
}

ci_mu2   <- perc_ci(boot_mu2, 0.95)
ci_p90_2 <- perc_ci(boot_p90_2, 0.95)
ci_p95_2 <- perc_ci(boot_p95_2, 0.95)

# Arrivals/day
counts2 <- aggregate(list(Count = d2$PatientType), list(Date = d2$Date), length)$Count
lambda2_hat <- mean(counts2)

prob_ge_k2_hat <- colMeans(matrix(as.numeric(outer(counts2, K, `>=`)),
                                  nrow = length(counts2)))

boot_lambda2 <- numeric(B)
boot_prob_ge_k2 <- matrix(NA_real_, nrow=B, ncol=length(K),
                          dimnames=list(NULL, paste0("P(N>=",K,")")))
for (b in 1:B) {
  Ns <- sample(counts2, size = length(counts2), replace = TRUE)
  boot_lambda2[b] <- mean(Ns)
  boot_prob_ge_k2[b, ] <- colMeans(matrix(as.numeric(outer(Ns, K, `>=`)),
                                          nrow = length(Ns)))
}
ci_lambda2    <- perc_ci(boot_lambda2, 0.95)
ci_prob_ge_k2 <- t(apply(boot_prob_ge_k2, 2, perc_ci, level=0.95))

## 4) Report ------------------------------------------------------------
cat("\n=== TYPE 1 (Normal durations; Poisson arrivals) ===\n")
cat(sprintf("Durations — mean (h):     %0.3f   95%% CI [%0.3f, %0.3f]\n",
            mu1_hat, ci_mu1[1], ci_mu1[2]))
cat(sprintf("Durations — 90th pct (h): %0.3f   95%% CI [%0.3f, %0.3f]\n",
            quantile(x1,0.90, names = FALSE, type = 6), ci_p90_1[1], ci_p90_1[2]))
cat(sprintf("Durations — 95th pct (h): %0.3f   95%% CI [%0.3f, %0.3f]\n",
            quantile(x1,0.95, names = FALSE, type = 6), ci_p95_1[1], ci_p95_1[2]))

cat(sprintf("\nDaily arrivals λ:         %0.2f   95%% CI [%0.2f, %0.2f]\n",
            lambda1_hat, ci_lambda1[1], ci_lambda1[2]))
for (j in seq_along(K)) {
  cat(sprintf("P(N ≥ %2d) per day:        %0.3f   95%% CI [%0.3f, %0.3f]\n",
              K[j], mean(counts1 >= K[j]), ci_prob_ge_k1[j,1], ci_prob_ge_k1[j,2]))
}

cat("\n=== TYPE 2 (Unknown distributions; nonparametric) ===\n")
cat(sprintf("Durations — mean (h):     %0.3f   95%% CI [%0.3f, %0.3f]\n",
            mu2_hat, ci_mu2[1], ci_mu2[2]))
cat(sprintf("Durations — 90th pct (h): %0.3f   95%% CI [%0.3f, %0.3f]\n",
            p90_2_hat, ci_p90_2[1], ci_p90_2[2]))
cat(sprintf("Durations — 95th pct (h): %0.3f   95%% CI [%0.3f, %0.3f]\n",
            p95_2_hat, ci_p95_2[1], ci_p95_2[2]))

cat(sprintf("\nDaily arrivals λ:         %0.2f   95%% CI [%0.2f, %0.2f]\n",
            lambda2_hat, ci_lambda2[1], ci_lambda2[2]))
for (j in seq_along(K)) {
  cat(sprintf("P(N ≥ %2d) per day:        %0.3f   95%% CI [%0.3f, %0.3f]\n",
              K[j], prob_ge_k2_hat[j], ci_prob_ge_k2[j,1], ci_prob_ge_k2[j,2]))
}

## (Optional) Monte Carlo stub — leave commented
# nr.sim <- 1000
# mc_stat <- numeric(nr.sim)
# for (i in 1:nr.sim) {
#   # simulate a “true” world and evaluate bias/coverage
# }
# mean(mc_stat); sd(mc_stat)



########################################################################################################################################
########################################################### Bootstrap Functions ########################################################
########################################################################################################################################
#######################
# Parametric 
#######################
# Poisson distribution - expected number of patients
parametric.bootstrap.poisson <- function(x, lambda=mean(x), alpha, B, capacities = c(16,17),percentiles = c(0.90, 0.95)) {
  n <- length(x)
  lambda_hat <- mean(x)
  sd_hat <- sqrt(lambda_hat)
  Q <- sqrt(n)*(lambda_hat-lambda) / sd_hat
  
  x_mean <- rep(0, times = B)
  x_median<- rep(0, times = B)
  x_sd <-rep(0, times = B)
  Q.star <- rep(0, times = B) 
  
  excessprobability<- matrix(0, nrow = B, ncol = length(capacities))
  colnames(excessprobability) <- paste0("Probability of exceeding ", capacities)
  
  quantiles <- matrix(0, nrow = B, ncol = length(percentiles))
  colnames(quantiles) <- paste0(percentiles*100, "th Quantile")
  
  for (b in 1:B) {
    x_star <- rpois(n, lambda = lambda_hat)
    
    x_mean[b] <- mean(x_star)
    x_median[b]<- median(x_star)
    x_sd[b] <-sqrt(x_mean[b])
    Q.star[b] <- sqrt(n)*(x_mean[b]-lambda_hat)/x_sd[b]
    
    
    for (j in seq_along(percentiles)) {
      P <- percentiles[j]
      quantiles[b, j] <- as.numeric(quantile(x_star, probs = P, names = FALSE))}
    
    for (j in seq_along(capacities)) {
      C <- capacities[j]
      excessprobability[b, j] <- mean(x_star > C)}
  }
  p.val <-round(mean(abs(Q.star) >= abs(Q)), 4)
  
  sample_estimates <- c(
    Mean = round(mean(x),4),
    Median = round(median(x),4),
    SD = round(sd(x),4),
    setNames(as.numeric(quantile(x, probs = percentiles)),
             paste0(percentiles*100, "th Quantile")),
    setNames(sapply(capacities, \(C) mean(x > C)),
             paste0("Probability of exceeding ", capacities))
  )
  
  bootstrap <- data.frame(
    "Bootstrap Mean"   = x_mean,
    "Bootstrap Median" = x_median,
    "Bootstrap Standard Deviation"= x_sd,
    quantiles,
    excessprobability,
    check.names = FALSE
  )
  
  estimates <- round(apply(bootstrap, 2, mean),4)
  se <- round(apply(bootstrap, 2, sd),4)
  ci <- t(apply(bootstrap, 2, quantile, probs = c(alpha/2, 1 - alpha/2)))
  ci <- round(ci, 4)
  
  
  results <- data.frame(
    Sample = round(sample_estimates,4),
    Bootstrap = estimates,
    Standard.Error = se,
    CI.lower = ci[, 1],
    CI.upper = ci[, 2],
    check.names = FALSE
  )
  invisible(list(results = results, p_value = p.val, bootstrap_draws = bootstrap))
}

# Normal distribution -  scan duration
parametric.bootstrap.normal <- function(x, mu =mean(x), alpha, B, capacities = c(16,17) ,percentiles = c(0.90, 0.95)) {
  n <- length(x)
  mu_hat <- mean(x)
  sd_hat <- sd(x)
  Q <- sqrt(n)*(mu_hat-mu) / sd_hat
  
  x_mean <- rep(0, times = B)
  x_median<- rep(0, times = B)
  x_sd <-rep(0, times = B)
  Q.star <- rep(0, times = B) 
  
  excessprobability<- matrix(0, nrow = B, ncol = length(capacities))
  colnames(excessprobability) <- paste0("Probability of exceeding ", capacities)
  
  quantiles <- matrix(0, nrow = B, ncol = length(percentiles))
  colnames(quantiles) <- paste0(percentiles*100, "th Quantile")
  
  for (b in 1:B) {
    x_star <- rnorm(n, mean = mu_hat, sd = sd_hat)
    
    x_mean[b] <- mean(x_star)
    x_median[b]<- median(x_star)
    x_sd[b] <-sd(x_star)
    Q.star[b] <- sqrt(n)*(x_mean[b]-mu_hat)/x_sd[b]
    
    for (j in seq_along(percentiles)) {
      P <- percentiles[j]
      quantiles[b, j] <- as.numeric(quantile(x_star, probs = P, names = FALSE))}
    
    for (j in seq_along(capacities)) {
      C <- capacities[j]
      excessprobability[b, j] <- mean(x_star > C)}
  }
  p.val <-round(mean(abs(Q.star) >= abs(Q)), 4)
  
  sample_estimates <- c(
    Mean = round(mean(x),4),
    Median = round(median(x),4),
    Standard.deviation = round(sd(x),4),
    Statistic = round(Q,4),
    setNames(as.numeric(quantile(x, probs = percentiles)),
             paste0(percentiles*100, "th Quantile")),
    setNames(sapply(capacities, \(C) mean(x > C)),
             paste0("Probability of exceeding ", capacities))
  )
  
  bootstrap <- data.frame(
    Mean   = x_mean,
    Median = x_median,
    Standard.deviation= x_sd,
    Statistic = Q.star,
    quantiles,
    excessprobability,
    check.names = FALSE
  )
  
  estimates <- round(apply(bootstrap, 2, mean),4)
  se <- round(apply(bootstrap, 2, sd),4)
  ci <- t(apply(bootstrap, 2, quantile, probs = c(alpha/2, 1 - alpha/2)))
  ci <- round(ci, 4)
  
  results <- data.frame(
    Sample = round(sample_estimates,4),
    Bootstrap = estimates,
    Standard.Error = se,
    CI.lower = ci[, 1],
    CI.upper = ci[, 2],
    check.names = FALSE
  )
  invisible(list(results = results, p_value = p.val, bootstrap_draws = bootstrap))
}

###################################
# Nonparametric
###################################
nonparametric.bootstrap<- function(x, mu =mean(x),alpha, B, capacities = c(16,17) ,percentiles = c(0.90, 0.95)){
  n <- length(x)
  x.bar <- mean(x)                    
  sd.hat <- sd(x)                  
  Q <- sqrt(n)*(x.bar-mu) / sd.hat
  
  x.star <- rep(0, times = B)         
  x.mean.star <- rep(0, times = B)         
  x.median.star<- rep(0, times = B)
  x.sd.star <- rep(0, times = B)         
  Q.star <- rep(0, times = B) 
  
  excessprobability<- matrix(0, nrow = B, ncol = length(capacities))
  colnames(excessprobability) <- paste0("Probability of exceeding ", capacities)
  
  quantiles <- matrix(0, nrow = B, ncol = length(percentiles))
  colnames(quantiles) <- paste0(percentiles*100, "th Quantile")
  
  for (b in 1:B) {
    J <- sample.int(n, size = n, replace = TRUE)        
    x.star <- x[J]                                      
    x.mean.star[b] <- mean(x.star)
    x.median.star[b]<- median(x.star)
    x.sd.star[b] <- sd(x.star)                           
    Q.star[b] <- sqrt(n)*(x.mean.star[b]-x.bar)/x.sd.star[b]
    
    for (j in seq_along(percentiles)) {
      P <- percentiles[j]
      quantiles[b, j] <- as.numeric(quantile(x.star, probs = P, names = FALSE))}
    
    for (j in seq_along(capacities)) {
      C <- capacities[j]
      excessprobability[b, j] <- mean(x.star > C)}
  }
  p.val <-round(mean(abs(Q.star) >= abs(Q)), 4)
  
  sample_estimates <- c(
    Mean = round(mean(x),4),
    Median = round(median(x),4),
    Standard.deviation = round(sd(x),4),
    Statistic = round(Q,4),
    setNames(as.numeric(quantile(x, probs = percentiles)),
             paste0(percentiles*100, "th Quantile")),
    setNames(sapply(capacities, \(C) mean(x > C)),
             paste0("Probability of exceeding ", capacities))
  )
  
  bootstrap <- data.frame(
    Mean   = x.mean.star,
    Median = x.median.star,
    Standard.deviation= x.sd.star,
    Statistic = Q.star,
    quantiles,
    excessprobability,
    check.names = FALSE
  )
  
  estimates <- round(apply(bootstrap, 2, mean),4)
  se <- round(apply(bootstrap, 2, sd),4)
  ci <- t(apply(bootstrap, 2, quantile, probs = c(alpha/2, 1 - alpha/2)))
  ci <- round(ci, 4)
  
  results <- data.frame(
    Sample = round(sample_estimates,4),
    Bootstrap = estimates,
    Standard.Error = se,
    CI.lower = ci[, 1],
    CI.upper = ci[, 2],
    check.names = FALSE
  )
  invisible(list(results = results, p_value = p.val, bootstrap_draws = bootstrap))
}

########################################################################################################################################
########################################################### Monte-Carlo Functions ########################################################
########################################################################################################################################
###################################
# Parametric Validation
###################################
validate.parametric.poisson <- function(x,S, B, alpha, capacities = c(15, 20, 25), percentiles = c(0.90, 0.95)){
  n = length(x)
  lambda <- mean(x)
  theta <- c(
    Mean   = lambda,
    Median = qpois(0.5, lambda),
    SD     = sqrt(lambda),
    setNames(qpois(percentiles, lambda),
             paste0(percentiles * 100, "th Quantile")),
    setNames(1 - ppois(capacities, lambda),
             paste0("Probability of exceeding ", capacities))
  )
  metric_names <- names(theta)
  K <- length(theta)
  
  sample_hat <- matrix(NA_real_, nrow = S, ncol = K, dimnames = list(NULL, metric_names))
  boot_hat   <- matrix(NA_real_, nrow = S, ncol = K, dimnames = list(NULL, metric_names))
  reject <- logical(S)
  
  for (s in 1:S) {
    y <- rpois(n, lambda = lambda)
    out <- parametric.bootstrap.poisson(
      y,
      lambda = lambda,   
      alpha = alpha,
      B = B,
      capacities = capacities,
      percentiles = percentiles
    )
    
    reject[s] <- (out$p_value <= alpha)
    bootstrap_results <- out$results
    
    # "Sample" column = plug-in estimate on y
    # "Bootstrap" column = mean of bootstrap distribution
    sample_hat[s, ] <- as.numeric(bootstrap_results[metric_names, "Sample"])
    boot_hat[s, ]   <- as.numeric(bootstrap_results[metric_names, "Bootstrap"])
  }
  
  bias_sample <- colMeans(sample_hat - matrix(theta, S, K, byrow = TRUE), na.rm = TRUE)
  rmse_sample <- sqrt(colMeans((sample_hat - matrix(theta, S, K, byrow = TRUE))^2, na.rm = TRUE))
  
  bias_boot <- colMeans(boot_hat - matrix(theta, S, K, byrow = TRUE), na.rm = TRUE)
  rmse_boot <- sqrt(colMeans((boot_hat - matrix(theta, S, K, byrow = TRUE))^2, na.rm = TRUE))
  
  rejection_rate <- mean(reject, na.rm = TRUE)
  validation_results<- data.frame(
    Metric = metric_names,
    True = as.numeric(theta),
    Bias_Sample = as.numeric(bias_sample),
    RMSE_Sample = as.numeric(rmse_sample),
    Bias_Bootstrap = as.numeric(bias_boot),
    RMSE_Bootstrap = as.numeric(rmse_boot),
    row.names = NULL,
    check.names = FALSE
  )
  
  list(results = validation_results, ERF = rejection_rate)
}



validate.parametric.normal <- function(x,S,B,alpha,capacities = c(30, 35, 40),percentiles = c(0.90, 0.95)) {
  n <- length(x)
  mu <- mean(x)
  sd_hat <- sd(x)   
  
  theta <- c(
    Mean   = mu,
    Median = mu,
    SD     = sd_hat,
    setNames(qnorm(percentiles, mean = mu, sd = sd_hat),
             paste0(percentiles * 100, "th Quantile")),
    setNames(1 - pnorm(capacities, mean = mu, sd = sd_hat),
             paste0("Probability of exceeding ", capacities))
  )
  
  metric_names <- names(theta)
  K <- length(theta)
  
  sample_hat <- matrix(NA_real_, nrow = S, ncol = K, dimnames = list(NULL, metric_names))
  boot_hat   <- matrix(NA_real_, nrow = S, ncol = K, dimnames = list(NULL, metric_names))
  reject <- logical(S)
  
  for (s in 1:S) {
    y <- rnorm(n, mean = mu, sd = sd_hat)
    
    out <- parametric.bootstrap.normal(
      y,
      mu = mu,          
      alpha = alpha,
      B = B,
      capacities = capacities,
      percentiles = percentiles
    )
    
    reject[s] <- (out$p_value <= alpha)
    bootstrap_results <- out$results
    
    sample_hat[s, ] <- as.numeric(bootstrap_results[metric_names, "Sample"])
    boot_hat[s, ]   <- as.numeric(bootstrap_results[metric_names, "Bootstrap"])
  }
  
  theta_mat <- matrix(theta, nrow = S, ncol = K, byrow = TRUE)
  
  bias_sample <- colMeans(sample_hat - theta_mat, na.rm = TRUE)
  rmse_sample <- sqrt(colMeans((sample_hat - theta_mat)^2, na.rm = TRUE))
  
  bias_boot <- colMeans(boot_hat - theta_mat, na.rm = TRUE)
  rmse_boot <- sqrt(colMeans((boot_hat - theta_mat)^2, na.rm = TRUE))
  
  rejection_rate <- mean(reject, na.rm = TRUE)
  validation_results<-data.frame(
    Metric = metric_names,
    True = as.numeric(theta),
    Bias_Sample = as.numeric(bias_sample),
    RMSE_Sample = as.numeric(rmse_sample),
    Bias_Bootstrap = as.numeric(bias_boot),
    RMSE_Bootstrap = as.numeric(rmse_boot),
    row.names = NULL,
    check.names = FALSE
  )
  list(results = validation_results, ERF = rejection_rate)
}

##########################
# Nonparametric Validation
##########################
validate.nonparametric <- function(x,S,B,alpha,capacities = c(30, 35, 40),percentiles = c(0.90, 0.95)) {
  n <- length(x)
  
  theta <- c(
    Mean   = mean(x),
    Median = median(x),
    SD     = sd(x),
    setNames(as.numeric(quantile(x, probs = percentiles, names = FALSE)),
             paste0(percentiles * 100, "th Quantile")),
    setNames(sapply(capacities, function(C) mean(x > C)),
             paste0("Probability of exceeding ", capacities))
  )
  
  metric_names <- names(theta)
  K <- length(theta)
  
  sample_hat <- matrix(NA_real_, nrow = S, ncol = K, dimnames = list(NULL, metric_names))
  boot_hat   <- matrix(NA_real_, nrow = S, ncol = K, dimnames = list(NULL, metric_names))
  
  cover <- matrix(NA_integer_, nrow = S, ncol = K, dimnames = list(NULL, metric_names))
  reject <- logical(S)
  
  for (s in 1:S) {
    y <- sample(x, size = n, replace = TRUE)
    
    out <- nonparametric.bootstrap(
      y,
      alpha = alpha,
      B = B,
      capacities = capacities,
      percentiles = percentiles
    )
    
    reject[s] <- (out$p_value <= alpha)
    res <- out$results  
    
    sample_hat[s, ] <- as.numeric(res[metric_names, "Sample"])
    boot_hat[s, ]   <- as.numeric(res[metric_names, "Bootstrap"])
    
    ci_low  <- as.numeric(res[metric_names, "CI.lower"])
    ci_high <- as.numeric(res[metric_names, "CI.upper"])
    
    cover[s, ] <- as.integer(theta >= ci_low & theta <= ci_high)
  }
  
  theta_mat <- matrix(theta, nrow = S, ncol = K, byrow = TRUE)
  
  bias_sample <- colMeans(sample_hat - theta_mat, na.rm = TRUE)
  rmse_sample <- sqrt(colMeans((sample_hat - theta_mat)^2, na.rm = TRUE))
  
  bias_boot <- colMeans(boot_hat - theta_mat, na.rm = TRUE)
  rmse_boot <- sqrt(colMeans((boot_hat - theta_mat)^2, na.rm = TRUE))
  
  coverage <- colMeans(cover, na.rm = TRUE)
  
  rejection_rate <- mean(reject, na.rm = TRUE)
  validation_results<-data.frame(
    Metric = metric_names,
    True = as.numeric(theta),
    Bias_Sample = as.numeric(bias_sample),
    RMSE_Sample = as.numeric(rmse_sample),
    Bias_Bootstrap = as.numeric(bias_boot),
    RMSE_Bootstrap = as.numeric(rmse_boot),
    Coverage = as.numeric(coverage),
    row.names = NULL,
    check.names = FALSE
  )
  list(results = validation_results, ERF = rejection_rate)
}

#############################################################
# Additional robustness check for the nonparametric resutls
#############################################################
nonparametric.bootstrap.robustness <- function(x,mu = mean(x),alpha,S,B,binom_size) {
  n <- length(x)
  m <- mean(x)
  s <- sd(x)
  
  bootstrap_t_reject <- function(x, mu0, alpha, B) {
    n <- length(x)
    xbar <- mean(x)
    sdx  <- sd(x)
    
    if (!is.finite(sdx) || sdx <= 0) return(NA)
    
    T_obs <- sqrt(n) * (xbar - mu0) / sdx
    
    T_star <- replicate(B, {
      xb <- sample(x, replace = TRUE)
      sdb <- sd(xb)
      if (!is.finite(sdb) || sdb <= 0) return(NA)
      sqrt(n) * (mean(xb) - xbar) / sdb
    })
    
    crit <- quantile(abs(T_star), 1 - alpha, na.rm = TRUE)
    abs(T_obs) > crit
  }
  
  # parameters for alternative DGPs
  p_bin <- min(max(m / binom_size, 1e-6), 1 - 1e-6)
  sdlog <- sqrt(log(1 + (s/m)^2))
  meanlog <- log(m) - 0.5 * sdlog^2
  shape_g <- (m / s)^2
  rate_g  <- m / s^2
  
  reject <- matrix(NA, S, 6)
  colnames(reject) <- c("Empirical", "Normal", "Poisson",
                        "Binomial", "Lognormal", "Gamma")
  
  for (i in 1:S) {
    
    data_list <- list(
      Empirical = sample(x, n, TRUE),
      Normal    = rnorm(n, m, s),
      Poisson   = rpois(n, m),
      Binomial  = rbinom(n, binom_size, p_bin),
      Lognormal = rlnorm(n, meanlog, sdlog),
      Gamma     = rgamma(n, shape_g, rate_g)
    )
    
    for (j in seq_along(data_list)) {
      reject[i, j] <- bootstrap_t_reject(
        data_list[[j]], mu0 = mu, alpha = alpha, B = B
      )
    }
  }
  
  data.frame(
    Distribution = colnames(reject),
    RejectionRate = colMeans(reject, na.rm = TRUE),
    TargetAlpha = alpha,
    row.names = NULL
  )
}


#####################################################################################################################
################################################## Visualization #############################################
#####################################################################################################################
histograms <- function(x,main = "Histogram",xlab = "",bins = 30,x_break_by = NULL,shade_alpha_90_95 = 0.30,shade_alpha_95 = 0.60) {
  x <- x[is.finite(x)]
  q90 <- as.numeric(quantile(x, 0.90, names = FALSE))
  q95 <- as.numeric(quantile(x, 0.95, names = FALSE))
  
  h <- hist(x, breaks = bins, plot = FALSE)
  df <- data.frame(
    xmin = h$breaks[-length(h$breaks)],
    xmax = h$breaks[-1],
    mid  = h$mids,
    density = h$density
  )
  
  df$TailRegion <- factor(
    ifelse(df$mid >= q95, "95%+ tail",
           ifelse(df$mid >= q90, "90–95% tail", "Below 90%")),
    levels = c("Below 90%", "90–95% tail", "95%+ tail")
  )
  
  p <- ggplot(df, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = density, fill = TailRegion)) +
    geom_rect(color = "white") +
    labs(title = main, x = xlab, y = "Density", fill = NULL) +
    theme_minimal() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.position = "bottom"
    ) +
    scale_fill_manual(
      values = c(
        "Below 90%"   = "grey90",
        "90–95% tail" = scales::alpha("black", shade_alpha_90_95),
        "95%+ tail"   = scales::alpha("black", shade_alpha_95)))
  
  if (!is.null(x_break_by)) {
    p <- p + scale_x_continuous(
      breaks = seq(floor(min(df$xmin)), ceiling(max(df$xmax)), by = x_break_by))}
  p
}


plot.exceedance.probability <- function(x,boot_fun,alpha,B,C_grid = NULL,grid_n = 200,
                                        main = "Exceedance curve with bootstrap band",xlab = "Threshold C",ylab = "P(X > C)") {
  x <- x[is.finite(x)]
  if (is.null(C_grid)) {
    rng <- range(x)
    C_grid <- seq(rng[1], rng[2], length.out = grid_n)
  }
  
  out <- boot_fun(
    x,
    alpha = alpha,
    B = B,
    capacities = C_grid
  )
  
  boot_df <- out$bootstrap_draws
  
  exc_cols <- grep("^Probability of exceeding ", names(boot_df), value = TRUE)
  if (length(exc_cols) == 0) stop("No exceedance columns found in bootstrap_draws.")
  
  caps <- as.numeric(sub("^Probability of exceeding\\s+", "", exc_cols))
  ord <- order(caps)
  caps <- caps[ord]
  exc_cols <- exc_cols[ord]
  
  boot_mat <- as.matrix(boot_df[, exc_cols, drop = FALSE])
  ci_low  <- apply(boot_mat, 2, quantile, probs = alpha/2, na.rm = TRUE)
  ci_high <- apply(boot_mat, 2, quantile, probs = 1 - alpha/2, na.rm = TRUE)
  
  p_hat <- vapply(caps, function(C) mean(x > C), numeric(1))
  
  plot_df <- data.frame(
    C = caps,
    p_hat = p_hat,
    ci_low = ci_low,
    ci_high = ci_high
  )
  
  ggplot(plot_df, aes(x = C)) +
    geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.20) +
    geom_line(aes(y = p_hat), linewidth = 1.1) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = main, x = xlab, y = ylab) +
    theme_minimal() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
}

plot.quantiles <- function(x,boot_fun,probs = seq(0.50, 0.99, by = 0.01),key_probs = c(0.90, 0.95),alpha,B,
                           main = "",xlab = "Quantile",ylab = "",x_zoom = c(0.75, 0.99)) {
  x <- x[is.finite(x)]
  out <- boot_fun(
    x,
    alpha = alpha,
    B = B,
    capacities = c(1),
    percentiles = probs
  )
  
  boot_df <- out$bootstrap_draws
  q_cols <- grep("th Quantile$", names(boot_df), value = TRUE)
  q_probs <- as.numeric(sub("th Quantile$", "", q_cols)) / 100
  
  keep <- q_probs %in% probs
  q_cols <- q_cols[keep]
  q_probs <- q_probs[keep]
  ord <- order(q_probs)
  q_cols <- q_cols[ord]
  q_probs <- q_probs[ord]
  
  boot_mat <- as.matrix(boot_df[, q_cols, drop = FALSE])
  
  q_hat  <- as.numeric(quantile(x, probs = q_probs, names = FALSE))
  q_low  <- apply(boot_mat, 2, quantile, probs = alpha/2, na.rm = TRUE)
  q_high <- apply(boot_mat, 2, quantile, probs = 1 - alpha/2, na.rm = TRUE)
  
  df <- data.frame(prob = q_probs, q = q_hat, low = q_low, high = q_high)
  dfz <- subset(df, prob >= x_zoom[1] & prob <= x_zoom[2])
  
  y_min <- min(dfz$low, na.rm = TRUE)
  y_max <- max(dfz$high, na.rm = TRUE)
  y_lab <- y_min - 0.06 * (y_max - y_min)
  
  ggplot(dfz, aes(x = prob)) +
    geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.12) +
    geom_line(aes(y = q), linewidth = 1.1) +
    geom_vline(xintercept = key_probs, linetype = "dashed", linewidth = 0.7, alpha = 0.45) +
    coord_cartesian(ylim = c(y_lab, y_max), clip = "off") +
    scale_x_continuous(breaks = seq(0.6, 1.0, by = 0.05)) +
    labs(title = main, x = xlab, y = ylab) +
    theme_minimal() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
}
library(metafor)

# ==========================================
# Simulation settings
# ==========================================
mu <- 2
sigma2 <- 0.5
tau2 <- 0.04

rep1 <- 10^5
# rep1 <- 10   # use this first for testing speed

set.seed(1)

sample_size <- c(10,50,100)
n_sizes <- length(sample_size)

n0 <- 5
n2 <- 10
N_total <- n0 + n2

sigma2_lnorm <- seq(0.1, 1, 0.1)

method_ids <- c(
  "shi_test",
  "balakrishnan_test",
  "proposed_log_test",
  "ideal"
)

method_labels <- c(
  "Shi test",
  "Balakrishnan test",
  "New test",
  "Ideal case"
)

n_methods <- length(method_ids)

# ==========================================
# Helper functions
# ==========================================

five <- function(x) {
  c(
    min(x),
    quantile(x, 0.25, names = FALSE),
    median(x),
    quantile(x, 0.75, names = FALSE),
    max(x)
  )
}

# --------------------------------------------------
# Normal-based estimator:
# Luo et al. mean estimator + Shi et al. SD estimator
# for S3: min, Q1, median, Q3, max
# --------------------------------------------------
est_normal_luo_shi <- function(q, n) {
  xi <- 2 * qnorm((n - 0.375) / (n + 0.25))
  eta <- 2 * qnorm((0.75 * n - 0.125) / (n + 0.25))
  
  w31 <- 2.2 / (2.2 + n^0.75)
  w32 <- 0.7 - 0.72 / n^0.55
  ws <- 1 / (1 + 0.07 * n^0.6)
  
  est_mean <- w31 * (q[1] + q[5]) / 2 +
    w32 * (q[2] + q[4]) / 2 +
    (1 - w31 - w32) * q[3]
  
  est_sd <- ws * (q[5] - q[1]) / xi +
    (1 - ws) * (q[4] - q[2]) / eta
  
  est_sd <- max(est_sd, sqrt(.Machine$double.eps))
  
  list(mean = est_mean, sd = est_sd)
}

# --------------------------------------------------
# Log-normal based estimator:
# Shi et al. log-normal estimator for S3
# --------------------------------------------------
est_lognormal_shi <- function(q, n, fallback) {
  if (any(!is.finite(q)) || any(q <= 0)) {
    return(fallback)
  }
  
  q_log <- log(q)
  
  xi <- 2 * qnorm((n - 0.375) / (n + 0.25))
  eta <- 2 * qnorm((0.75 * n - 0.125) / (n + 0.25))
  
  w31 <- 2.2 / (2.2 + n^0.75)
  w32 <- 0.7 - 0.72 / n^0.55
  ws <- 1 / (1 + 0.07 * n^0.6)
  
  sm_l_est <- w31 * (q_log[1] + q_log[5]) / 2 +
    w32 * (q_log[2] + q_log[4]) / 2 +
    (1 - w31 - w32) * q_log[3]
  
  sd_l_est <- ws * (q_log[5] - q_log[1]) / xi +
    (1 - ws) * (q_log[4] - q_log[2]) / eta
  
  v_l_est <- sd_l_est^2 / (1 + 0.28 / log(n)^2)
  
  est_mean <- exp(sm_l_est + v_l_est / 2) /
    (1 + 0.405 / n * v_l_est + 0.315 / n * v_l_est^2)
  
  est_var <- exp(2 * sm_l_est + 2 * v_l_est) /
    (1 + 1.62 / n * v_l_est + 5.04 / n * v_l_est^2) -
    exp(2 * sm_l_est + v_l_est) /
    (1 + 1.62 / n * v_l_est + 1.26 / n * v_l_est^2)
  
  if (!is.finite(est_mean) || !is.finite(est_var) || est_var <= 0) {
    return(fallback)
  }
  
  list(mean = est_mean, sd = sqrt(est_var))
}

# --------------------------------------------------
# Shi et al. skewness test for S3
# TRUE = not skewed, use normal-based estimator
# FALSE = skewed, use log-normal estimator
# --------------------------------------------------
test_shi_s3 <- function(q, n) {
  if ((q[5] - q[1]) <= 0 || (q[4] - q[2]) <= 0) {
    return(TRUE)
  }
  
  cv3 <- 3 / sqrt(n) - 40 / n^3
  
  t3 <- max(
    2.65 * log(0.6 * n) *
      abs((q[1] + q[5] - 2 * q[3]) / (q[5] - q[1])) / sqrt(n),
    abs((q[2] + q[4] - 2 * q[3]) / (q[4] - q[2]))
  )
  
  abs(t3) <= cv3
}

# --------------------------------------------------
# Balakrishnan et al. skewness test
# Implemented via estimated Pearson median skewness:
# skew = (estimated mean - median) / estimated SD
# TRUE = not skewed, use normal-based estimator
# FALSE = skewed, use log-normal estimator
# --------------------------------------------------
test_balakrishnan_s3 <- function(q, n, alpha = 0.05) {
  
  if (n<=100) {wa=1} else {wa=1-0.7*(log10(n)-2)}  
  if (n<=200) {wb=0.1} else {wb=1}
  if (n<=200) {kn=0.55} else {kn=1.65}
  
  skew_est <- kn*((q[1]-2*q[3]+q[5])*wa+(q[2]-2*q[3]+q[4])*(1-wa))/((q[5]-q[1])*wb+(q[4]-q[2])*(1-wb))
  cv=1/((15*0.05+0.6)^0.4/100*n+(26*0.05+0.8)^0.4)
  abs(skew_est) <= cv
}

# --------------------------------------------------
# Proposed log skewness test
# TRUE = not skewed, use normal-based estimator
# FALSE = skewed, use log-normal estimator
# --------------------------------------------------
test_proposed_log_s3 <- function(q, n) {
  if ((q[5] - q[4]) <= 0 || (q[2] - q[1]) <= 0) {
    return(TRUE)
  }
  
  cv3_log <- 3.23 / log(n^1.25 - 3.1) + 8.45 / (n + 1.5)
  t3_log <- log((q[5] - q[4]) / (q[2] - q[1]))
  
  abs(t3_log) <= cv3_log
}

# --------------------------------------------------
# Hybrid estimator:
# skewness test decides normal-based vs log-normal based
# --------------------------------------------------
hybrid_estimate <- function(q, n, test_name) {
  est_n <- est_normal_luo_shi(q, n)
  
  if (test_name == "shi") {
    pass <- test_shi_s3(q, n)
  } else if (test_name == "balakrishnan") {
    pass <- test_balakrishnan_s3(q, n)
  } else if (test_name == "proposed_log") {
    pass <- test_proposed_log_s3(q, n)
  } else {
    stop("Unknown test_name.")
  }
  
  if (pass || any(q <= 0)) {
    return(est_n)
  } else {
    return(est_lognormal_shi(q, n, fallback = est_n))
  }
}

fit_meta <- function(yi, vi) {
  vi <- pmax(vi, .Machine$double.eps)
  suppressWarnings(
    rma.uni(yi = yi, vi = vi, method = "DL")
  )
}

# ==========================================
# Storage
# ==========================================
mean_mat <- matrix(NA_real_, nrow = n_sizes, ncol = n_methods)
bias_mat <- matrix(NA_real_, nrow = n_sizes, ncol = n_methods)
mse_mat <- matrix(NA_real_, nrow = n_sizes, ncol = n_methods)
cp_mat <- matrix(NA_real_, nrow = n_sizes, ncol = n_methods)
ll_mat <- matrix(NA_real_, nrow = n_sizes, ncol = n_methods)

colnames(mean_mat) <- method_ids
colnames(bias_mat) <- method_ids
colnames(mse_mat) <- method_ids
colnames(cp_mat) <- method_ids
colnames(ll_mat) <- method_ids

# ==========================================
# Main simulation
# ==========================================
for (x in seq_len(n_sizes)) {
  n <- sample_size[x]
  
  beta_store <- matrix(NA_real_, nrow = rep1, ncol = n_methods)
  cover_store <- matrix(FALSE, nrow = rep1, ncol = n_methods)
  length_store <- matrix(NA_real_, nrow = rep1, ncol = n_methods)
  
  colnames(beta_store) <- method_ids
  colnames(cover_store) <- method_ids
  colnames(length_store) <- method_ids
  
  for (i in seq_len(rep1)) {
    mui <- rnorm(N_total, mu, sqrt(tau2))
    
    mu_lnorm <- log(mui[(n0 + 1):N_total]) - sigma2_lnorm / 2
    
    yi_mat <- matrix(NA_real_, nrow = N_total, ncol = n_methods)
    vi_mat <- matrix(NA_real_, nrow = N_total, ncol = n_methods)
    colnames(yi_mat) <- method_ids
    colnames(vi_mat) <- method_ids
    
    for (k in seq_len(N_total)) {
      if (k <= n0) {
        dat <- rnorm(n, mui[k], sqrt(sigma2))
      } else {
        l_idx <- k - n0
        dat <- rlnorm(n, mu_lnorm[l_idx], sqrt(sigma2_lnorm[l_idx]))
      }
      
      q <- five(dat)
      
      # 1) Shi test hybrid
      est_shi <- hybrid_estimate(q, n, "shi")
      yi_mat[k, "shi_test"] <- est_shi$mean
      vi_mat[k, "shi_test"] <- est_shi$sd^2 / n
      
      # 2) Balakrishnan test hybrid
      est_bal <- hybrid_estimate(q, n, "balakrishnan")
      yi_mat[k, "balakrishnan_test"] <- est_bal$mean
      vi_mat[k, "balakrishnan_test"] <- est_bal$sd^2 / n
      
      # 3) Proposed log test hybrid
      est_prop <- hybrid_estimate(q, n, "proposed_log")
      yi_mat[k, "proposed_log_test"] <- est_prop$mean
      vi_mat[k, "proposed_log_test"] <- est_prop$sd^2 / n
      
      # 4) Ideal case
      yi_mat[k, "ideal"] <- mean(dat)
      vi_mat[k, "ideal"] <- var(dat) / n
    }
    
    for (m in seq_len(n_methods)) {
      res <- fit_meta(yi_mat[, m], vi_mat[, m])
      
      beta_store[i, m] <- as.numeric(res$beta)
      cover_store[i, m] <- (mu <= res$ci.ub && mu >= res$ci.lb)
      length_store[i, m] <- res$ci.ub - res$ci.lb
    }
  }
  
  mean_mat[x, ] <- colMeans(beta_store, na.rm = TRUE)
  bias_mat[x, ] <- mean_mat[x, ] - mu
  
  for (m in seq_len(n_methods)) {
    mse_mat[x, m] <- mean((beta_store[, m] - mu)^2, na.rm = TRUE)
  }
  
  cp_mat[x, ] <- colMeans(cover_store, na.rm = TRUE)
  ll_mat[x, ] <- colMeans(length_store, na.rm = TRUE)
  
  cat("Finished n =", n, "\n")
}

# ==========================================
# Output summary table
# ==========================================
summary_out <- data.frame(n = sample_size)

for (m in seq_len(n_methods)) {
  id <- method_ids[m]
  summary_out[[paste0("mean_", id)]] <- mean_mat[, m]
  summary_out[[paste0("bias_", id)]] <- bias_mat[, m]
  summary_out[[paste0("mse_", id)]] <- mse_mat[, m]
  summary_out[[paste0("cp_", id)]] <- cp_mat[, m]
  summary_out[[paste0("ll_", id)]] <- ll_mat[, m]
}

print(summary_out)
write.csv(summary_out, "simulated_meta_summary.csv", row.names = FALSE)

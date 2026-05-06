library(meta)
library(metafor)
library(estmeansd)

mu <- 2
sigma2 <- 0.5
tau2 <- 0.04
rep1 <- 10^5 # 建议先改为 10 测试运行速度

set.seed(1)

# 优化函数：去掉 names 属性，提高计算速度
five <- function(x) {
  c(min(x), quantile(x, 0.25, names = FALSE), median(x), quantile(x, 0.75, names = FALSE), max(x))
}

# 偏态检验 1 (标准检验)
test_s3 <- function(x, n, cv3) {
  t3 <- max(2.65 * log(0.6 * n) * abs((x[1] + x[5] - 2 * x[3]) / (x[5] - x[1])) / sqrt(n),
            abs((x[2] + x[4] - 2 * x[3]) / (x[4] - x[2])))
  return(abs(t3) <= cv3) 
}

# 偏态检验 2 (对数检验)
test_s3_log <- function(x, cv3_log) {
  t3_log <- log((x[5] - x[4]) / (x[2] - x[1]))
  return(abs(t3_log) <= cv3_log)
}

sample_size <- seq(50,500,50)
n_sizes <- length(sample_size)

# 预分配外层循环的存储向量
r1_b3 <- r2_b3 <- r_qe_b3 <- r_mln_b3 <- r3_b3 <- r4_b3 <- r4_b3_log <- numeric(n_sizes)
cp13 <- cp23 <- cp_qe3 <- cp_mln3 <- cp33 <- cp43 <- cp43_log <- numeric(n_sizes)
ll13 <- ll23 <- ll_qe3 <- ll_mln3 <- ll33 <- ll43 <- ll43_log <- numeric(n_sizes)
mse13 <- mse23 <- mse_qe3 <- mse_mln3 <- mse33 <- mse43 <- mse43_log <- numeric(n_sizes)

sigma2_lnorm <- seq(0.1, 1, 0.1)

for (x in 1:n_sizes) {
  n <- sample_size[x]
  n0 <- 5
  n2 <- 10
  N_total <- n0 + n2 
  
  cv3 <- 3 / sqrt(n) - 40 / n^3
  cv3_log <- 3.23 / (log(n^(1.25) - 3.1)) + 8.45 / (n + 1.5)
  xi <- 2 * qnorm((n - 0.375) / (n + 0.25))
  eta <- 2 * qnorm((0.75 * n - 0.125) / (n + 0.25))
  w31 <- 2.2 / (2.2 + n^0.75)
  w32 <- 0.7 - 0.72 / n^0.55
  ws <- 1 / (1 + 0.07 * n^0.6)
  
  # 为内层循环预分配内存
  r1_beta3 <- r2_beta3 <- r_qe_beta3 <- r_mln_beta3 <- r3_beta3 <- r4_beta3 <- r4_beta3_log <- numeric(rep1)
  cover13 <- cover23 <- cover_qe3 <- cover_mln3 <- cover33 <- cover43 <- cover43_log <- logical(rep1)
  l13 <- l23 <- l_qe3 <- l_mln3 <- l33 <- l43 <- l43_log <- numeric(rep1)
  
  for (i in 1:rep1) {
    mui <- rnorm(N_total, mu, sqrt(tau2))
    mu_lnorm <- log(mui[(n0 + 1):N_total] * exp(-sigma2_lnorm / 2))
    
    y_ideal <- v_ideal <- numeric(N_total)
    y_opt1 <- v_opt1 <- numeric(N_total)
    y_opt2 <- v_opt2 <- numeric(N_total)       # McGrath Box-Cox
    y_opt_qe <- v_opt_qe <- numeric(N_total)   # McGrath QE
    y_opt_mln <- v_opt_mln <- numeric(N_total) # Cai MLN
    y_opt4 <- v_opt4 <- numeric(N_total)
    y_opt5 <- v_opt5 <- numeric(N_total)
    
    for (k in 1:N_total) {
      if (k <= n0) {
        data <- rnorm(n, mui[k], sqrt(sigma2))
      } else {
        l_idx <- k - n0
        data <- rlnorm(n, mu_lnorm[l_idx], sqrt(sigma2_lnorm[l_idx]))
      }
      
      y_ideal[k] <- mean(data)
      v_ideal[k] <- var(data) / n
      
      q <- five(data)
      q_log <- log(q)
      
      pass1 <- test_s3(q, n, cv3)
      pass2 <- test_s3_log(q, cv3_log)
      
      # 5.1 正态估计方法
      sm_n <- w31 * (q[1] + q[5]) / 2 + w32 * (q[2] + q[4]) / 2 + (1 - w31 - w32) * q[3]
      sd_n <- ws * (q[5] - q[1]) / xi + (1 - ws) * (q[4] - q[2]) / eta
      
      # 5.2 McGrath Box-Cox 估计方法 (加入 suppressWarnings 避免刷屏)
      mcg_bc <- tryCatch({
        suppressWarnings(bc.mean.sd(min.val = q[1], q1.val = q[2], med.val = q[3], q3.val = q[4], max.val = q[5], n = n))
      }, error = function(e) list(est.mean = sm_n, est.sd = sd_n)) 
      sm_m_bc <- mcg_bc$est.mean
      sd_m_bc <- mcg_bc$est.sd
      
      # 5.3 McGrath Quantile Estimation (QE) 估计方法
      mcg_qe <- tryCatch({
        suppressWarnings(qe.mean.sd(min.val = q[1], q1.val = q[2], med.val = q[3], q3.val = q[4], max.val = q[5], n = n))
      }, error = function(e) list(est.mean = sm_n, est.sd = sd_n))
      sm_m_qe <- mcg_qe$est.mean
      sd_m_qe <- mcg_qe$est.sd
      
      # 5.4 Cai's MLN (Box-Cox) 估计方法
      cai_mln <- tryCatch({
        suppressWarnings(mln.mean.sd(min.val = q[1], q1.val = q[2], med.val = q[3], q3.val = q[4], max.val = q[5], n = n))
      }, error = function(e) list(est.mean = sm_n, est.sd = sd_n))
      sm_c_mln <- cai_mln$est.mean
      sd_c_mln <- cai_mln$est.sd
      
      # 5.5 我们的对数正态估计方法
      sm_l_est <- w31 * (q_log[1] + q_log[5]) / 2 + w32 * (q_log[2] + q_log[4]) / 2 + (1 - w31 - w32) * q_log[3]
      v_l_est <- (ws * (q_log[5] - q_log[1]) / xi + (1 - ws) * (q_log[4] - q_log[2]) / eta)^2 / (1 + 0.28 / log(n)^2)
      sm_ln <- exp(sm_l_est + v_l_est / 2) / (1 + 0.405 / n * v_l_est + 0.315 / n * v_l_est^2)
      sd_ln <- sqrt(exp(2 * sm_l_est + 2 * v_l_est) / (1 + 1.62 / n * v_l_est + 5.04 / n * v_l_est^2) - 
                      exp(2 * sm_l_est + v_l_est) / (1 + 1.62 / n * v_l_est + 1.26 / n * v_l_est^2))
      
      # 6. 根据检验结果选择估计量 (Hybrid Strategies)
      y_opt1[k] <- sm_n
      v_opt1[k] <- sd_n^2 / n
      
      y_opt2[k] <- ifelse(pass1, sm_n, sm_m_bc)
      v_opt2[k] <- ifelse(pass1, sd_n^2, sd_m_bc^2) / n
      
      y_opt_qe[k] <- ifelse(pass1, sm_n, sm_m_qe)
      v_opt_qe[k] <- ifelse(pass1, sd_n^2, sd_m_qe^2) / n
      
      y_opt_mln[k] <- ifelse(pass1, sm_n, sm_c_mln)
      v_opt_mln[k] <- ifelse(pass1, sd_n^2, sd_c_mln^2) / n
      
      y_opt4[k] <- ifelse(pass1, sm_n, sm_ln)
      v_opt4[k] <- ifelse(pass1, sd_n^2, sd_ln^2) / n
      
      y_opt5[k] <- ifelse(pass2, sm_n, sm_ln)
      v_opt5[k] <- ifelse(pass2, sd_n^2, sd_ln^2) / n
    }
    
    # Meta-analyses
    result1 <- rma(y_opt1, v_opt1, method = "DL") 
    result2 <- rma(y_opt2, v_opt2, method = "DL") 
    result_qe <- rma(y_opt_qe, v_opt_qe, method = "DL")
    result_mln <- rma(y_opt_mln, v_opt_mln, method = "DL")
    result3 <- rma(y_ideal, v_ideal, method = "DL")
    result4 <- rma(y_opt4, v_opt4, method = "DL") 
    result5 <- rma(y_opt5, v_opt5, method = "DL")
    
    # 保存结果
    r1_beta3[i] <- result1$beta
    cover13[i] <- (mu <= result1$ci.ub & mu >= result1$ci.lb)
    l13[i] <- result1$ci.ub - result1$ci.lb
    
    r2_beta3[i] <- result2$beta
    cover23[i] <- (mu <= result2$ci.ub & mu >= result2$ci.lb)
    l23[i] <- result2$ci.ub - result2$ci.lb
    
    r_qe_beta3[i] <- result_qe$beta
    cover_qe3[i] <- (mu <= result_qe$ci.ub & mu >= result_qe$ci.lb)
    l_qe3[i] <- result_qe$ci.ub - result_qe$ci.lb
    
    r_mln_beta3[i] <- result_mln$beta
    cover_mln3[i] <- (mu <= result_mln$ci.ub & mu >= result_mln$ci.lb)
    l_mln3[i] <- result_mln$ci.ub - result_mln$ci.lb
    
    r3_beta3[i] <- result3$beta
    cover33[i] <- (mu <= result3$ci.ub & mu >= result3$ci.lb)
    l33[i] <- result3$ci.ub - result3$ci.lb
    
    r4_beta3[i] <- result4$beta
    cover43[i] <- (mu <= result4$ci.ub & mu >= result4$ci.lb)
    l43[i] <- result4$ci.ub - result4$ci.lb
    
    r4_beta3_log[i] <- result5$beta
    cover43_log[i] <- (mu <= result5$ci.ub & mu >= result5$ci.lb)
    l43_log[i] <- result5$ci.ub - result5$ci.lb
  }
  
  # 汇总统计
  r1_b3[x] <- mean(r1_beta3)
  r2_b3[x] <- mean(r2_beta3)
  r_qe_b3[x] <- mean(r_qe_beta3)
  r_mln_b3[x] <- mean(r_mln_beta3)
  r3_b3[x] <- mean(r3_beta3)
  r4_b3[x] <- mean(r4_beta3)
  r4_b3_log[x] <- mean(r4_beta3_log)
  
  mse13[x] <- (r1_b3[x] - mu)^2 + var(r1_beta3)
  mse23[x] <- (r2_b3[x] - mu)^2 + var(r2_beta3)
  mse_qe3[x] <- (r_qe_b3[x] - mu)^2 + var(r_qe_beta3)
  mse_mln3[x] <- (r_mln_b3[x] - mu)^2 + var(r_mln_beta3)
  mse33[x] <- (r3_b3[x] - mu)^2 + var(r3_beta3)
  mse43[x] <- (r4_b3[x] - mu)^2 + var(r4_beta3)
  mse43_log[x] <- (r4_b3_log[x] - mu)^2 + var(r4_beta3_log)
  
  ll13[x] <- mean(l13)
  ll23[x] <- mean(l23)
  ll_qe3[x] <- mean(l_qe3)
  ll_mln3[x] <- mean(l_mln3)
  ll33[x] <- mean(l33)
  ll43[x] <- mean(l43)
  ll43_log[x] <- mean(l43_log)
  
  cp13[x] <- mean(cover13)
  cp23[x] <- mean(cover23)
  cp_qe3[x] <- mean(cover_qe3)
  cp_mln3[x] <- mean(cover_mln3)
  cp33[x] <- mean(cover33)
  cp43[x] <- mean(cover43)
  cp43_log[x] <- mean(cover43_log)
}

# ==========================================
# 绘图部分 (包含 6 种方法，动态调整坐标轴)
# ==========================================
postscript("simulated_meta.eps", horizontal = FALSE, width = 10, height = 6)

par(mfrow = c(2, 2), mar = c(4, 3.5, 0, 1), oma = c(0, 0, 0, 0), mgp = c(2.5, 0.8, 0)) 

# 1. Bias 图
min_bias <- min(c(r2_b3, r_qe_b3, r_mln_b3, r3_b3, r4_b3, r4_b3_log) - mu, na.rm = TRUE)
max_bias <- max(c(r2_b3, r_qe_b3, r_mln_b3, r3_b3, r4_b3, r4_b3_log) - mu, na.rm = TRUE)
ylim_bias <- c(min(0, min_bias * 1.1), max(0, max_bias * 1.1))

plot(sample_size, r2_b3 - mu, 'b', ylim = ylim_bias, col = "purple", pch = 15, cex = 1.2,
     xlab = "", ylab = "Bias", main = "", 
     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2)
points(sample_size, r_qe_b3 - mu, 'b', col = "darkgreen", pch = 18, cex = 1.2)
points(sample_size, r_mln_b3 - mu, 'b', col = "orange", pch = 8, cex = 1.2)
points(sample_size, r3_b3 - mu, 'b', lty = 2, col = "black", pch = 3, cex = 1.2)
#points(sample_size, r4_b3 - mu, 'b', col = "blue", pch = 17, cex = 1.2)
points(sample_size, r4_b3_log - mu, 'b', col = "red", pch = 16, cex = 1.2)
abline(h = 0, col = "gray50", lty = 2) 


# 2. MSE 图
max_mse <- max(c(mse23, mse_qe3, mse_mln3, mse33, mse43, mse43_log), na.rm = TRUE)
plot(sample_size, mse23, 'b', ylim = c(0.002, max_mse * 1.05), col = "purple", pch = 15, cex = 1.2,
     xlab = "", ylab = "MSE", main = "", 
     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2)
points(sample_size, mse_qe3, 'b', col = "darkgreen", pch = 18, cex = 1.2)
points(sample_size, mse_mln3, 'b', col = "orange", pch = 8, cex = 1.2)
points(sample_size, mse33, 'b', lty = 2, col = "black", pch = 3, cex = 1.2)
#points(sample_size, mse43, 'b', col = "blue", pch = 17, cex = 1.2)
points(sample_size, mse43_log, 'b', col = "red", pch = 16, cex = 1.2)

legend("topright", 
       c('McGrath BC', 'McGrath QE', 'Cai MLN', 'Our method','Ideal case'),
       lty = c(1, 1, 1, 1, 2), 
       pch = c(15, 18, 8, 16, 3), 
       col = c('purple', 'darkgreen', 'orange', 'red', 'black'),
       lwd = 1, cex = 0.85, bg = 'gray90')

# 3. Coverage Probability 图
min_cp <- min(c(cp23, cp_qe3, cp_mln3, cp33, cp43, cp43_log), na.rm = TRUE)
ylim_cp <- c(min(0.85, min_cp * 0.98), 0.96) 

plot(sample_size, cp23, 'b', ylim = ylim_cp, col = "purple", pch = 15, cex = 1.2,
     xlab = expression(n), ylab = "Coverage probability", main = "", 
     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2)
points(sample_size, cp_qe3, 'b', col = "darkgreen", pch = 18, cex = 1.2)
points(sample_size, cp_mln3, 'b', col = "orange", pch = 8, cex = 1.2)
points(sample_size, cp33, 'b', lty = 2, col = "black", pch = 3, cex = 1.2)
#points(sample_size, cp43, 'b', col = "blue", pch = 17, cex = 1.2)
points(sample_size, cp43_log, 'b', col = "red", pch = 16, cex = 1.2)
abline(h = 0.95, col = "gray50", lty = 2) 

# 4. CI Length 图
min_ll <- min(c(ll23, ll_qe3, ll_mln3, ll33, ll43, ll43_log), na.rm = TRUE)
max_ll <- max(c(ll23, ll_qe3, ll_mln3, ll33, ll43, ll43_log), na.rm = TRUE)
plot(sample_size, ll23, 'b', ylim = c(0.2, 0.28), col = "purple", pch = 15, cex = 1.2,
     xlab = expression(n), ylab = "Length of 95% CI", main = "", 
     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2)
points(sample_size, ll_qe3, 'b', col = "darkgreen", pch = 18, cex = 1.2)
points(sample_size, ll_mln3, 'b', col = "orange", pch = 8, cex = 1.2)
points(sample_size, ll33, 'b', lty = 2, col = "black", pch = 3, cex = 1.2)
#points(sample_size, ll43, 'b', col = "blue", pch = 17, cex = 1.2)
points(sample_size, ll43_log, 'b', col = "red", pch = 16, cex = 1.2)

dev.off()
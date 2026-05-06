# 设置随机种子以保证结果可重复
set.seed(2000)

# 加载必要的包
library(metamedian)
library(metafor)

###############################################################
# 1. 数据预处理
###############################################################
dat.phq9_raw2 <- dat.phq9_raw

# 修正 0 值
a1 <- which(dat.phq9_raw2$min.g1 == 0)
a2 <- which(dat.phq9_raw2$q1.g1 == 0)
dat.phq9_raw2$min.g1[a1] <- dat.phq9_raw2$min.g1[a1] + 0.01
dat.phq9_raw2$q1.g1[a2] <- dat.phq9_raw2$q1.g1[a2] + 0.26

###############################################################
# 2. 偏度检验与亚组划分
###############################################################
description <- describe_studies(dat.phq9_raw2) 
skew_indicator <- description$skew_test_g1$reject
skew_indicator1 <- skew_indicator

# 手动调整部分研究的偏态标记
skew_indicator1[c(8, 51)] <- TRUE
skew_indicator1[c(11, 23, 26)] <- FALSE

# 划分亚组
nonskewed_subgroup1 <- dat.phq9_raw2[is.na(skew_indicator1) | !skew_indicator1, ]
skewed_subgroup1 <- dat.phq9_raw2[!is.na(skew_indicator1) & skew_indicator1, ]

###############################################################
# 3. 非偏态亚组 (Non-skewed) - 使用 Luo 方法估计
###############################################################
res_luo_nonskewed1 <- metamedian::metamean(
  data = nonskewed_subgroup1, 
  mean_method = "luo", 
  sd_method = "wan/shi_normal"
)

res_yang_nonskewed1 <- metamedian::metamean(
  data = nonskewed_subgroup1, 
  mean_method = "yang"
)

###############################################################
# 4. 偏态亚组 (Skewed) - 分别使用 QE, BC, MLN 方法估计
###############################################################
res_qe_skewed1 <- metamedian::metamean(
  data = skewed_subgroup1, 
  mean_method = "qe"
)

res_bc_skewed1 <- metamedian::metamean(
  data = skewed_subgroup1, 
  mean_method = "bc"
)

res_mln_skewed1 <- metamedian::metamean(
  data = skewed_subgroup1, 
  mean_method = "mln"
)

res_le_skewed1 <- metamedian::metamean(
  data = skewed_subgroup1, 
  mean_method = "shi_lognormal"
)

###############################################################
# 4'. 所有study - 分别使用 QE, BC, MLN 方法估计
###############################################################
res_qe_all <- metamedian::metamean(
  data = dat.phq9_raw2, 
  mean_method = "qe"
)

res_bc_all <- metamedian::metamean(
  data = dat.phq9_raw2, 
  mean_method = "bc"
)

res_mln_all <- metamedian::metamean(
  data = dat.phq9_raw2, 
  mean_method = "mln"
)

###############################################################
# 5. 总体分析与亚组数据准备 (Non-skewed[Luo] + Skewed[BC])
###############################################################

# 提取并组合数据：Non-skewed (Luo) + Skewed (BC)
dat_comb_bc <- data.frame(
  author = c(nonskewed_subgroup1$author, skewed_subgroup1$author),
  yi = c(res_luo_nonskewed1$yi, res_bc_skewed1$yi),
  vi = c(res_luo_nonskewed1$vi, res_bc_skewed1$vi),
  subgroup = c(rep("Non-skewed", length(res_luo_nonskewed1$yi)),
               rep("Skewed", length(res_bc_skewed1$yi))),
  method = c(rep("Luo", length(res_luo_nonskewed1$yi)),
             rep("McGrath-BC", length(res_bc_skewed1$yi)))
)

# 总体合并 (Overall)
res_comb_bc <- rma(yi = yi, vi = vi, data = dat_comb_bc, method = "REML")

# 亚组分别合并 (为了在森林图中画出各自亚组的菱形)
res_nonskewed <- rma(yi = yi, vi = vi, data = dat_comb_bc, subset = (subgroup == "Non-skewed"), method = "REML")
res_skewed <- rma(yi = yi, vi = vi, data = dat_comb_bc, subset = (subgroup == "Skewed"), method = "REML")



###############################################################
# 6. 绘制精美亚组森林图 (完美排版版 - 修复顶部横线重叠)
###############################################################

# ---------------- 1. 绝对行号分配 (从下到上，杜绝重叠) ----------------
k_skewed <- res_skewed$k
k_nonskewed <- res_nonskewed$k

# --- 总体 (Overall) 区域 ---
y_overall_het <- 1      # 第1行：总体异质性文字
y_overall_dia <- 2      # 第2行：总体合并菱形

# --- 偏态亚组 (Skewed) 区域 ---
y_skewed_het <- 4       # 第4行：偏态亚组异质性文字 (与总体隔开1行)
y_skewed_dia <- 5       # 第5行：偏态亚组菱形
rows_skewed <- 7:(7 + k_skewed - 1)  # 第7行开始：偏态组的具体研究
y_skewed_title <- max(rows_skewed) + 1 # 研究上方1行：偏态组标题

# --- 非偏态亚组 (Non-skewed) 区域 ---
y_nonskewed_het <- y_skewed_title + 2  # 标题上方隔1行：非偏态亚组异质性文字
y_nonskewed_dia <- y_nonskewed_het + 1 # 异质性上方1行：非偏态亚组菱形
rows_nonskewed <- (y_nonskewed_dia + 2):(y_nonskewed_dia + 2 + k_nonskewed - 1) # 非偏态组的具体研究
y_nonskewed_title <- max(rows_nonskewed) + 1 # 研究上方1行：非偏态组标题

# 【关键修改】：将顶部的留白从 +2 增加到 +3（或+4），把天花板抬高，给 header 和横线留出充足空间
y_max <- y_nonskewed_title + 3

# ---------------- 2. 准备异质性数学表达式 ----------------
# 使用 bquote 生成标准格式: Heterogeneity: τ² = 0.00, I² = 0.0%
expr_overall_het <- bquote("Heterogeneity: " ~ tau^2 == .(sprintf("%.2f", res_comb_bc$tau2)) * ", " ~ I^2 == .(sprintf("%.1f", res_comb_bc$I2)) * "%")
expr_skewed_het <- bquote("Heterogeneity: " ~ tau^2 == .(sprintf("%.2f", res_skewed$tau2)) * ", " ~ I^2 == .(sprintf("%.1f", res_skewed$I2)) * "%")
expr_nonskewed_het <- bquote("Heterogeneity: " ~ tau^2 == .(sprintf("%.2f", res_nonskewed$tau2)) * ", " ~ I^2 == .(sprintf("%.1f", res_nonskewed$I2)) * "%")

# ---------------- 3. 开始绘图 ----------------
dev.new(width = 10, height = 15) 

par(mgp = c(1.2, 0.1, 0), tcl = -0.15)
# 绘制基础森林图
sav <- forest(res_comb_bc,
              ylim = c(0, y_max),      # 使用抬高后的 y_max
              xlim = c(-26, 26),       # 控制整个画布宽度
              alim = c(-13, 20),       # 严格控制画出来的坐标轴范围
              rows = c(rows_nonskewed, rows_skewed), 
              slab = dat_comb_bc$author, 
              header = c("Study", "Effect size [95% CI]"),
              addfit = FALSE,          # 关闭默认的总体菱形，全手动添加
              xlab = "PHQ-9 score",
              cex=0.45,
              col = "darkblue", border = "darkblue")

# 获取 forest 自动计算的字体大小，保证我们手动加的字和原图一样大
cex_val <- sav$cex 

# --- 添加总体 (Overall) 菱形与异质性 ---
addpoly(res_comb_bc, row = y_overall_dia, mlab = "Overall effect", col = "darkblue", border = "darkblue", cex = cex_val)
text(sav$xlim[1], y_overall_het, pos = 4, cex = cex_val, expr_overall_het)

# --- 添加偏态亚组 (Skewed) 菱形、异质性与标题 ---
addpoly(res_skewed, row = y_skewed_dia, mlab = "Subgroup effect", col = "darkred", border = "darkred", cex = cex_val)
text(sav$xlim[1], y_skewed_het, pos = 4, cex = cex_val, expr_skewed_het)
text(sav$xlim[1], y_skewed_title, "Skewed studies", pos = 4, font = 2, cex = cex_val)

# --- 添加非偏态亚组 (Non-skewed) 菱形、异质性与标题 ---
addpoly(res_nonskewed, row = y_nonskewed_dia, mlab = "Subgroup effect", col = "darkred", border = "darkred", cex = cex_val)
text(sav$xlim[1], y_nonskewed_het, pos = 4, cex = cex_val, expr_nonskewed_het)
text(sav$xlim[1], y_nonskewed_title, "Normal studies", pos = 4, font = 2, cex = cex_val)

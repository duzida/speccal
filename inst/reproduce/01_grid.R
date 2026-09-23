## 验收 1：用 speccal 重跑论文主网格，与 4_结果/4.1_规格网格/规格曲线_v2.csv 逐项核对
## 用法：Rscript 01_grid.R [暴露列表, 逗号分隔; 默认全部 10 个]  [cores]
## 例：OMP_NUM_THREADS=1 Rscript 01_grid.R WBC,AGR 8
suppressMessages({library(survey); library(dplyr); devtools::load_all(Sys.getenv("SPECCAL_DIR", "."), quiet = TRUE)})   # 在包根目录运行
args <- commandArgs(trailingOnly = TRUE)
B <- Sys.getenv("NH_BASE", "/Users/hp/hp_data/nhanes")
source(file.path(B, "3_分析脚本/spec_core.R"))           # prep() 与论文轴定义（中文标签）
d <- prep(); options(survey.lonely.psu = "adjust")
EXS <- if (length(args) >= 1) strsplit(args[1], ",")[[1]] else EXPO
CORES <- if (length(args) >= 2) as.integer(args[2]) else 1L

## 论文的六条轴 -> speccal 轴（英文键名；核对时映射回中文标签）
ax <- spec_axes(
  outcome    = c(cdc_any = "y_cdc_any", cdc_modsev = "y_cdc_modsev", cdc_severe = "y_cdc_severe", cdc_allsite = "y_allsites",
                 efp1 = "y_efp1", efp2 = "y_efp2", efp3 = "y_efp3", efp4 = "y_efp4",
                 st4_pairs = "y_st4_pairs", st4_t20 = "y_st4_t20", st4_jaw10 = "y_st4_jaw10"),
  covariates = list(crude = character(0), demo = COV[["+人口学"]], socio = COV[["+行为社会"]], full = COV[["全模型"]]),
  coding     = c("per_sd", "log2", "q4_vs_q1"),
  missing    = list(imputed = TRUE, complete = quote(complete_case == 1)),
  sample     = list(all = TRUE, wbc_normal = quote(WBC >= 3.5 & WBC <= 11), age40 = quote(Age >= 40)),
  weighting  = c("design", "unweighted"))
print(ax)
des <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WT_TOTAL, data = d, nest = TRUE)
t0 <- Sys.time()
g <- spec_fit(ax, d, des, exposure = EXS, cores = CORES)
cat("拟合", nrow(g), "个设定，用时", round(difftime(Sys.time(), t0, units = "mins"), 1), "分钟\n")

## 映射回论文标签并核对
MAP <- list(
  outcome = c(cdc_any = "CDC/AAP 任何", cdc_modsev = "CDC/AAP 中重度", cdc_severe = "CDC/AAP 重度", cdc_allsite = "CDC 全位点变体",
              efp1 = "EFP 1mm", efp2 = "EFP 2mm", efp3 = "EFP 3mm", efp4 = "EFP 4mm",
              st4_pairs = "EFP IV期(对颌对<10)", st4_t20 = "EFP IV期(余留<20)", st4_jaw10 = "EFP IV期(单颌<10)"),
  covset  = c(crude = "粗模型", demo = "+人口学", socio = "+行为社会", full = "全模型"),
  coding  = c(per_sd = "每 SD", log2 = "每翻倍(log2)", q4_vs_q1 = "Q4 vs Q1"),
  missing = c(imputed = "MICE 插补", complete = "完整病例"),
  sample  = c(all = "全样本", wbc_normal = "WBC 3.5-11", age40 = "年龄>=40"),
  weight  = c(design = "加权", unweighted = "未加权"))
k <- as.data.frame(g) %>% transmute(exposure, outcome = MAP$outcome[outcome], covset = MAP$covset[covset],
  scale = MAP$coding[coding], sample = MAP$sample[sample], missing = MAP$missing[missing], weight = MAP$weight[weight],
  or, lo, hi, p, rd, rd_se, p_rd, rr, lrr_se, p_rr, base_risk, prevalence, n, design_df, converged, iter, max_abs_b, max_abs_se, sep_flag)
ref <- read.csv(file.path(B, "4_结果/4.1_规格网格/规格曲线_v2.csv"), fileEncoding = "UTF-8") %>% filter(exposure %in% EXS)
key <- c("exposure", "outcome", "covset", "scale", "sample", "missing", "weight")
m <- inner_join(ref, k, by = key, suffix = c("_ref", "_pkg"))
cat("参照行", nrow(ref), "| 包输出行", nrow(k), "| 匹配", nrow(m), "\n")
stopifnot(nrow(m) == nrow(ref), nrow(k) == nrow(ref))
num <- c("or", "lo", "hi", "p", "rd", "rd_se", "p_rd", "rr", "lrr_se", "p_rr", "base_risk", "prevalence", "n", "design_df", "iter", "max_abs_b", "max_abs_se")
rel <- sapply(num, function(v) { a <- m[[paste0(v, "_ref")]]; b <- m[[paste0(v, "_pkg")]]
  ok <- is.finite(a) & is.finite(b); max(abs(a[ok] - b[ok]) / pmax(abs(a[ok]), 1e-12)) })
print(signif(rel, 3))
lg <- c(converged = mean(m$converged_ref == m$converged_pkg), sep_flag = mean(m$sep_flag_ref == m$sep_flag_pkg))
print(lg)
if (all(rel < 1e-6) && all(lg == 1)) cat("\n验收通过：speccal 重现论文主网格（相对误差 <1e-6）\n") else stop("与台账不一致")

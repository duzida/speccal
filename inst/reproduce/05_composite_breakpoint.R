## 验收 5：composite_constraint / composite_coefficients / breakpoint_boot 与论文存档核对
##   约束检验 -> 4_结果/4.3_方差与约束/比值约束检验.csv（1,056 次）；成分系数 -> 成分系数_CDCany_全模型.csv；
##   拐点点估计 -> 4_结果/4.4_拐点/拐点v2_点估计.csv；Rao–Wu 200 次 -> 拐点bootstrap_RaoWu.csv（种子 20260922）；
##   naive bootstrap 前 10 次 -> 拐点bootstrap_v2.csv（种子 70000+b）
## 用法：OMP_NUM_THREADS=1 Rscript inst/reproduce/05_composite_breakpoint.R [cores]
suppressMessages({library(survey); library(dplyr); devtools::load_all(Sys.getenv("SPECCAL_DIR", "."), quiet = TRUE)})
args <- commandArgs(trailingOnly = TRUE); CORES <- if (length(args) >= 1) as.integer(args[1]) else 1L
B <- Sys.getenv("NH_BASE", "/Users/hp/hp_data/nhanes")
source(file.path(B, "3_分析脚本/spec_core.R")); d <- prep(); options(survey.lonely.psu = "adjust")
raw <- read.csv(file.path(B, "2_数据/10.NHANES/thrid_test/NHANES_0914_raw_thrid.csv"), fileEncoding = "UTF-8")
d <- d %>% mutate(NEU = neutrophils_count, LYM = lymphocytes_count, MON = monocytes_count, PLT = platelets_count,
                  NONNEU = WBC_count - neutrophils_count, ALB = Albumin, GLOB = Globulin)
des <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WT_TOTAL, data = d, nest = TRUE)
W <- list(NLR = c(NEU = 1, LYM = -1), dNLR = c(NEU = 1, NONNEU = -1), PLR = c(PLT = 1, LYM = -1), LMR = c(LYM = 1, MON = -1),
          SII = c(PLT = 1, NEU = 1, LYM = -1), SIRI = c(NEU = 1, MON = 1, LYM = -1), PIV = c(NEU = 1, PLT = 1, MON = 1, LYM = -1), AGR = c(ALB = 1, GLOB = -1))
ax <- spec_axes(
  outcome    = c(cdc_any = "y_cdc_any", cdc_modsev = "y_cdc_modsev", cdc_severe = "y_cdc_severe", cdc_allsite = "y_allsites",
                 efp1 = "y_efp1", efp2 = "y_efp2", efp3 = "y_efp3", efp4 = "y_efp4",
                 st4_pairs = "y_st4_pairs", st4_t20 = "y_st4_t20", st4_jaw10 = "y_st4_jaw10"),
  covariates = list(crude = character(0), demo = COV[["+人口学"]], socio = COV[["+行为社会"]], full = COV[["全模型"]]),
  coding     = "log2",
  sample     = list(all = TRUE, wbc_normal = quote(WBC >= 3.5 & WBC <= 11), age40 = quote(Age >= 40)))
MAP <- list(outcome = c(cdc_any = "CDC/AAP 任何", cdc_modsev = "CDC/AAP 中重度", cdc_severe = "CDC/AAP 重度", cdc_allsite = "CDC 全位点变体",
                        efp1 = "EFP 1mm", efp2 = "EFP 2mm", efp3 = "EFP 3mm", efp4 = "EFP 4mm",
                        st4_pairs = "EFP IV期(对颌对<10)", st4_t20 = "EFP IV期(余留<20)", st4_jaw10 = "EFP IV期(单颌<10)"),
            covset = c(crude = "粗模型", demo = "+人口学", socio = "+行为社会", full = "全模型"),
            sample = c(all = "全样本", wbc_normal = "WBC 3.5-11", age40 = "年龄>=40"))
## ---- 约束检验 ----
t0 <- Sys.time()
CT <- bind_rows(lapply(names(W), function(ex) as.data.frame(composite_constraint(d, des, ex, W[[ex]], ax, cores = CORES))))
cat("约束检验", nrow(CT), "次，用时", round(difftime(Sys.time(), t0, units = "mins"), 1), "分钟\n")
ref <- read.csv(file.path(B, "4_结果/4.3_方差与约束/比值约束检验.csv"), fileEncoding = "UTF-8")
k <- CT %>% transmute(exposure = composite, outcome = MAP$outcome[outcome], covset = MAP$covset[covset], sample = MAP$sample[sample],
                      k_free, p_constraint, F, p_composite, sep_flag)
m <- inner_join(ref, k, by = c("exposure", "outcome", "covset", "sample"), suffix = c("_ref", "_pkg"))
stopifnot(nrow(m) == 1056, all(m$k_free_ref == m$k_free_pkg))
e1 <- max(abs(m$F_ref - m$F_pkg) / m$F_ref, abs(m$p_constraint_ref - m$p_constraint_pkg), abs(m$p_composite_ref - m$p_composite_pkg))
cat(sprintf("约束检验 1,056 次：F / P 最大差 %.2e；sep_flag 一致 %.3f\n", e1, mean(m$sep_flag_ref == m$sep_flag_pkg))); stopifnot(e1 < 1e-8)
## ---- 成分系数 ----
CC <- bind_rows(lapply(names(W), function(ex) composite_coefficients(d, des, ex, W[[ex]], outcome = "y_cdc_any", covariates = COV[["全模型"]])))
rc <- read.csv(file.path(B, "4_结果/4.3_方差与约束/成分系数_CDCany_全模型.csv"))
rc$constituent <- sub("^l", "", rc$constituent)
m2 <- inner_join(rc, CC, by = c("composite", "constituent"), suffix = c("_ref", "_pkg")); stopifnot(nrow(m2) == 20)
e2 <- max(abs(m2$implied_ref - m2$implied_pkg), abs(m2$free_ref - m2$free_pkg), abs(m2$free_se_ref - m2$free_se_pkg), abs(m2$free_p_ref - m2$free_p_pkg))
cat(sprintf("成分系数 20 行：最大差 %.2e\n", e2)); stopifnot(e2 < 1e-8)
## ---- 拐点：点估计、Rao–Wu 200 次、naive 前 10 次 ----
PE <- read.csv(file.path(B, "4_结果/4.4_拐点/拐点v2_点估计.csv")); RW <- read.csv(file.path(B, "4_结果/4.4_拐点/拐点bootstrap_RaoWu.csv"))
V2 <- read.csv(file.path(B, "4_结果/4.4_拐点/拐点bootstrap_v2.csv"))
for (e in c("SII", "dNLR", "AGR")) {
  t0 <- Sys.time()
  bp <- breakpoint_boot(d, des, e, "y_cdc_any", COV[["全模型"]], transform = log2, grid = c(0.05, 0.95, 41), B = 200, resampling = "rao-wu",
                        seed = 20260922, cores = CORES)
  o <- PE[PE$exposure == e, ]
  stopifnot(abs(bp$observed$bp - o$bp) < 1e-9, abs(bp$observed$ddev - o$ddev) < 1e-6, abs(bp$observed$p_slope_change - o$p_seg) < 1e-8)
  r <- RW[RW$exposure == e, ]; r <- r[order(r$b), ]
  stopifnot(nrow(r) == 200, max(abs(bp$boot$bp - r$bp)) < 1e-9, all(bp$boot$at_lo == r$at_lo), all(bp$boot$at_hi == r$at_hi), max(abs(bp$boot$ddev - r$ddev)) < 1e-6)
  cat(sprintf("%-5s 点估计 %.4g、Rao–Wu 200 次与存档一致（%.1f 分钟）；区间 %.3f-%.3f，下界 %.1f%% 上界 %.1f%%\n", e, bp$observed$bp,
              as.numeric(difftime(Sys.time(), t0, units = "mins")), bp$summary$lo, bp$summary$hi, 100 * bp$summary$at_lower_bound, 100 * bp$summary$at_upper_bound))
}
## naive：bp_bootstrap2.R 用 set.seed(70000+b) 后按层重抽 PSU；这里只核对前 10 次
for (e in c("SII", "dNLR", "AGR")) {
  bn <- breakpoint_boot(d, des, e, "y_cdc_any", COV[["全模型"]], transform = log2, grid = c(0.05, 0.95, 41), B = 10, resampling = "naive",
                        seed = 70000, cores = CORES)
  r <- V2[V2$exposure == e & V2$b <= 10, ]; r <- r[order(r$b), ]
  stopifnot(nrow(r) == 10, max(abs(bn$boot$bp - r$bp)) < 1e-9, all(bn$boot$at_lo == r$at_lo), max(abs(bn$boot$ddev - r$ddev)) < 1e-4)
  cat(sprintf("%-5s naive bootstrap 前 10 次与存档一致\n", e))
}
cat("\n验收通过：composite_constraint / composite_coefficients / breakpoint_boot 重现存档\n")

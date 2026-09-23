## 验收 3：spec_defensible / spec_keep / spec_summary 与论文台账核对
##   存档的 300 次 FL 零分布（FL零分布_合并.rds）转成 spec_null 对象，不重跑；观察网格由 01 的同一定义重算（1.7 分钟）。
##   核对：5_图表/Fig2_source_data.csv（各水平假阳性率与区间）、4_结果/4.2_置换与零分布/合理子集_FL推断_300.csv（显著比例、置换 P）、
##         4_结果/4.2_置换与零分布/有效独立设定数.csv（n_eff 与区间）、results_ledger.json（share_by_exposure, Janus）
## 用法：OMP_NUM_THREADS=1 Rscript inst/reproduce/03_defensible_summary.R [cores]
suppressMessages({library(survey); library(dplyr); devtools::load_all(Sys.getenv("SPECCAL_DIR", "."), quiet = TRUE)})
args <- commandArgs(trailingOnly = TRUE); CORES <- if (length(args) >= 1) as.integer(args[1]) else 1L
B <- Sys.getenv("NH_BASE", "/Users/hp/hp_data/nhanes")
source(file.path(B, "3_分析脚本/spec_core.R")); d <- prep(); options(survey.lonely.psu = "adjust")
ax <- spec_axes(
  outcome    = c(cdc_any = "y_cdc_any", cdc_modsev = "y_cdc_modsev", cdc_severe = "y_cdc_severe", cdc_allsite = "y_allsites",
                 efp1 = "y_efp1", efp2 = "y_efp2", efp3 = "y_efp3", efp4 = "y_efp4",
                 st4_pairs = "y_st4_pairs", st4_t20 = "y_st4_t20", st4_jaw10 = "y_st4_jaw10"),
  covariates = list(crude = character(0), demo = COV[["+人口学"]], socio = COV[["+行为社会"]], full = COV[["全模型"]]),
  coding     = c("per_sd", "log2", "q4_vs_q1"),
  missing    = list(imputed = TRUE, complete = quote(complete_case == 1)),
  sample     = list(all = TRUE, wbc_normal = quote(WBC >= 3.5 & WBC <= 11), age40 = quote(Age >= 40)),
  weighting  = c("design", "unweighted"))
des <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WT_TOTAL, data = d, nest = TRUE)
g <- spec_fit(ax, d, des, exposure = EXPO, cores = CORES, quiet = TRUE)

## 存档零分布 -> spec_null（标签映射回英文键；子网格 = 插补+加权）
INV <- list(
  outcome = c("CDC/AAP 任何" = "cdc_any", "CDC/AAP 中重度" = "cdc_modsev", "CDC/AAP 重度" = "cdc_severe", "CDC 全位点变体" = "cdc_allsite",
              "EFP 1mm" = "efp1", "EFP 2mm" = "efp2", "EFP 3mm" = "efp3", "EFP 4mm" = "efp4",
              "EFP IV期(对颌对<10)" = "st4_pairs", "EFP IV期(余留<20)" = "st4_t20", "EFP IV期(单颌<10)" = "st4_jaw10"),
  covset = c("粗模型" = "crude", "+人口学" = "demo", "+行为社会" = "socio", "全模型" = "full"),
  scale  = c("每 SD" = "per_sd", "每翻倍(log2)" = "log2", "Q4 vs Q1" = "q4_vs_q1"),
  sample = c("全样本" = "all", "WBC 3.5-11" = "wbc_normal", "年龄>=40" = "age40"))
N <- readRDS(file.path(B, "4_结果/4.2_置换与零分布/FL零分布_合并.rds"))
perm <- data.frame(rep = N$rep, exposure = N$exposure, outcome = unname(INV$outcome[N$outcome]), covset = unname(INV$covset[N$covset]),
                   coding = unname(INV$scale[N$scale]), missing = "imputed", sample = unname(INV$sample[N$sample]), weight = "design",
                   estimate = log(N$or), p = N$p, rr = N$rr, p_rr = N$p_rr, sep_flag = N$sep_flag, stringsAsFactors = FALSE)
stopifnot(!anyNA(perm$outcome), !anyNA(perm$covset), !anyNA(perm$coding), !anyNA(perm$sample))
nul <- structure(list(perm = perm, scheme = "fl", reps = sort(unique(perm$rep)), seed = 500000, stratum = "SDMVSTRA",
                      restrict = list(missing = "imputed", weight = "design"), axis_cols = axis_cols(g), exposures = EXPO), class = "spec_null")
print(nul)

## ---- spec_defensible vs Fig2_source_data.csv ----
def <- spec_defensible(nul, threshold = 0.10, sensitivity = c(0.075, 0.15), within = list(covset = c("socio", "full")))
print(def)
F2 <- read.csv(file.path(B, "5_图表/Fig2_source_data.csv"), fileEncoding = "UTF-8")
F2$level <- sub("  \\(.*$", "", F2$lab)
LAB <- c(Crude = "crude", Demographics = "demo", "+ Socio-behavioural" = "socio", Full = "full",
         "CDC/AAP any" = "cdc_any", "CDC/AAP mod-severe" = "cdc_modsev", "CDC/AAP severe" = "cdc_severe", "CDC/AAP all-site" = "cdc_allsite",
         "EFP/AAP 1 mm" = "efp1", "EFP/AAP 2 mm" = "efp2", "EFP/AAP 3 mm" = "efp3", "EFP/AAP 4 mm" = "efp4",
         "Stage IV, pairs < 10" = "st4_pairs", "Stage IV, teeth < 20" = "st4_t20", "Stage IV, jaw < 10" = "st4_jaw10",
         "Per SD" = "per_sd", "Per doubling" = "log2", "Q4 vs Q1" = "q4_vs_q1",
         "Full sample" = "all", "WBC 3.5-11" = "wbc_normal", "Age 40 y or older" = "age40")
F2$key <- unname(LAB[F2$level]); stopifnot(!anyNA(F2$key))
D <- as.data.frame(def)
cmp <- merge(F2[, c("key", "FPR", "lo", "hi")], D[, c("level", "fpr", "lo", "hi")], by.x = "key", by.y = "level", suffixes = c("_ref", "_pkg"))
stopifnot(nrow(cmp) == nrow(F2))
e2 <- max(abs(cmp$FPR - 100 * cmp$fpr), abs(cmp$lo_ref - 100 * cmp$lo_pkg), abs(cmp$hi_ref - 100 * cmp$hi_pkg))
cat(sprintf("假阳性率与区间：%d 个水平，最大绝对差 %.2e 个百分点\n", nrow(cmp), e2)); stopifnot(e2 < 1e-8)
keep <- spec_keep(def, g)
stopifnot(sum(keep) == 162 * 10, all(g$covset[keep] %in% c("socio", "full")), !any(g$outcome[keep] %in% c("efp1", "efp2")))
excl <- D$level[!D$defensible]; cat("剔除:", paste(excl, collapse = ", "), "\n"); stopifnot(setequal(excl, c("crude", "demo", "efp1", "efp2")))

## ---- spec_summary（合理子集）vs 合理子集_FL推断_300.csv / 有效独立设定数.csv ----
sm <- spec_summary(g, nul, keep = keep); print(sm)
ref <- read.csv(file.path(B, "4_结果/4.2_置换与零分布/合理子集_FL推断_300.csv"), fileEncoding = "UTF-8")
m <- merge(sm, ref, by = "exposure")
e3 <- max(abs(100 * m$sig_modal - m$sig), abs(m$p_perm - m$FL_P), abs(m$direction - m$dom))
cat(sprintf("合理子集显著比例 / 置换 P / 模态方向：最大差 %.2e\n", e3)); stopifnot(e3 < 1e-8)
NE <- read.csv(file.path(B, "4_结果/4.2_置换与零分布/有效独立设定数.csv"), fileEncoding = "UTF-8")
m2 <- merge(sm, NE[NE$set == "defensible (162)", ], by = "exposure")
e4 <- max(abs(m2$n_eff_var - m2$n_eff_variance), abs(m2$n_eff_liji - m2$n_eff_LiJi), abs(100 * m2$ci_lo - m2$ci_lo.y), abs(100 * m2$ci_hi - m2$ci_hi.y),
          abs(100 * m2$null_mean_sig - m2$null_mean_sig.y), abs(100 * m2$null_sd_sig - m2$null_sd_share))
cat(sprintf("有效独立设定数（方差法、Li-Ji）与近似区间：最大差 %.2e\n", e4)); stopifnot(e4 < 1e-6)
## 全子网格 396
sub <- g$missing == "imputed" & g$weight == "design"
sm2 <- spec_summary(g, nul, keep = sub)
m3 <- merge(sm2, NE[NE$set == "all (396)", ], by = "exposure")
e5 <- max(abs(m3$n_eff_var - m3$n_eff_variance), abs(m3$n_eff_liji - m3$n_eff_LiJi), abs(100 * m3$sig_modal - m3$obs_share))
cat(sprintf("全子网格 396：n_eff 与观察比例最大差 %.2e\n", e5)); stopifnot(e5 < 1e-6)
## 全网格 1,584：与台账 share_by_exposure 与效应振动
L <- jsonlite::fromJSON(file.path(B, "4_结果/results_ledger.json"))
sm3 <- spec_summary(g); m4 <- merge(sm3, L$share_by_exposure, by = "exposure")
e6 <- max(abs(100 * m4$sig.x - m4$sig.y), abs(100 * m4$sig_rr.x - m4$sig_rr.y), abs(100 * m4$sig_rd.x - m4$sig_rd.y), abs(100 * m4$modal_share - m4$dir_all))
cat(sprintf("全网格显著比例（OR/RR/RD）与模态侧比例：最大差 %.2e（台账保留 4 位小数）\n", e6)); stopifnot(e6 < 1e-3)
VE <- read.csv(file.path(B, "4_结果/4.3_方差与约束/效应振动指标.csv"), fileEncoding = "UTF-8")
m5 <- merge(sm3, VE[VE$grid == "全网格", ], by = "exposure")
e7 <- max(abs(m5$relative_or - m5$ROR), abs(m5$relative_p - m5$RP)); cat(sprintf("相对 OR / 相对 P：最大差 %.2e\n", e7)); stopifnot(e7 < 1e-6)
stopifnot(all(m5$janus == (m5$Janus == "是")))
cat("\n验收通过：spec_defensible / spec_keep / spec_summary 重现台账\n")

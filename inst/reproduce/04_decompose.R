## 验收 4：spec_decompose 与论文台账核对
##   主效应三尺度份额与均方比 -> results_ledger.json$vd；两两交互 -> 4_结果/4.3_方差与约束/两两交互.csv；
##   等基数敏感性（定义轴 11 选 4，330 种）-> 等基数敏感性.csv；零分布校准 -> 与脚本内直接按 FL零分布_合并 逐次 anova 的结果对照（自洽）
## 用法：OMP_NUM_THREADS=1 Rscript inst/reproduce/04_decompose.R [cores]
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
L <- jsonlite::fromJSON(file.path(B, "4_结果/results_ledger.json"))

## ---- 主效应，三尺度 ----
for (k in c(or = "estimate", rr = "lrr", rd = "rd")) {
  v <- spec_decompose(g, scale = k)
  ref <- L$vd[[names(which(c(or = "estimate", rr = "lrr", rd = "rd") == k))]]
  m <- merge(v$main, ref, by = "exposure")
  ms <- sapply(m$exposure, function(e) { s <- v$shares[v$shares$exposure == e, ]; s$mean_sq[s$term == "outcome"] / s$mean_sq[s$term == "covset"] })
  err <- max(abs(100 * m$outcome - m$def), abs(100 * m$covset - m$cov), abs(ms[m$exposure] - m$ms_ratio) / m$ms_ratio)
  cat(sprintf("%-8s 定义/协变量份额与均方比：最大差 %.2e（台账 4 位小数）\n", k, err)); stopifnot(err < 1e-3)
}
## ---- 两两交互 ----
vi <- spec_decompose(g, interactions = TRUE)
INT <- read.csv(file.path(B, "4_结果/4.3_方差与约束/两两交互.csv"), fileEncoding = "UTF-8", check.names = FALSE)
canon <- function(t) sapply(strsplit(gsub("coding", "scale", t), ":", fixed = TRUE), function(v) paste(sort(v), collapse = ":"))   # 交互项按因子名排序，消除公式顺序差异
names(INT) <- c("exposure", canon(names(INT)[-1]))
for (e in EXPO) {
  s <- vi$shares[vi$shares$exposure == e, ]; s$term2 <- canon(s$term)
  r <- INT[INT$exposure == e, ]
  common <- intersect(s$term2, names(r)); stopifnot(length(common) == 22)
  err <- max(abs(100 * s$share[match(common, s$term2)] - as.numeric(r[1, common])))
  stopifnot(err < 1e-8)
}
cat("两两交互（10 暴露 × 22 项）：与存档一致\n")
## ---- 等基数敏感性：定义轴 11 选 4 ----
ve <- spec_decompose(g, equal_cardinality = list(axis = "outcome", k = 4))
EC <- read.csv(file.path(B, "4_结果/4.3_方差与约束/等基数敏感性.csv"), fileEncoding = "UTF-8") %>% filter(pool == "all 11 definitions")
for (e in EXPO) {
  sh <- ve$equal_cardinality_draws[[e]]; r <- EC[EC$exposure == e, ]
  err <- max(abs(100 * median(sh[, "outcome"]) - r$def_median), abs(100 * median(sh[, "covset"]) - r$cov_median),
             abs(100 * mean(sh[, "outcome"] > sh[, "covset"]) - r$pct_def_gt_cov))
  stopifnot(nrow(sh) == 330, err < 1e-8)
}
cat("等基数敏感性（330 种组合）：定义/协变量份额中位数与占优比例与存档一致\n")
## ---- 零分布校准（自洽核对）----
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
                   estimate = log(N$or), p = N$p, stringsAsFactors = FALSE)
nul <- structure(list(perm = perm, scheme = "fl", reps = 1:300, seed = 500000, stratum = "SDMVSTRA",
                      restrict = list(missing = "imputed", weight = "design"), axis_cols = axis_cols(g), exposures = EXPO), class = "spec_null")
sub <- g$missing == "imputed" & g$weight == "design"
vn <- spec_decompose(g, nul, keep = sub)
chk <- N %>% filter(exposure == "dNLR") %>% mutate(lo = log(or)) %>% group_by(rep) %>%
  group_modify(~{ a <- anova(lm(lo ~ outcome + covset + scale + sample, data = .x)); ss <- a[["Sum Sq"]]; tibble(def = ss[1] / sum(ss)) }) %>% ungroup()
q <- vn$null[vn$null$exposure == "dNLR" & vn$null$axis == "outcome", ]
stopifnot(abs(q$null_median - median(chk$def)) < 1e-12, abs(q$null_q95 - quantile(chk$def, .95)) < 1e-12)
cat(sprintf("零分布校准（dNLR 定义份额）：观察 %.1f%%，零分布中位 %.1f%%、95 分位 %.1f%%，P = %.3f —— 与直接计算一致\n",
            100 * q$observed, 100 * q$null_median, 100 * q$null_q95, q$p_perm))
cat("\n验收通过：spec_decompose 重现台账（主效应三尺度、两两交互、等基数敏感性、零分布校准）\n")
